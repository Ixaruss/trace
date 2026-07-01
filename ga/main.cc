#include <cstddef>
#include <iomanip>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <functional>
#include <string>
#include <unordered_set>
#include <cstring>
#include <chrono>
#include <numeric>
// Only include common.h for Global, Bin, Config, and FeasibleGrid structures.
// Local coverage logic has been moved to radar.h — do not re-embed it here.
#include "../common.h"
#include "h5.h"
#include "radar.h"
#include "config.h"

using namespace std;

// ============================================================
// CONSTANTS / DEFAULTS
// (all overridden at runtime when a config file is present)
// ============================================================
const int NUM_CHROMOSOMES = 100;
static constexpr double EARTH_R     = 6371000.0;
static constexpr double K_ATMOS     = 4.0 / 3.0;
static constexpr double K_EFF_R     = K_ATMOS * EARTH_R;
static constexpr double INV_2_KEFFR = 1.0 / (2.0 * K_EFF_R);
static constexpr double DEG2RAD     = 3.14159265358979 / 180.0;
static constexpr double RAD2DEG     = 180.0 / 3.14159265358979;

// Default fitness weights — used only when no config is loaded
static constexpr double DEFAULT_W1 = 0.9;   // unique coverage reward
static constexpr double DEFAULT_W2 = 0.75;  // low-overlap reward
static constexpr double DEFAULT_W3 = 0.9;   // high-overlap penalty
static constexpr double DEFAULT_W4 = 0.9;   // spacing-gap penalty
static constexpr double DEFAULT_X  = 0.0;   // lambda (radar-count cost, run_auto)

// ============================================================
// ENUMS
// ============================================================
enum class SelectionMethod     { TOURNAMENT, ROULETTE_WHEEL, RANK_BASED, SUS };
enum class CrossoverMethod     { SINGLE_POINT, MULTI_POINT, UNIFORM };
enum class MutationMethod      { SINGLE_POINT, SPATIAL_NUDGE, GAUSSIAN_CREEP };
enum class AutoCrossoverMethod { RANDOM_SUBSET, LENGTH_NORMALIZED, POSITIONAL_PADDING };

// ============================================================
// FLEET CONFIG STRUCTS
// ============================================================

struct RadarConfig {
    double rangeKm;
    int    count;
};

struct FlexRadarConfig {
    double rangeKm;
    int    minCount;
    int    maxCount;
};

// Struct representing a configuration case for the hyperparameter search
struct GAConfig {
    std::string name;
    int popSize;
    int generations;
    double mutationRate;
    double crossoverRate;
    SelectionMethod selection;
    CrossoverMethod crossover;
    MutationMethod mutation;
};

// Struct to keep track of benchmark results for the final leaderboard
struct BenchmarkResult {
    std::string name;
    std::string selection;
    std::string crossover;
    std::string mutation;
    int popSize;
    int generations;
    double mutationRate;
    double crossoverRate;
    double bestDistance;
    double totalTimeMs;
    double evalTimeMs;
    double evoTimeMs;
};

// ============================================================
// CHROMOSOME TYPES
// ============================================================

struct Chromosome {
    std::vector<GridPoint> radar;
    double fitness = 0.0;

    explicit Chromosome(int length) : radar(length), fitness(0.0) {}
    Chromosome() : fitness(0.0) {}
};

struct RadarGene {
    int r, c;
    int typeIndex;   // -1 = NULL sentinel for positional-padding crossover
};

struct AutoChromosome {
    std::vector<RadarGene> radars;
    double fitness = 0.0;

    int countType(int t) const {
        return (int)std::count_if(radars.begin(), radars.end(),
            [t](const RadarGene& g) { return g.typeIndex == t; });
    }
    int totalCount() const { return (int)radars.size(); }
};

// ============================================================
// THREAD POOL
// ============================================================
class ThreadPool {
private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;
    std::mutex queueMutex;
    std::condition_variable cv;
    bool stop;

public:
    ThreadPool() : stop(false) {
        size_t threads = std::thread::hardware_concurrency();
        if (threads == 0) threads = 4;

        for (size_t i = 0; i < threads; ++i) {
            workers.emplace_back([this]() {
                while (true) {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(this->queueMutex);
                        this->cv.wait(lock, [this]() {
                            return this->stop || !this->tasks.empty();
                        });
                        if (this->stop && this->tasks.empty()) return;
                        task = std::move(this->tasks.front());
                        this->tasks.pop();
                    }
                    task();
                }
            });
        }
    }

    template<class F, class... Args>
    auto enqueue(F&& f, Args&&... args)
        -> std::future<typename std::invoke_result<F, Args...>::type>
    {
        using return_type = typename std::invoke_result<F, Args...>::type;
        auto task = std::make_shared<std::packaged_task<return_type()>>(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...)
        );
        std::future<return_type> res = task->get_future();
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            if (stop) throw std::runtime_error("Enqueue requested on stopped ThreadPool context");
            tasks.emplace([task]() { (*task)(); });
        }
        cv.notify_one();
        return res;
    }

    ~ThreadPool() {
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            stop = true;
        }
        cv.notify_all();
        for (std::thread& worker : workers)
            if (worker.joinable()) worker.join();
    }
};

// ============================================================
// GENETIC ALGORITHM CLASS
// ============================================================
class GaEngine {
    Global g;
    vector<Chromosome> population;
    ThreadPool& pool;
    GaConf conf;
    bool confLoaded = false;
    double totalFitnessEvalTimeMs = 0.0;
    double totalEvolutionaryTimeMs = 0.0;

    void reloadConfig() {
        if (GaConf::is_available()) {
            conf = GaConf::load();
            confLoaded = true;
            applyConfig();
        }
    }
    // --- Hyperparameters (set from config or defaults) ---
    int    popSize       = NUM_CHROMOSOMES;
    double mutationRate  = 0.08;
    double crossoverRate = 0.85;
    double lambda        = DEFAULT_X;
    Weights weights      = { DEFAULT_W1, DEFAULT_W2, DEFAULT_W3, DEFAULT_W4, DEFAULT_X };

    // Output paths (overridden by config when present)
    std::string coverageOutput = "final_optimized_coverage.png";
    std::string scatterOutput  = "final_optimized_scatter.png";

    // Physical spec for the homogenous evaluator (all radars share these values)
    RadarConf homoParams;

    // --- Fixed-length operator selectors ---
    SelectionMethod    currentSelection  = SelectionMethod::TOURNAMENT;
    CrossoverMethod    currentCrossover  = CrossoverMethod::UNIFORM;
    MutationMethod     currentMutation   = MutationMethod::SPATIAL_NUDGE;

    // --- Variable-length operator selectors ---
    SelectionMethod     autoSelection    = SelectionMethod::TOURNAMENT;
    AutoCrossoverMethod autoCrossover    = AutoCrossoverMethod::RANDOM_SUBSET;
    double autoMutationRate  = 0.10;
    double autoCrossoverRate = 0.85;

    std::mt19937 rng;

    // =========================================================
    // SECTION 0 — CONFIG LOADING
    // =========================================================

    // Apply all tunable parameters from a successfully loaded GaConf.
    // Operator enums are left to setOperators() / setAutoOperators() so the
    // caller can still override them after construction if desired.
    void applyConfig() {
        const auto& p = conf.parameters;

        popSize      = (p.numPopulation  > 0) ? p.numPopulation  : popSize;
        mutationRate = (p.mutationRate   > 0) ? p.mutationRate   : mutationRate;
        crossoverRate= (p.crossoverRate  > 0) ? p.crossoverRate  : crossoverRate;
        weights      = p.weights;
        lambda       = p.weights.x;     // x doubles as the radar-count cost λ

        if (!p.output.empty())        coverageOutput = p.output;
        if (!p.scatterOutput.empty()) scatterOutput  = p.scatterOutput;

        // Physical radar params for the homogenous run
        if (conf.homo.is_case) {
            const auto& r = conf.homo.radar;
            if (r.range        > 0) homoParams.range        = r.range * 1000.0; // km → m
            if (r.antennaHeight> 0) homoParams.antennaHeight = r.antennaHeight;
            if (r.target       > 0) homoParams.target  = r.target;
            homoParams.minAngle = r.minAngle;  // float, always valid
            homoParams.maxAngle = r.maxAngle;
        }

        std::cout << "[ga] Config loaded — "
                  << "pop=" << popSize
                  << "  mutRate=" << mutationRate
                  << "  cxRate="  << crossoverRate
                  << "  λ="       << lambda
                  << "  output=\"" << coverageOutput << "\"\n";
    }

    // =========================================================
    // SECTION 1 — FIXED-LENGTH GA (run_homogenous / run_heterogenous)
    // =========================================================

    Chromosome generate_random_chromosome(const FeasibleGrid& f, int length) {
        int max_cells = (int)f.cells.size();
        if (length > max_cells) length = max_cells;

        std::uniform_int_distribution<int> gen(0, max_cells - 1);
        Chromosome cn(length);
        std::unordered_set<int> used;

        for (int i = 0; i < length; i++) {
            int x;
            do { x = gen(rng); } while (!used.insert(x).second);
            cn.radar[i] = {f.cells[x].r, f.cells[x].c};
        }
        return cn;
    }

    // Homogenous fitness — weights and radar physical spec come from the engine.
    // Both are passed by value so worker threads have private copies (no races).
    static void evaluate_fitness(
        Global& g, const FeasibleGrid& base_grid,
        Chromosome& cn, int index,
        const Weights& w, const RadarConf& rp, bool writePng, vector<RadarConf> radars = {},
        const FeasibleGrid& aoi_grid = {})   // full AOI (allCells=WITHALL); fitness now scored against this
    {
        std::string png_path = writePng
                   ? "chromosome_" + std::to_string(index) + "_coverage.png" : "";

         std::vector<int> res;
        if(radars.empty()) {
                res = calculate_radar_coverage_overlap(
                   g, base_grid, cn.radar,
                   5.0,
                   png_path,
                   rp.antennaHeight,
                   rp.minAngle,
                   rp.target,
                   rp.minAngle,
                   rp.range,
                   {},
                   aoi_grid);
        }else
        {
            res = calculate_radar_coverage_overlap(
                g, base_grid, cn.radar, 5.0, png_path,0,0,0,0,0, radars, aoi_grid);
        }

               // total denominator matches whatever calculate_radar_coverage_overlap scored
               // result[] against: full AOI when provided, feasible-only otherwise.
               int    total        = !aoi_grid.cells.empty()
                                    ? (int)aoi_grid.cells.size()
                                    : (int)base_grid.cells.size();
               double unique       = (double)res[1] / total;
               double low_overlap  = (double)res[2] / total;
               double high_overlap = (double)res[3] / total;
               double spacing_gaps = (double)res[0] / total;
               double overlap_factor = (low_overlap + high_overlap) / unique;
               bool f_type_four_weights    = false;

               if (f_type_four_weights == true){
                   cn.fitness = (unique       * w.w1)
                              + (low_overlap  * w.w2)
                              - (high_overlap * w.w3)
                              - (spacing_gaps * w.w4);
               }
               else{
                   double target_error = abs(0.25 - overlap_factor);
                   double closeness_reward = std::max<double>(0.0, 1.0 - target_error);
                   cn.fitness = (0.70 * unique) + (0.30 * closeness_reward);
               }
    }

    // Heterogenous fitness — weighted cell-count formula (case 2).
    // Weights struct kept for interface consistency; case-2 formula uses
    // fixed multipliers by design, but w4 is reused for the proximity penalty.


    // --- Parallel evaluation (homogenous) ---
    void parallelEvaluatePopulation(
        std::vector<Chromosome>& pop, Global& globalRef,
        const FeasibleGrid& base_grid, bool writePng = false,
        const FeasibleGrid& aoi_grid = {})
    {
        if(conf.homo.is_case) {
            auto startTime = std::chrono::high_resolution_clock::now();
           std::vector<std::future<void>> futures;
           std::vector<RadarConf> empty(0);
           futures.reserve(pop.size());
           for (size_t i = 0; i < pop.size(); i++)
               futures.push_back(pool.enqueue(
                   evaluate_fitness,
                   std::ref(globalRef), std::ref(base_grid),
                   std::ref(pop[i]), (int)i,
                   weights, homoParams, writePng, empty, std::ref(aoi_grid)));
           for (auto& f : futures) f.wait();
           auto endTime = std::chrono::high_resolution_clock::now();
           std::chrono::duration<double, std::milli> duration = endTime - startTime;
           totalFitnessEvalTimeMs += duration.count();
        }else {
            auto startTime = std::chrono::high_resolution_clock::now();
           std::vector<std::future<void>> futures;
           futures.reserve(pop.size());
           for (size_t i = 0; i < pop.size(); i++)
               futures.push_back(pool.enqueue(
                   evaluate_fitness,
                   std::ref(globalRef), std::ref(base_grid),
                   std::ref(pop[i]), (int)i,
                   weights, homoParams, writePng, std::ref(conf.hetero.radars), std::ref(aoi_grid)));  // homoParams carries correct rangeM
           for (auto& f : futures) f.wait();
           auto endTime = std::chrono::high_resolution_clock::now();
           std::chrono::duration<double, std::milli> duration = endTime - startTime;
           totalFitnessEvalTimeMs += duration.count();
        }

    }

    // --- Parallel evaluation (heterogenous) ---


    GridPoint snapToFeasible(int r, int c, const FeasibleGrid& f) {
        int minDistSq = INT_MAX;
        GridPoint best = {f.cells.front().r, f.cells.front().c};
        for (const auto& cell : f.cells) {
            int dr = cell.r - r, dc = cell.c - c;
            int d  = dr*dr + dc*dc;
            if (d < minDistSq) {
                minDistSq = d;
                best = {cell.r, cell.c};
                if (d == 0) break;
            }
        }
        return best;
    }

    // --- SELECTION (fixed-length) ---

    Chromosome tournamentSelection(const std::vector<Chromosome>& pop) {
        std::uniform_int_distribution<int> dist(0, popSize - 1);
        Chromosome best = pop[dist(rng)];
        for (int i = 0; i < 2; ++i) {
            Chromosome comp = pop[dist(rng)];
            if (comp.fitness > best.fitness) best = comp;
        }
        return best;
    }

    Chromosome rouletteWheelSelection(const std::vector<Chromosome>& pop) {
        double minFit = pop[0].fitness;
        for (const auto& ind : pop) if (ind.fitness < minFit) minFit = ind.fitness;
        double offset = (minFit < 0) ? std::abs(minFit) + 1.0 : 0.0;

        double total = 0.0;
        for (const auto& ind : pop) total += (ind.fitness + offset);

        std::uniform_real_distribution<double> dist(0.0, total);
        double slice = dist(rng), accum = 0.0;
        for (const auto& ind : pop) {
            accum += (ind.fitness + offset);
            if (accum >= slice) return ind;
        }
        return pop.back();
    }

    Chromosome rankBasedSelection(const std::vector<Chromosome>& sorted_pop) {
        int n = (int)sorted_pop.size();
        int totalRankSum = (n * (n + 1)) / 2;
        std::uniform_int_distribution<int> dist(1, totalRankSum);
        int target = dist(rng), rankSum = 0;
        for (int i = 0; i < n; ++i) {
            rankSum += (i + 1);
            if (rankSum >= target) return sorted_pop[i];
        }
        return sorted_pop.back();
    }

    std::vector<Chromosome> stochasticUniversalSampling(
        const std::vector<Chromosome>& pop, int numToSelect)
    {
        double minFit = pop[0].fitness;
        for (const auto& ind : pop) if (ind.fitness < minFit) minFit = ind.fitness;
        double offset = (minFit < 0) ? std::abs(minFit) + 1.0 : 0.0;

        double total = 0.0;
        for (const auto& ind : pop) total += (ind.fitness + offset);

        double ptr_dist  = total / numToSelect;
        std::uniform_real_distribution<double> distR(0.0, ptr_dist);
        double start_ptr = distR(rng);

        std::vector<Chromosome> selected;
        double accum = 0.0;
        size_t idx = 0;
        for (int i = 0; i < numToSelect; ++i) {
            double ptr = start_ptr + i * ptr_dist;
            while (idx < pop.size()-1 &&
                   accum + (pop[idx].fitness + offset) < ptr) {
                accum += (pop[idx].fitness + offset);
                idx++;
            }
            selected.push_back(pop[idx]);
        }
        return selected;
    }

    Chromosome executeSelection(const std::vector<Chromosome>& pop) {
        if (currentSelection == SelectionMethod::TOURNAMENT)     return tournamentSelection(pop);
        if (currentSelection == SelectionMethod::ROULETTE_WHEEL) return rouletteWheelSelection(pop);
        return rankBasedSelection(pop);
    }

    // --- CROSSOVER (fixed-length) ---

    std::pair<Chromosome, Chromosome> executeCrossover(
        const Chromosome& p1, const Chromosome& p2)
    {
        std::uniform_real_distribution<double> distR(0.0, 1.0);
        if (distR(rng) > crossoverRate) return {p1, p2};

        int len = (int)p1.radar.size();
        Chromosome c1(0), c2(0);

        if (currentCrossover == CrossoverMethod::SINGLE_POINT) {
            std::uniform_int_distribution<int> distP(1, len - 2);
            int cp = distP(rng);
            c1.radar.insert(c1.radar.end(), p1.radar.begin(), p1.radar.begin() + cp);
            c1.radar.insert(c1.radar.end(), p2.radar.begin() + cp, p2.radar.end());
            c2.radar.insert(c2.radar.end(), p2.radar.begin(), p2.radar.begin() + cp);
            c2.radar.insert(c2.radar.end(), p1.radar.begin() + cp, p1.radar.end());
        }
        else if (currentCrossover == CrossoverMethod::MULTI_POINT) {
            std::uniform_int_distribution<int> distP(1, len - 2);
            int cp1 = distP(rng), cp2 = distP(rng);
            if (cp1 > cp2) std::swap(cp1, cp2);
            for (int i = 0; i < len; ++i) {
                if (i >= cp1 && i < cp2) {
                    c1.radar.push_back(p2.radar[i]); c2.radar.push_back(p1.radar[i]);
                } else {
                    c1.radar.push_back(p1.radar[i]); c2.radar.push_back(p2.radar[i]);
                }
            }
        }
        else { // UNIFORM
            for (int i = 0; i < len; ++i) {
                if (distR(rng) < 0.5) {
                    c1.radar.push_back(p1.radar[i]); c2.radar.push_back(p2.radar[i]);
                } else {
                    c1.radar.push_back(p2.radar[i]); c2.radar.push_back(p1.radar[i]);
                }
            }
        }
        return {c1, c2};
    }

    // --- MUTATION (fixed-length) ---

    void executeMutation(
        Chromosome& ind, const FeasibleGrid& f,
        std::normal_distribution<double>& gauss)
    {
        std::uniform_real_distribution<double> distR(0.0, 1.0);
        std::uniform_int_distribution<int> distF(0, (int)f.cells.size() - 1);
        std::uniform_int_distribution<int> distNudge(-5, 5);

        for (auto& gene : ind.radar) {
            if (distR(rng) < mutationRate) {
                if (currentMutation == MutationMethod::SINGLE_POINT) {
                    gene = {f.cells[distF(rng)].r, f.cells[distF(rng)].c};
                }
                else if (currentMutation == MutationMethod::SPATIAL_NUDGE) {
                    gene = snapToFeasible(gene.r + distNudge(rng), gene.c + distNudge(rng), f);
                }
                else { // GAUSSIAN_CREEP
                    int rOff = (int)std::round(gauss(rng));
                    int cOff = (int)std::round(gauss(rng));
                    gene = snapToFeasible(gene.r + rOff, gene.c + cOff, f);
                }
            }
        }
    }

    // =========================================================
    // SECTION 2 — VARIABLE-LENGTH GA (run_auto)
    // =========================================================

    static int64_t posKey(int r, int c) {
        return static_cast<int64_t>(r) * 1'000'000LL + c;
    }

    static std::vector<GridPoint> extractPositions(const std::vector<RadarGene>& genes) {
        std::vector<GridPoint> pts;
        pts.reserve(genes.size());
        for (const auto& gene : genes)
            if (gene.typeIndex >= 0) pts.push_back({gene.r, gene.c});
        return pts;
    }

    static std::vector<double> extractRangesM(
        const std::vector<RadarGene>& genes,
        const std::vector<FlexRadarConfig>& configs)
    {
        std::vector<double> ranges;
        ranges.reserve(genes.size());
        for (const auto& gene : genes)
            if (gene.typeIndex >= 0 && gene.typeIndex < (int)configs.size())
                ranges.push_back(configs[gene.typeIndex].rangeKm * 1000.0);
        return ranges;
    }

    AutoChromosome generateAutoChromosome(
        const FeasibleGrid& f, const std::vector<FlexRadarConfig>& configs)
    {
        AutoChromosome cn;
        std::unordered_set<int64_t> occupied;
        std::uniform_int_distribution<int> distF(0, (int)f.cells.size() - 1);

        for (int t = 0; t < (int)configs.size(); t++) {
            int count = std::uniform_int_distribution<int>(
                configs[t].minCount, configs[t].maxCount)(rng);

            for (int k = 0; k < count; k++) {
                for (int attempt = 0; attempt < 500; attempt++) {
                    const auto& cell = f.cells[distF(rng)];
                    int64_t key = posKey(cell.r, cell.c);
                    if (!occupied.count(key)) {
                        cn.radars.push_back({cell.r, cell.c, t});
                        occupied.insert(key);
                        break;
                    }
                }
            }
        }
        return cn;
    }

    // Auto fitness — same layered-count formula as heterogenous, plus λ·n penalty.
    // Weights::w4 scales the proximity penalty; weights::x is λ (radar-count cost).
    // static void evaluate_auto_fitness(
    //     Global& g, const FeasibleGrid& base_grid,
    //     AutoChromosome& cn, int index,
    //     const std::vector<FlexRadarConfig>& configs,
    //     const Weights& w, bool writePng)
    // {
    //     auto positions = extractPositions(cn.radars);
    //     auto ranges    = extractRangesM(cn.radars, configs);

    //     std::string png_path = writePng
    //         ? "chromosome_" + std::to_string(index) + "_coverage.png" : "";

    //     auto res = calculate_radar_coverage_overlap(
    //         g, base_grid, positions, 5.0, png_path, ranges);

    //     const double penaltyScale = (w.w4 > 0.0) ? w.w4 * 5000.0 / DEFAULT_W4 : 5000.0;
    //     int n = cn.totalCount();
    //     cn.fitness = (res[1] * 1.0) + (res[2] * 2.0) + (res[3] * 3.0)
    //                - (res[0] * penaltyScale)
    //                - w.x * n;           // λ sourced from weights.x
    // }

    // void parallelEvaluateAutoPopulation(
    //     std::vector<AutoChromosome>& pop, Global& globalRef,
    //     const FeasibleGrid& base_grid,
    //     const std::vector<FlexRadarConfig>& configs,
    //     bool writePng = false)
    // {
    //     std::vector<std::future<void>> futures;
    //     futures.reserve(pop.size());
    //     for (size_t i = 0; i < pop.size(); i++)
    //         futures.push_back(pool.enqueue(
    //             evaluate_auto_fitness,
    //             std::ref(globalRef), std::ref(base_grid),
    //             std::ref(pop[i]), (int)i,
    //             std::ref(configs), weights, writePng));
    //     for (auto& f : futures) f.wait();
    // }

    // --- REPAIR ---
    void repairChromosome(
        AutoChromosome& ind, const FeasibleGrid& f,
        const std::vector<FlexRadarConfig>& configs)
    {
        // 1. Remove NULL sentinels
        ind.radars.erase(
            std::remove_if(ind.radars.begin(), ind.radars.end(),
                [](const RadarGene& g) { return g.typeIndex < 0; }),
            ind.radars.end());

        // 2. Deduplicate by position
        {
            std::unordered_set<int64_t> seen;
            ind.radars.erase(
                std::remove_if(ind.radars.begin(), ind.radars.end(),
                    [&](const RadarGene& g) {
                        return !seen.insert(posKey(g.r, g.c)).second;
                    }),
                ind.radars.end());
        }

        // 3. Trim types that exceed maxCount
        for (int t = 0; t < (int)configs.size(); t++) {
            int excess = ind.countType(t) - configs[t].maxCount;
            while (excess > 0) {
                std::vector<int> byType;
                for (int i = 0; i < (int)ind.radars.size(); i++)
                    if (ind.radars[i].typeIndex == t) byType.push_back(i);
                if (byType.empty()) break;
                int rem = byType[
                    std::uniform_int_distribution<int>(0, (int)byType.size()-1)(rng)];
                ind.radars.erase(ind.radars.begin() + rem);
                excess--;
            }
        }

        // 4. Fill types below minCount
        {
            std::unordered_set<int64_t> occupied;
            for (const auto& gene : ind.radars)
                occupied.insert(posKey(gene.r, gene.c));

            std::uniform_int_distribution<int> distF(0, (int)f.cells.size() - 1);

            for (int t = 0; t < (int)configs.size(); t++) {
                int deficit = configs[t].minCount - ind.countType(t);
                while (deficit > 0) {
                    bool placed = false;
                    for (int attempt = 0; attempt < 500; attempt++) {
                        const auto& cell = f.cells[distF(rng)];
                        int64_t key = posKey(cell.r, cell.c);
                        if (!occupied.count(key)) {
                            ind.radars.push_back({cell.r, cell.c, t});
                            occupied.insert(key);
                            placed = true;
                            break;
                        }
                    }
                    if (!placed) break;
                    deficit--;
                }
            }
        }
    }

    // --- MUTATION (variable-length) ---

    void mutate_add(AutoChromosome& ind, const FeasibleGrid& f,
                    const std::vector<FlexRadarConfig>& configs)
    {
        std::vector<int> eligible;
        for (int t = 0; t < (int)configs.size(); t++)
            if (ind.countType(t) < configs[t].maxCount) eligible.push_back(t);
        if (eligible.empty()) return;

        int typeIdx = eligible[
            std::uniform_int_distribution<int>(0, (int)eligible.size()-1)(rng)];

        std::unordered_set<int64_t> occupied;
        for (const auto& gene : ind.radars)
            if (gene.typeIndex >= 0) occupied.insert(posKey(gene.r, gene.c));

        std::uniform_int_distribution<int> distF(0, (int)f.cells.size()-1);
        for (int attempt = 0; attempt < 300; attempt++) {
            const auto& cell = f.cells[distF(rng)];
            if (!occupied.count(posKey(cell.r, cell.c))) {
                ind.radars.push_back({cell.r, cell.c, typeIdx});
                return;
            }
        }
    }

    void mutate_remove(AutoChromosome& ind,
                       const std::vector<FlexRadarConfig>& configs)
    {
        std::vector<int> removable;
        for (int i = 0; i < (int)ind.radars.size(); i++) {
            int t = ind.radars[i].typeIndex;
            if (t >= 0 && ind.countType(t) > configs[t].minCount)
                removable.push_back(i);
        }
        if (removable.empty()) return;

        int idx = removable[
            std::uniform_int_distribution<int>(0, (int)removable.size()-1)(rng)];
        ind.radars.erase(ind.radars.begin() + idx);
    }

    void mutate_move(AutoChromosome& ind, const FeasibleGrid& f) {
        if (ind.radars.empty()) return;

        auto& gene = ind.radars[
            std::uniform_int_distribution<int>(0, (int)ind.radars.size()-1)(rng)];
        if (gene.typeIndex < 0) return;

        std::uniform_int_distribution<int> distNudge(-10, 10);
        GridPoint snapped = snapToFeasible(
            gene.r + distNudge(rng), gene.c + distNudge(rng), f);
        gene.r = snapped.r;
        gene.c = snapped.c;
    }

    void executeMutationAuto(AutoChromosome& ind, const FeasibleGrid& f,
                             const std::vector<FlexRadarConfig>& configs)
    {
        if (std::uniform_real_distribution<double>(0.0, 1.0)(rng) > autoMutationRate)
            return;
        switch (std::uniform_int_distribution<int>(0, 2)(rng)) {
            case 0: mutate_add(ind, f, configs);    break;
            case 1: mutate_remove(ind, configs);    break;
            case 2: mutate_move(ind, f);            break;
        }
    }

    // --- CROSSOVER (variable-length) ---

    std::pair<AutoChromosome, AutoChromosome> crossover_randomSubset(
        const AutoChromosome& p1, const AutoChromosome& p2)
    {
        auto splitParent = [&](const AutoChromosome& parent)
            -> std::pair<std::vector<RadarGene>, std::vector<RadarGene>>
        {
            int n = (int)parent.radars.size();
            int k = (n > 0) ? std::uniform_int_distribution<int>(0, n)(rng) : 0;
            std::vector<int> idx(n);
            std::iota(idx.begin(), idx.end(), 0);
            std::shuffle(idx.begin(), idx.end(), rng);
            std::vector<RadarGene> chosen, rest;
            for (int i = 0; i < n; i++)
                (i < k ? chosen : rest).push_back(parent.radars[idx[i]]);
            return {chosen, rest};
        };

        auto [a_chosen, a_rest] = splitParent(p1);
        auto [b_chosen, b_rest] = splitParent(p2);

        AutoChromosome c1, c2;
        c1.radars.insert(c1.radars.end(), a_chosen.begin(), a_chosen.end());
        c1.radars.insert(c1.radars.end(), b_chosen.begin(), b_chosen.end());
        c2.radars.insert(c2.radars.end(), a_rest.begin(),   a_rest.end());
        c2.radars.insert(c2.radars.end(), b_rest.begin(),   b_rest.end());
        return {c1, c2};
    }

    std::pair<AutoChromosome, AutoChromosome> crossover_lengthNormalized(
        const AutoChromosome& p1, const AutoChromosome& p2)
    {
        int avgLen = ((int)p1.radars.size() + (int)p2.radars.size() + 1) / 2;
        std::vector<RadarGene> pool;
        pool.insert(pool.end(), p1.radars.begin(), p1.radars.end());
        pool.insert(pool.end(), p2.radars.begin(), p2.radars.end());
        int take = std::min(avgLen, (int)pool.size());

        auto buildChild = [&]() {
            std::shuffle(pool.begin(), pool.end(), rng);
            AutoChromosome child;
            child.radars.assign(pool.begin(), pool.begin() + take);
            return child;
        };
        return {buildChild(), buildChild()};
    }

    std::pair<AutoChromosome, AutoChromosome> crossover_positionalPadding(
        const AutoChromosome& p1, const AutoChromosome& p2)
    {
        int maxLen = (int)std::max(p1.radars.size(), p2.radars.size());
        if (maxLen == 0) return {p1, p2};

        const RadarGene null_gene = {0, 0, -1};
        std::vector<RadarGene> pad1(p1.radars), pad2(p2.radars);
        while ((int)pad1.size() < maxLen) pad1.push_back(null_gene);
        while ((int)pad2.size() < maxLen) pad2.push_back(null_gene);

        int cut = std::uniform_int_distribution<int>(1, maxLen)(rng);
        AutoChromosome c1, c2;
        for (int i = 0; i < maxLen; i++) {
            if (i < cut) { c1.radars.push_back(pad1[i]); c2.radars.push_back(pad2[i]); }
            else          { c1.radars.push_back(pad2[i]); c2.radars.push_back(pad1[i]); }
        }

        auto stripNull = [](std::vector<RadarGene>& v) {
            v.erase(std::remove_if(v.begin(), v.end(),
                [](const RadarGene& g) { return g.typeIndex < 0; }), v.end());
        };
        stripNull(c1.radars);
        stripNull(c2.radars);
        return {c1, c2};
    }

    std::pair<AutoChromosome, AutoChromosome> executeCrossoverAuto(
        const AutoChromosome& p1, const AutoChromosome& p2)
    {
        if (std::uniform_real_distribution<double>(0.0, 1.0)(rng) > autoCrossoverRate)
            return {p1, p2};
        switch (std::uniform_int_distribution<int>(0, 2)(rng)) {
            case 0:  return crossover_randomSubset(p1, p2);
            case 1:  return crossover_lengthNormalized(p1, p2);
            default: return crossover_positionalPadding(p1, p2);
        }
    }

    // --- SELECTION (variable-length) ---

    AutoChromosome tournamentSelectionAuto(const std::vector<AutoChromosome>& pop) {
        std::uniform_int_distribution<int> dist(0, (int)pop.size()-1);
        AutoChromosome best = pop[dist(rng)];
        for (int i = 0; i < 2; ++i) {
            const AutoChromosome& comp = pop[dist(rng)];
            if (comp.fitness > best.fitness) best = comp;
        }
        return best;
    }

    AutoChromosome rouletteWheelSelectionAuto(const std::vector<AutoChromosome>& pop) {
        double minFit = pop[0].fitness;
        for (const auto& ind : pop) if (ind.fitness < minFit) minFit = ind.fitness;
        double offset = (minFit < 0) ? std::abs(minFit) + 1.0 : 0.0;

        double total = 0.0;
        for (const auto& ind : pop) total += (ind.fitness + offset);

        double slice = std::uniform_real_distribution<double>(0.0, total)(rng);
        double accum = 0.0;
        for (const auto& ind : pop) {
            accum += (ind.fitness + offset);
            if (accum >= slice) return ind;
        }
        return pop.back();
    }

    AutoChromosome rankBasedSelectionAuto(const std::vector<AutoChromosome>& sorted_pop) {
        int n = (int)sorted_pop.size();
        int totalRankSum = (n * (n + 1)) / 2;
        int target = std::uniform_int_distribution<int>(1, totalRankSum)(rng);
        int rankSum = 0;
        for (int i = 0; i < n; ++i) {
            rankSum += (i + 1);
            if (rankSum >= target) return sorted_pop[i];
        }
        return sorted_pop.back();
    }

    std::vector<AutoChromosome> stochasticUniversalSamplingAuto(
        const std::vector<AutoChromosome>& pop, int numToSelect)
    {
        double minFit = pop[0].fitness;
        for (const auto& ind : pop) if (ind.fitness < minFit) minFit = ind.fitness;
        double offset = (minFit < 0) ? std::abs(minFit) + 1.0 : 0.0;

        double total = 0.0;
        for (const auto& ind : pop) total += (ind.fitness + offset);

        double ptr_dist  = total / numToSelect;
        double start_ptr = std::uniform_real_distribution<double>(0.0, ptr_dist)(rng);

        std::vector<AutoChromosome> selected;
        double accum = 0.0;
        size_t idx   = 0;
        for (int i = 0; i < numToSelect; ++i) {
            double ptr = start_ptr + i * ptr_dist;
            while (idx < pop.size()-1 &&
                   accum + (pop[idx].fitness + offset) < ptr) {
                accum += (pop[idx].fitness + offset);
                idx++;
            }
            selected.push_back(pop[idx]);
        }
        return selected;
    }

    AutoChromosome executeSelectionAuto(const std::vector<AutoChromosome>& pop) {
        if (autoSelection == SelectionMethod::TOURNAMENT)     return tournamentSelectionAuto(pop);
        if (autoSelection == SelectionMethod::ROULETTE_WHEEL) return rouletteWheelSelectionAuto(pop);
        return rankBasedSelectionAuto(pop);
    }

    // =========================================================
    // SECTION 3 — EVOLUTION ENTRY POINTS
    // =========================================================

public:
    explicit GaEngine(ThreadPool& p) : pool(p) {
        rng.seed((unsigned int)
            std::chrono::system_clock::now().time_since_epoch().count());

        // Load config if available; otherwise the defaults declared above are used
        if (GaConf::is_available()) {
            conf = GaConf::load();
            confLoaded = true;
            applyConfig();
        } else {
            std::cout << "[ga] No config file found — using built-in defaults.\n";
        }

        std::cout << "[ga] Synchronously pre-loading all .bin matrices to memory...\n";
        g.bin.init();
    }

    // Override GA operators for fixed-length runs.
    // Called from main() or caller code after construction.
    void setOperators(SelectionMethod sm, CrossoverMethod cm,
                      MutationMethod mm, double mRate, double cRate)
    {
        currentSelection = sm;
        currentCrossover = cm;
        currentMutation  = mm;
        // Only override rates if no config was loaded (config takes priority)
        if (!confLoaded) {
            mutationRate  = mRate;
            crossoverRate = cRate;
        }
    }

    // Override GA operators for variable-length runs.
    void setAutoOperators(SelectionMethod sm, AutoCrossoverMethod cm,
                          double mRate, double cRate)
    {
        autoSelection = sm;
        autoCrossover = cm;
        if (!confLoaded) {
            autoMutationRate  = mRate;
            autoCrossoverRate = cRate;
        }
    }

    // Explicit lambda setter; ignored when config provides weights.x
    void setLambda(double l) {
        if (!confLoaded) lambda = l;
    }

    const vector<Chromosome>& get_population() const { return population; }

    // generate random inputs for heatmap
    void run_native_random_simulation(int num_experiments) {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

        std::default_random_engine generator(seed);
        std::uniform_real_distribution<double> dist_range(100.0, 150.0);
        std::uniform_real_distribution<double> dist_antenna(10.0, 50.0);
        std::uniform_real_distribution<double> dist_target(500.0, 5000.0);
        std::uniform_real_distribution<double> dist_min_angle(0.0, 5.0);
        std::uniform_real_distribution<double> dist_max_angle(40.0, 60.0);
        std::uniform_real_distribution<double> dist_num_radars(3.0, 5.0);

        std::cout << "--> Starting " << num_experiments << " Native Randomized Homogeneous Simulations..." << std::endl;
        ThreadPool pool;
        GaEngine ga(pool);

        ga.setOperators(SelectionMethod::TOURNAMENT, CrossoverMethod::UNIFORM,
                        MutationMethod::SPATIAL_NUDGE,
                        /*mRate=*/0.08, /*cRate=*/0.85);

        string CONFIG_FILE = "gaconf.json";

        for (int i = 1; i <= num_experiments; ++i) {
            std::cout << "\n==================================================" << std::endl;
            std::cout << " RUN SIMULATION [" << i << " / " << num_experiments << "]" << std::endl;
            std::cout << "==================================================" << std::endl;

            // 2. Read the current configuration file into a Json::Value object
            Json::Value root;
            std::ifstream file_input(CONFIG_FILE);
            if (!file_input.is_open()) {
                std::cerr << "❌ Error: Could not open config.json!" << std::endl;
                return;
            }

            Json::CharReaderBuilder reader;
            std::string errs;
            if (!Json::parseFromStream(reader, file_input, &root, &errs)) {
                std::cerr << "❌ JSON Parsing Error: " << errs << std::endl;
                return;
            }
            file_input.close(); // Close immediately after reading

            double rand_range    = dist_range(generator);
            double rand_antenna  = dist_antenna(generator);
            double rand_target   = dist_target(generator);
            double rand_min_angle = dist_min_angle(generator);
            double rand_max_angle = dist_max_angle(generator);
            double rand_num_radars = dist_num_radars(generator);

            root["heterogenous"]["case"] = false;

            root["homogenous"]["case"] = true;
            root["homogenous"]["chromosome-length"]  = int(rand_num_radars);
            root["homogenous"]["radar"]["range"]     = int(rand_range);
            root["homogenous"]["radar"]["antenna"]   = int(rand_antenna);
            root["homogenous"]["radar"]["target"]    = int(rand_target);
            root["homogenous"]["radar"]["min-angle"] = int(rand_min_angle);
            root["homogenous"]["radar"]["max-angle"] = int(rand_max_angle);

            // Change output file names per run dynamically so they don't overwrite each other
            root["output"]         = "result_run_" + std::to_string(i) + ".png";
            root["scatter-output"] = "scatter_run_" + std::to_string(i) + ".png";

            // 5. Write the updated JSON structure back to config.json
            std::ofstream file_output(CONFIG_FILE);
            Json::StreamWriterBuilder writer;
            file_output << Json::writeString(writer, root);
            file_output.close();

            std::cout << ">> Config updated with random parameters:" << std::endl;
            std::cout << "   - Range: " << int(rand_range) << " km\n"
                      << "   - Antenna: " << int(rand_antenna) << " m\n"
                      << "   - No. of Radars: " << int(rand_num_radars) << " \n"
                      << "   - Target Threat: " << int(rand_target) << " m\n"
                      << "   - Min Angle: " << int(rand_min_angle) << " \n"
                      << "   - Max Angle: " << int(rand_max_angle) << " " << std::endl;
            ga.reloadConfig();
            ga.run_homogenous();
            std::cout << ">> Run [" << i << "] completed and saved successfully." << std::endl;
        }
    }
    // ------------------------------------------------------------------
    // run_homogenous
    // All radars share a single fixed range (set in radar.h defaults).
    // chromosomeLength and generations sourced from config when available.
    // ------------------------------------------------------------------
    double run_homogenous(int pathLength   = 2,
                        int generations  = 50,
                        bool opt         = false,
                        int  num_top_k   = 15)
    {
        // Config overrides caller-supplied values when present
        if (confLoaded && conf.homo.is_case) {
            if (conf.homo.chromosomeLength > 0)
                pathLength  = conf.homo.chromosomeLength;
            if (conf.parameters.numGenerations > 0){
                generations = conf.parameters.numGenerations;
            }
            if (conf.parameters.save)
                opt = true;
            if (conf.parameters.saveSolutions > 0)
                num_top_k = conf.parameters.saveSolutions;
        }

        population.clear();
        population.reserve(popSize);

        std::cout << "[ga] [homogenous] chromosomes=" << pathLength
                  << "  generations=" << generations
                  << "  popSize=" << popSize << "\n";
        std::cout << "[ga] Extracting baseline Area of Operation boundaries...\n";
        auto f = Bin::getFeasibleGrid(g.bin, WRITE);
        auto aoiGrid = Bin::getFeasibleGrid(g.bin, 0, 0, WITHALL); // unfiltered AOI, same AABB, for PNG background
         auto initStart = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < popSize; i++)
            population.push_back(generate_random_chromosome(f, pathLength));
        auto initEnd = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> initDuration = initEnd - initStart;

        std::cout << "[ga] Dispatching parallel coverage evaluation for Generation 0...\n";
        parallelEvaluatePopulation(population, g, f, false, aoiGrid);

        std::normal_distribution<double> gaussDist(0.0, 10.0);

        for (int gen = 1; gen <= generations; ++gen) {

            auto stepStart = std::chrono::high_resolution_clock::now();
            std::sort(population.begin(), population.end(),
                [](const Chromosome& a, const Chromosome& b) {
                    return a.fitness < b.fitness; });

            std::vector<Chromosome> nextGen;
            nextGen.reserve(popSize);
            nextGen.push_back(population.back()); // elitism

            bool useSUS = (currentSelection == SelectionMethod::SUS);
            std::vector<Chromosome> matingPool;
            size_t poolIdx = 0;
            if (useSUS)
                matingPool = stochasticUniversalSampling(population, popSize - 1);

            while ((int)nextGen.size() < popSize) {
                Chromosome p1(0), p2(0);
                if (useSUS) {
                    p1 = matingPool[poolIdx % matingPool.size()]; poolIdx++;
                    p2 = matingPool[poolIdx % matingPool.size()]; poolIdx++;
                } else {
                    p1 = executeSelection(population);
                    p2 = executeSelection(population);
                }

                auto [c1, c2] = executeCrossover(p1, p2);
                executeMutation(c1, f, gaussDist);
                executeMutation(c2, f, gaussDist);

                nextGen.push_back(c1);
                if ((int)nextGen.size() < popSize) nextGen.push_back(c2);
            }

            auto stepEnd = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::milli> stepDuration = stepEnd - stepStart;
            totalEvolutionaryTimeMs += stepDuration.count();

            parallelEvaluatePopulation(nextGen, g, f, false, aoiGrid);
            population = std::move(nextGen);

            auto best = std::max_element(population.begin(), population.end(),
                [](const Chromosome& a, const Chromosome& b) {
                    return a.fitness < b.fitness; });
            std::cout << "  -> Gen " << gen
                      << " Best Fitness: " << best->fitness << "\n";
        }

        std::sort(population.begin(), population.end(),
            [](const Chromosome& a, const Chromosome& b) {
                return a.fitness > b.fitness; });

        if (opt) {
            vector<GridPoint> top_k;
            for (int i = 0; i < conf.parameters.saveSolutions && i < (int)population.size(); i++)
                for (const auto& pt : population[i].radar)
                    top_k.push_back(pt);
            write_h5_indiviual(g, f, "data.h5", top_k,homoParams.antennaHeight,
            homoParams.target,
            homoParams.minAngle,
            homoParams.maxAngle,
            homoParams.range);
        }

        auto& finalBest = population[0];
        std::cout << "\n[ga] Evolution Complete. Best Lifetime Score: "
                  << finalBest.fitness << "\n";

        calculate_radar_coverage_overlap(
            g, f, finalBest.radar,
            5.0, coverageOutput,
            homoParams.antennaHeight,
            homoParams.minAngle,
            homoParams.target,
            homoParams.maxAngle,
            homoParams.range,
            {},
            aoiGrid);
        Bin::writeRadarScatterPng(finalBest.radar, scatterOutput);
        return finalBest.fitness;
    }

    // ------------------------------------------------------------------
    // run_heterogenous
    // Fixed fleet; each chromosome locus maps to a specific radar range.
    // Fleet sourced from config (hetero.radars) when available.
    // ------------------------------------------------------------------
    void run_heterogenous(std::vector<RadarConfig> fleetConfig = {},
                          int generations = 50)
    {
        if (confLoaded && conf.hetero.is_case) {
            if (conf.parameters.numGenerations > 0)
                generations = conf.parameters.numGenerations;

            // Build fleet from config if none supplied by caller
            if (fleetConfig.empty()) {
                for (const auto& rc : conf.hetero.radars)
                    fleetConfig.push_back({ (double)rc.range, rc.quantity });
            }
        }

        if (fleetConfig.empty()) {
            std::cerr << "[ga] run_heterogenous: no fleet config supplied and "
                         "config unavailable — aborting.\n";
            return ;
        }

        // Flatten fleet into a per-locus range array (metres)
        std::vector<double> flattenedRangesM;
        for (const auto& config : fleetConfig)
            for (int i = 0; i < config.count; ++i)
                flattenedRangesM.push_back(config.rangeKm * 1000.0);

        int pathLength = (int)flattenedRangesM.size();

        std::cout << "[ga] [heterogenous] total radar loci=" << pathLength
                  << "  generations=" << generations
                  << "  popSize=" << popSize << "\n";

        population.clear();
        population.reserve(popSize);

        std::cout << "[ga] Extracting baseline Area of Operation boundaries...\n";
        auto f = Bin::getFeasibleGrid(g.bin, WRITE);
        auto aoiGrid = Bin::getFeasibleGrid(g.bin, 0, 0, WITHALL); // unfiltered AOI, same AABB, for PNG background

        for (int i = 0; i < popSize; i++)
            population.push_back(generate_random_chromosome(f, pathLength));

        std::cout << "[ga] Dispatching parallel coverage evaluation for Generation 0...\n";
        parallelEvaluatePopulation(population, g, f, false, aoiGrid);

        std::normal_distribution<double> gaussDist(0.0, 10.0);

        for (int gen = 1; gen <= generations; ++gen) {
            std::sort(population.begin(), population.end(),
                [](const Chromosome& a, const Chromosome& b) {
                    return a.fitness < b.fitness; });

            std::vector<Chromosome> nextGen;
            nextGen.reserve(popSize);
            nextGen.push_back(population.back()); // elitism

            bool useSUS = (currentSelection == SelectionMethod::SUS);
            std::vector<Chromosome> matingPool;
            size_t poolIdx = 0;
            if (useSUS)
                matingPool = stochasticUniversalSampling(population, popSize - 1);

            while ((int)nextGen.size() < popSize) {
                Chromosome p1(0), p2(0);
                if (useSUS) {
                    p1 = matingPool[poolIdx % matingPool.size()]; poolIdx++;
                    p2 = matingPool[poolIdx % matingPool.size()]; poolIdx++;
                } else {
                    p1 = executeSelection(population);
                    p2 = executeSelection(population);
                }

                auto [c1, c2] = executeCrossover(p1, p2);
                executeMutation(c1, f, gaussDist);
                executeMutation(c2, f, gaussDist);

                nextGen.push_back(c1);
                if ((int)nextGen.size() < popSize) nextGen.push_back(c2);
            }

            parallelEvaluatePopulation(nextGen, g, f, false, aoiGrid);
            population = std::move(nextGen);

            auto best = std::max_element(population.begin(), population.end(),
                [](const Chromosome& a, const Chromosome& b) {
                    return a.fitness < b.fitness; });
            std::cout << "  -> Gen " << gen
                      << " Best Fitness: " << best->fitness << "\n";
        }

        auto finalBest = std::max_element(population.begin(), population.end(),
            [](const Chromosome& a, const Chromosome& b) {
                return a.fitness < b.fitness; });
        std::cout << "\n[ga] Evolution Complete. Best Lifetime Score: "
                  << finalBest->fitness << "\n";

        calculate_radar_coverage_overlap(g, f, finalBest->radar,
                                         5.0, coverageOutput,0,0,0,0,0, conf.hetero.radars,
                                         aoiGrid);
        Bin::writeRadarScatterPng(finalBest->radar, scatterOutput);

    }

    // ------------------------------------------------------------------
    // run_auto
    // Variable-length, variable-count, heterogeneous fleet.
    // Fleet sourced from config (autoG.radars) when available.
    // ------------------------------------------------------------------
    // void run_auto(std::vector<FlexRadarConfig> configs = {},
    //               int generations = 50)
    // {
    //     if (confLoaded && conf.autoG.is_case) {
    //         if (conf.parameters.numGenerations > 0)
    //             generations = conf.parameters.numGenerations;

    //         // Build flex fleet from config if none supplied by caller.
    //         // Each RadarConf.quantity is used as maxCount; minCount is
    //         // clamped to at least 0 and at most quantity.
    //         if (configs.empty()) {
    //             int numTypes = (int)conf.autoG.radars.size();
    //             for (const auto& rc : conf.autoG.radars) {
    //                 int maxC = rc.quantity;
    //                 // Distribute the global min/max evenly across types
    //                 int minC = (numTypes > 0)
    //                     ? std::max(0, conf.autoG.minRadars / numTypes)
    //                     : 0;
    //                 minC = std::min(minC, maxC); // sanity clamp
    //                 configs.push_back({ (double)rc.range, minC, maxC });
    //             }
    //         }
    //     }

    //     if (configs.empty()) {
    //         std::cerr << "[ga] run_auto: no fleet config supplied and "
    //                      "config unavailable — aborting.\n";
    //         return;
    //     }

    //     // Print fleet spec
    //     int totalMin = 0, totalMax = 0;
    //     std::cout << "[ga] [auto] Fleet specification:\n";
    //     for (int t = 0; t < (int)configs.size(); t++) {
    //         std::cout << "  Type " << t << " : " << configs[t].rangeKm
    //                   << " km  [" << configs[t].minCount
    //                   << " – "    << configs[t].maxCount << " units]\n";
    //         totalMin += configs[t].minCount;
    //         totalMax += configs[t].maxCount;
    //     }
    //     std::cout << "  Total radar range : [" << totalMin << ", " << totalMax << "]\n";
    //     std::cout << "  Lambda (cost/radar): " << weights.x << "\n";
    //     std::cout << "  generations=" << generations
    //               << "  popSize=" << popSize << "\n\n";

    //     std::cout << "[ga] Extracting feasible grid...\n";
    //     auto f = Bin::getFeasibleGrid(g.bin, WRITE);

    //     // Seed population
    //     std::vector<AutoChromosome> autoPop;
    //     autoPop.reserve(popSize);
    //     for (int i = 0; i < popSize; i++)
    //         autoPop.push_back(generateAutoChromosome(f, configs));

    //     std::cout << "[ga] Evaluating Generation 0...\n";
    //     parallelEvaluateAutoPopulation(autoPop, g, f, configs, false);

    //     for (int gen = 1; gen <= generations; ++gen) {
    //         std::sort(autoPop.begin(), autoPop.end(),
    //             [](const AutoChromosome& a, const AutoChromosome& b) {
    //                 return a.fitness < b.fitness; });

    //         std::vector<AutoChromosome> nextGen;
    //         nextGen.reserve(popSize);
    //         nextGen.push_back(autoPop.back()); // elitism

    //         bool useSUS = (autoSelection == SelectionMethod::SUS);
    //         std::vector<AutoChromosome> matingPool;
    //         size_t poolIdx = 0;
    //         if (useSUS)
    //             matingPool = stochasticUniversalSamplingAuto(autoPop, popSize - 1);

    //         while ((int)nextGen.size() < popSize) {
    //             AutoChromosome p1, p2;
    //             if (useSUS) {
    //                 p1 = matingPool[poolIdx % matingPool.size()]; poolIdx++;
    //                 p2 = matingPool[poolIdx % matingPool.size()]; poolIdx++;
    //             } else {
    //                 p1 = executeSelectionAuto(autoPop);
    //                 p2 = executeSelectionAuto(autoPop);
    //             }

    //             auto [c1, c2] = executeCrossoverAuto(p1, p2);

    //             repairChromosome(c1, f, configs);
    //             repairChromosome(c2, f, configs);

    //             executeMutationAuto(c1, f, configs);
    //             executeMutationAuto(c2, f, configs);

    //             repairChromosome(c1, f, configs);
    //             repairChromosome(c2, f, configs);

    //             nextGen.push_back(c1);
    //             if ((int)nextGen.size() < popSize) nextGen.push_back(c2);
    //         }

    //         parallelEvaluateAutoPopulation(nextGen, g, f, configs, false);
    //         autoPop = std::move(nextGen);

    //         auto best = std::max_element(autoPop.begin(), autoPop.end(),
    //             [](const AutoChromosome& a, const AutoChromosome& b) {
    //                 return a.fitness < b.fitness; });

    //         std::string breakdown;
    //         for (int t = 0; t < (int)configs.size(); t++)
    //             breakdown += " [T" + std::to_string(t) + ":"
    //                        + std::to_string(best->countType(t)) + "]";

    //         std::cout << "  Gen " << gen
    //                   << " | Fitness: " << best->fitness
    //                   << " | N=" << best->totalCount()
    //                   << breakdown << "\n";
    //     }

    //     auto finalBest = std::max_element(autoPop.begin(), autoPop.end(),
    //         [](const AutoChromosome& a, const AutoChromosome& b) {
    //             return a.fitness < b.fitness; });

    //     std::cout << "\n[ga] Evolution complete.\n";
    //     std::cout << "  Final fitness  : " << finalBest->fitness << "\n";
    //     std::cout << "  Total radars   : " << finalBest->totalCount() << "\n";
    //     for (int t = 0; t < (int)configs.size(); t++)
    //         std::cout << "  Type " << t << " (" << configs[t].rangeKm << " km) : "
    //                   << finalBest->countType(t) << " units\n";

    //     auto positions = extractPositions(finalBest->radars);
    //     auto ranges    = extractRangesM(finalBest->radars, configs);
    //     calculate_radar_coverage_overlap(g, f, positions,
    //                                      5.0, coverageOutput, ranges);
    //     Bin::writeRadarScatterPng(positions, scatterOutput);
    // }
    int run_benchmark() {
        BoundingBox areaOfOps = { 34.0522, 34.1022, -118.2437, -118.1937 };

        // 2. Initialize the shared multi-threaded hardware pool
        ThreadPool sharedPool;
        GaEngine gaEngine(sharedPool);
        // 3. Define 16 diverse test case profiles combining different techniques and parameters
        std::vector<GAConfig> sweepConfigs = {
            // Control & Baseline
            {"Baseline Control",     120, 100, 0.08, 0.85, SelectionMethod::TOURNAMENT,     CrossoverMethod::SINGLE_POINT, MutationMethod::SPATIAL_NUDGE},

            // Lightweight / Fast Profiles
            {"Lightweight (Fast)",    50,  50, 0.05, 0.70, SelectionMethod::TOURNAMENT,     CrossoverMethod::SINGLE_POINT, MutationMethod::SPATIAL_NUDGE},
            {"Micro-Sweep (SUS)",     40,  60, 0.08, 0.80, SelectionMethod::SUS,            CrossoverMethod::MULTI_POINT,  MutationMethod::GAUSSIAN_CREEP},

            // Intensive Searches (Large populations and generations)
            {"Intensive Search",     300, 200, 0.05, 0.80, SelectionMethod::TOURNAMENT,     CrossoverMethod::MULTI_POINT,  MutationMethod::SPATIAL_NUDGE},
            {"Intensive SUS",        250, 150, 0.08, 0.85, SelectionMethod::SUS,            CrossoverMethod::UNIFORM,      MutationMethod::GAUSSIAN_CREEP},

            // Exploration-Heavy (High mutation, uniform crossover, rank-based)
            {"Uniform Exploration",  150, 120, 0.15, 0.90, SelectionMethod::RANK_BASED,     CrossoverMethod::UNIFORM,      MutationMethod::GAUSSIAN_CREEP},
            {"SUS High-Mutation",    120, 100, 0.20, 0.75, SelectionMethod::SUS,            CrossoverMethod::UNIFORM,      MutationMethod::GAUSSIAN_CREEP},

            // Exploitation-Heavy (Low mutation, high crossover)
            {"Exploitation Focus",   150, 100, 0.02, 0.95, SelectionMethod::TOURNAMENT,     CrossoverMethod::MULTI_POINT,  MutationMethod::SPATIAL_NUDGE},

            // Roulette Wheel Variations
            {"Roulette & Uniform",   120, 100, 0.08, 0.85, SelectionMethod::ROULETTE_WHEEL, CrossoverMethod::UNIFORM,      MutationMethod::GAUSSIAN_CREEP},
            {"Roulette & Nudge",     120, 100, 0.06, 0.80, SelectionMethod::ROULETTE_WHEEL, CrossoverMethod::MULTI_POINT,  MutationMethod::SPATIAL_NUDGE},

            // Rank-Based Variations
            {"Rank & Gaussian",      120, 100, 0.08, 0.85, SelectionMethod::RANK_BASED,     CrossoverMethod::MULTI_POINT,  MutationMethod::GAUSSIAN_CREEP},
            {"Rank & Nudge",         120, 100, 0.06, 0.80, SelectionMethod::RANK_BASED,     CrossoverMethod::SINGLE_POINT, MutationMethod::SPATIAL_NUDGE},

            // Single Point Reset Mutation Variations
            {"Tournament & Reset",   120, 100, 0.10, 0.80, SelectionMethod::TOURNAMENT,     CrossoverMethod::SINGLE_POINT, MutationMethod::SINGLE_POINT},
            {"SUS & Reset",          120, 100, 0.10, 0.80, SelectionMethod::SUS,            CrossoverMethod::SINGLE_POINT, MutationMethod::SINGLE_POINT},
            {"Rank & Reset",         120, 100, 0.10, 0.80, SelectionMethod::RANK_BASED,     CrossoverMethod::SINGLE_POINT, MutationMethod::SINGLE_POINT},

            // High Crossover Exploration
            {"High-Cross SUS Nudge", 120, 100, 0.05, 0.95, SelectionMethod::SUS,            CrossoverMethod::SINGLE_POINT, MutationMethod::SPATIAL_NUDGE}
        };

        std::vector<BenchmarkResult> leaderboard;
        leaderboard.reserve(sweepConfigs.size());

        std::cout << "========================================================================\n";
        std::cout << "        STARTING HYPERPARAMETER OPTIMIZATION SWEEP (16 Test Profiles)    \n";
        std::cout << "========================================================================\n";
        std::cout << "Running evaluation sweep over C++ ThreadPool...\n\n";

        for (size_t idx = 0; idx < sweepConfigs.size(); ++idx) {
            const auto& c = sweepConfigs[idx];
            std::cout << " [*] Running Test " << std::setw(2) << (idx + 1) << " / " << sweepConfigs.size()
                      << ": " << std::left << std::setw(25) << c.name << " ... " << std::flush;

            gaEngine.setOperators(c.selection, c.crossover, c.mutation, c.mutationRate, c.crossoverRate);

            auto start = std::chrono::high_resolution_clock::now();
            double resultDistance = gaEngine.run_homogenous();
            auto end = std::chrono::high_resolution_clock::now();

            std::chrono::duration<double, std::milli> elapsed = end - start;

            leaderboard.push_back({
                c.name,
                toString(c.selection),
                toString(c.crossover),
                toString(c.mutation),
                c.popSize,
                c.generations,
                c.mutationRate,
                c.crossoverRate,
                resultDistance,
                elapsed.count(),
                gaEngine.getFitnessEvalTime(),
                gaEngine.getEvolutionaryTime()
            });

            std::cout << "Done. (Distance: " << std::fixed << std::setprecision(3) << resultDistance
                      << " km, Time: " << std::setprecision(1) << elapsed.count() << " ms)\n";
        }

        // Sort leaderboard by best (shortest) distance ascending
        std::sort(leaderboard.begin(), leaderboard.end(), [](const BenchmarkResult& a, const BenchmarkResult& b) {
            return a.bestDistance < b.bestDistance;
            });

        // Print Ranked Performance Results Leaderboard
        std::cout << "\n\n";
        std::cout << "=========================================================================================================================================\n";
        std::cout << "                                                 HYPERPARAMETER OPTIMIZATION LEADERBOARD                                                 \n";
        std::cout << "=========================================================================================================================================\n";
        std::cout << " Rank | " << std::left << std::setw(23) << "Config Profile Name"
                  << " | " << std::setw(11) << "Selection"
                  << " | " << std::setw(12) << "Crossover"
                  << " | " << std::setw(9) << "Mutation"
                  << " | " << std::right << std::setw(4) << "Pop"
                  << " | " << std::setw(4) << "Gen"
                  << " | " << std::setw(5) << "MutR"
                  << " | " << std::setw(5) << "XOver"
                  << " | " << std::setw(15) << "Best fitness (km)"
                  << " | " << std::setw(10) << "Total (ms)"
                  << " |\n";
        std::cout << "------|-------------------------|-------------|--------------|-----------|------|------|-------|-------|----------------|------------|\n";

        for (size_t i = 0; i < leaderboard.size(); ++i) {
            const auto& r = leaderboard[i];
            std::cout << "  " << std::right << std::setw(2) << (i + 1) << "  | "
                      << std::left << std::setw(23) << r.name << " | "
                      << std::setw(11) << r.selection << " | "
                      << std::setw(12) << r.crossover << " | "
                      << std::setw(9) << r.mutation << " | "
                      << std::right << std::setw(4) << r.popSize << " | "
                      << std::setw(4) << r.generations << " | "
                      << std::setw(5) << std::fixed << std::setprecision(2) << r.mutationRate << " | "
                      << std::setw(5) << r.crossoverRate << " | "
                      << std::setw(15) << std::setprecision(4) << r.bestDistance << " | "
                      << std::setw(10) << std::setprecision(1) << r.totalTimeMs << " |\n";
        }
        std::cout << "=========================================================================================================================================\n";
        std::cout << "Note: Leaderboard is sorted primarily by Best Distance (ascending, shortest is best).\n\n";

        return 0;
    }

    std::string toString(SelectionMethod sm) {
        switch(sm) {
            case SelectionMethod::TOURNAMENT:     return "Tournament";
            case SelectionMethod::ROULETTE_WHEEL: return "Roulette";
            case SelectionMethod::RANK_BASED:    return "Rank-Based";
            case SelectionMethod::SUS:           return "SUS";
        }
        return "Unknown";
    }

    std::string toString(CrossoverMethod cm) {
        switch(cm) {
            case CrossoverMethod::SINGLE_POINT: return "Single-Point";
            case CrossoverMethod::MULTI_POINT:  return "Multi-Point";
            case CrossoverMethod::UNIFORM:      return "Uniform";
        }
        return "Unknown";
    }

    std::string toString(MutationMethod mm) {
        switch(mm) {
            case MutationMethod::SINGLE_POINT:  return "Reset";
            case MutationMethod::SPATIAL_NUDGE:  return "Nudge";
            case MutationMethod::GAUSSIAN_CREEP: return "Gaussian";
        }
         return "Unknown";
    }
    double getFitnessEvalTime() const { return totalFitnessEvalTimeMs; }
    double getEvolutionaryTime() const { return totalEvolutionaryTimeMs; }
};


int main() {
    ThreadPool pool;
    GaEngine ga(pool);
    bool run_benchmark = false;
    bool gen_simdata = false;
    int num_experiments  = 100;
    if(run_benchmark) {
        ga.run_benchmark();
    }
    if(gen_simdata) {
        ga.run_native_random_simulation(num_experiments);
    } else {
    // ga.run_native_random_simulation(100);

    // Always configure operator enums here.
    // Numeric rates (mutationRate / crossoverRate) are only used by
    // setOperators if no config was loaded; otherwise config values win.
    ga.setOperators(SelectionMethod::TOURNAMENT, CrossoverMethod::UNIFORM,
                    MutationMethod::SPATIAL_NUDGE,
                    /*mRate=*/0.08, /*cRate=*/0.85);

    // --- Dispatch based on config when present ---
    if (GaConf::is_available()) {
        GaConf cfg = GaConf::load();

        if (cfg.homo.is_case) {
            // Parameters (pop, rates, weights, output) already applied in
            // GaEngine constructor via applyConfig(); just call the runner.
            ga.run_homogenous();
        }
        else if (cfg.hetero.is_case) {
            // Fleet built inside run_heterogenous() from cfg.hetero.radars
            ga.run_heterogenous();
        }
        else if (cfg.autoG.is_case) {
            ga.setAutoOperators(SelectionMethod::TOURNAMENT,
                                AutoCrossoverMethod::RANDOM_SUBSET,
                                /*mRate=*/0.10, /*cRate=*/0.85);
            // Fleet built inside run_auto() from cfg.autoG.radars
            // ga.run_auto();
        }
        else {
            std::cerr << "[ga] Config found but no is_case flag is set — "
                         "falling back to default homogenous run.\n";
            ga.run_homogenous(/*pathLength=*/2, /*generations=*/50,
                              /*opt=*/false, /*num_top_k=*/15);
        }

    }
    else {
        cerr << "config is not available. provide a config" << endl;
        return -1;
    }
}
    return 0;
}
