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
// CONSTANTS
// ============================================================
const int NUM_CHROMOSOMES = 100;
static constexpr double EARTH_R     = 6371000.0;
static constexpr double K_ATMOS     = 4.0 / 3.0;
static constexpr double K_EFF_R     = K_ATMOS * EARTH_R;
static constexpr double INV_2_KEFFR = 1.0 / (2.0 * K_EFF_R);
static constexpr double DEG2RAD     = 3.14159265358979 / 180.0;
static constexpr double RAD2DEG     = 180.0 / 3.14159265358979;

// ============================================================
// ENUMS
// ============================================================
enum class SelectionMethod    { TOURNAMENT, ROULETTE_WHEEL, RANK_BASED, SUS };

// Fixed-length crossover (run_homogenous / run_heterogenous)
enum class CrossoverMethod    { SINGLE_POINT, MULTI_POINT, UNIFORM };
enum class MutationMethod     { SINGLE_POINT, SPATIAL_NUDGE, GAUSSIAN_CREEP };

// Variable-length crossover (run_auto)
enum class AutoCrossoverMethod { RANDOM_SUBSET, LENGTH_NORMALIZED, POSITIONAL_PADDING };

// ============================================================
// FLEET CONFIG STRUCTS
// ============================================================

// Fixed count per type — used by run_heterogenous
struct RadarConfig {
    double rangeKm;
    int    count;
};

// Variable count per type — used by run_auto
struct FlexRadarConfig {
    double rangeKm;
    int    minCount;
    int    maxCount;
};

// ============================================================
// CHROMOSOME TYPES
// ============================================================

// Fixed-length chromosome — for run_homogenous and run_heterogenous
struct Chromosome {
    std::vector<GridPoint> radar;
    double fitness = 0.0;

    explicit Chromosome(int length) : radar(length), fitness(0.0) {}
    Chromosome() : fitness(0.0) {}
};

// Gene for variable-length chromosome — carries its own type label
// typeIndex == -1 is a NULL sentinel used only during positional-padding
// crossover; always stripped before evaluation.
struct RadarGene {
    int r, c;
    int typeIndex;
};

// Variable-length chromosome — for run_auto
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
    vector<Chromosome> population;   // used by run_homogenous / run_heterogenous
    ThreadPool& pool;
    GaConf conf;
    // --- Shared hyperparameters ---
    int    popSize       = NUM_CHROMOSOMES;
    double mutationRate  = 0.08;
    double crossoverRate = 0.85;
    double lambda        = 0.0;  // radar-count cost coefficient for run_auto

    // --- Fixed-length operator selectors ---
    SelectionMethod    currentSelection   = SelectionMethod::TOURNAMENT;
    CrossoverMethod    currentCrossover   = CrossoverMethod::UNIFORM;
    MutationMethod     currentMutation    = MutationMethod::SPATIAL_NUDGE;

    // --- Variable-length operator selectors ---
    SelectionMethod    autoSelection      = SelectionMethod::TOURNAMENT;
    AutoCrossoverMethod autoCrossover     = AutoCrossoverMethod::RANDOM_SUBSET;
    double             autoMutationRate   = 0.10;
    double             autoCrossoverRate  = 0.85;

    std::mt19937 rng;

    // =========================================================
    // SECTION 1 — FIXED-LENGTH GA (run_homogenous / run_heterogenous)
    // =========================================================

    // --- Chromosome generation ---

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

    // --- Fitness evaluation ---

    // Homogenous: all radars share the same fixed range set in radar.h defaults.
    // Fitness is normalised by total feasible cell count.
    static void evaluate_fitness(
        Global& g, const FeasibleGrid& base_grid,
        Chromosome& cn, int index, bool writePng)
    {
        std::string png_path = writePng
            ? "chromosome_" + std::to_string(index) + "_coverage.png" : "";
        auto res = calculate_radar_coverage_overlap(g, base_grid, cn.radar, 5.0, png_path);

        int    total        = (int)base_grid.cells.size();
        double unique       = (double)res[1] / total;
        double low_overlap  = (double)res[2] / total;
        double high_overlap = (double)res[3] / total;
        double spacing_gaps = (double)res[0] / total;

        // Maximise layered coverage, penalise proximity violations
        const double w1 = 0.9, w2 = 0.75, w3 = 0.9, w4 = 0.9;
        cn.fitness = (unique * w1) + (low_overlap * w2)
                   - (high_overlap * w3) - (spacing_gaps * w4);
    }

    // Heterogenous: each radar locus carries its own range from radarRanges[i].
    // Fitness is raw cell-count weighted (case2 formula).
    static void evaluate_fitness_hetero(
        Global& g, const FeasibleGrid& base_grid,
        Chromosome& cn, int index,
        const std::vector<double>& radarRanges, bool writePng)
    {
        std::string png_path = writePng
            ? "chromosome_" + std::to_string(index) + "_coverage.png" : "";
        auto res = calculate_radar_coverage_overlap(
            g, base_grid, cn.radar, 5.0, png_path, radarRanges);

        cn.fitness = (res[1] * 1.0) + (res[2] * 2.0) + (res[3] * 3.0)
                   - (res[0] * 5000.0);
    }

    // --- Parallel evaluation overloads ---

    void parallelEvaluatePopulation(
        std::vector<Chromosome>& pop, Global& globalRef,
        const FeasibleGrid& base_grid, bool writePng = false)
    {
        std::vector<std::future<void>> futures;
        futures.reserve(pop.size());
        for (size_t i = 0; i < pop.size(); i++)
            futures.push_back(pool.enqueue(
                evaluate_fitness,
                std::ref(globalRef), std::ref(base_grid),
                std::ref(pop[i]), (int)i, writePng));
        for (auto& f : futures) f.wait();
    }

    void parallelEvaluatePopulation(
        std::vector<Chromosome>& pop, Global& globalRef,
        const FeasibleGrid& base_grid,
        const std::vector<double>& radarRanges, bool writePng = false)
    {
        std::vector<std::future<void>> futures;
        futures.reserve(pop.size());
        for (size_t i = 0; i < pop.size(); i++)
            futures.push_back(pool.enqueue(
                evaluate_fitness_hetero,
                std::ref(globalRef), std::ref(base_grid),
                std::ref(pop[i]), (int)i,
                std::ref(radarRanges), writePng));
        for (auto& f : futures) f.wait();
    }

    // --- Snap to nearest feasible cell ---
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

    // --- Helpers ---

    static int64_t posKey(int r, int c) {
        return static_cast<int64_t>(r) * 1'000'000LL + c;
    }

    // Extract plain GridPoints (for PNG helpers that expect vector<GridPoint>)
    static std::vector<GridPoint> extractPositions(const std::vector<RadarGene>& genes) {
        std::vector<GridPoint> pts;
        pts.reserve(genes.size());
        for (const auto& gene : genes)
            if (gene.typeIndex >= 0) pts.push_back({gene.r, gene.c});
        return pts;
    }

    // Build a radarRanges vector aligned with extractPositions output
    // so calculate_radar_coverage_overlap (hetero path) can be reused.
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

    // --- Chromosome generation ---

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

    // --- Fitness evaluation ---

    // F(k) = single_cov×1 + dual_cov×2 + triple_cov×3
    //       - proximity_violations × 5000
    //       - λ × n          (n = total radar count)
    static void evaluate_auto_fitness(
        Global& g, const FeasibleGrid& base_grid,
        AutoChromosome& cn, int index,
        const std::vector<FlexRadarConfig>& configs,
        double lambda, bool writePng)
    {
        // Extract aligned positions + per-radar ranges, then reuse radar.h's hetero path
        auto positions = extractPositions(cn.radars);
        auto ranges    = extractRangesM(cn.radars, configs);

        std::string png_path = writePng
            ? "chromosome_" + std::to_string(index) + "_coverage.png" : "";

        auto res = calculate_radar_coverage_overlap(
            g, base_grid, positions, 5.0, png_path, ranges);

        int n = cn.totalCount();
        cn.fitness = (res[1] * 1.0) + (res[2] * 2.0) + (res[3] * 3.0)
                   - (res[0] * 5000.0)
                   - lambda * n;
    }

    void parallelEvaluateAutoPopulation(
        std::vector<AutoChromosome>& pop, Global& globalRef,
        const FeasibleGrid& base_grid,
        const std::vector<FlexRadarConfig>& configs,
        double lam, bool writePng = false)
    {
        std::vector<std::future<void>> futures;
        futures.reserve(pop.size());
        for (size_t i = 0; i < pop.size(); i++)
            futures.push_back(pool.enqueue(
                evaluate_auto_fitness,
                std::ref(globalRef), std::ref(base_grid),
                std::ref(pop[i]), (int)i,
                std::ref(configs), lam, writePng));
        for (auto& f : futures) f.wait();
    }

    // --- REPAIR ---
    // Steps: strip NULLs → deduplicate positions → trim excess → fill deficit
    void repairChromosome(
        AutoChromosome& ind, const FeasibleGrid& f,
        const std::vector<FlexRadarConfig>& configs)
    {
        // 1. Remove NULL sentinels (from positional-padding crossover)
        ind.radars.erase(
            std::remove_if(ind.radars.begin(), ind.radars.end(),
                [](const RadarGene& g) { return g.typeIndex < 0; }),
            ind.radars.end());

        // 2. Deduplicate by position (keep first occurrence per cell)
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

        // 4. Fill types that are below minCount
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
                    if (!placed) break; // grid saturated — extremely rare
                    deficit--;
                }
            }
        }
    }

    // --- MUTATION (variable-length) ---

    // ADD: insert one radar whose type is still under its maxCount
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

    // REMOVE: drop one radar whose type count is still above its minCount
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

    // MOVE: nudge one existing radar to a nearby feasible cell
    void mutate_move(AutoChromosome& ind, const FeasibleGrid& f) {
        if (ind.radars.empty()) return;

        std::uniform_int_distribution<int> distGene(0, (int)ind.radars.size()-1);
        std::uniform_int_distribution<int> distNudge(-10, 10);

        auto& gene = ind.radars[distGene(rng)];
        if (gene.typeIndex < 0) return;

        GridPoint snapped = snapToFeasible(
            gene.r + distNudge(rng),
            gene.c + distNudge(rng), f);
        gene.r = snapped.r;
        gene.c = snapped.c;
    }

    // Execute a randomly chosen mutation (add / remove / move)
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

    // --- CROSSOVER (variable-length) — method chosen at random each call ---

    // Method 1: RANDOM_SUBSET — partition each parent, swap halves
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

    // Method 2: LENGTH_NORMALIZED — sample from pooled genes to a fixed average length
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

    // Method 3: POSITIONAL_PADDING — pad shorter parent with NULLs, single-point cut
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

        // Strip NULL sentinels explicitly (repairChromosome does this too)
        auto stripNull = [](std::vector<RadarGene>& v) {
            v.erase(std::remove_if(v.begin(), v.end(),
                [](const RadarGene& g) { return g.typeIndex < 0; }), v.end());
        };
        stripNull(c1.radars);
        stripNull(c2.radars);

        return {c1, c2};
    }

    // Pick one of the three crossover methods at random each call
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
        return rankBasedSelectionAuto(pop); // RANK_BASED (SUS handled in run_auto)
    }

    // =========================================================
    // SECTION 3 — EVOLUTION ENTRY POINTS
    // =========================================================

public:
    explicit GaEngine(ThreadPool& p) : pool(p) {
        rng.seed((unsigned int)
            std::chrono::system_clock::now().time_since_epoch().count());
        std::cout << "[ga] Synchronously pre-loading all .bin matrices to memory...\n";
        g.bin.init();
    }

    // Configure fixed-length GA operators
    void setOperators(SelectionMethod sm, CrossoverMethod cm,
                      MutationMethod mm, double mRate, double cRate)
    {
        currentSelection = sm;
        currentCrossover = cm;
        currentMutation  = mm;
        mutationRate     = mRate;
        crossoverRate    = cRate;
    }

    // Configure variable-length GA operators
    void setAutoOperators(SelectionMethod sm, AutoCrossoverMethod cm,
                          double mRate, double cRate)
    {
        autoSelection    = sm;
        autoCrossover    = cm;
        autoMutationRate = mRate;
        autoCrossoverRate = cRate;
    }

    // λ = 0    → radars are free; algorithm maximises coverage only
    // λ = 10   → mild preference for fewer radars
    // λ = 1000 → strong pressure to minimise radar count
    void setLambda(double l) { lambda = l; }

    const vector<Chromosome>& get_population() const { return population; }

    // ------------------------------------------------------------------
    // run_homogenous
    // All radars share a single fixed range (set in radar.h defaults).
    // ------------------------------------------------------------------
    void run_homogenous(int pathLength, int generations, bool opt = 0, signed int num_top_k = 15) {
        population.clear();
        population.reserve(popSize);

        std::cout << "[ga] Extracting baseline Area of Operation boundaries...\n";
        auto f = Bin::getFeasibleGrid(g.bin, WRITE);

        for (int i = 0; i < popSize; i++)
            population.push_back(generate_random_chromosome(f, pathLength));

        std::cout << "[ga] Dispatching parallel coverage evaluation for Generation 0...\n";
        parallelEvaluatePopulation(population, g, f, false);

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

            parallelEvaluatePopulation(nextGen, g, f, false);
            population = std::move(nextGen);

            auto best = std::max_element(population.begin(), population.end(),
                [](const Chromosome& a, const Chromosome& b) {
                    return a.fitness < b.fitness; });
            std::cout << "  -> Gen " << gen
                      << " Best Fitness: " << best->fitness << "\n";
        }
        std::sort(population.begin(), population.end(), [](const Chromosome& a, const Chromosome& b) {
            return a.fitness > b.fitness;
        });

        if (opt == 1) {
            vector<GridPoint> top_k;
            for(int i = 0; i < num_top_k; i++)
                for(const auto& pt : population[i].radar)
                    top_k.push_back(pt);
            write_h5_indiviual(g, f, "data.h5", top_k);
        }
        auto finalBest = population[0];
        std::cout << "\n[ga] Evolution Complete. Best Lifetime Score: "
                  << finalBest.fitness << "\n";

        calculate_radar_coverage_overlap(g, f, finalBest.radar,
                                         5.0, "final_optimized_coverage.png");
        Bin::writeRadarScatterPng(finalBest.radar, "final_optimized_scatter.png");
    }

    // ------------------------------------------------------------------
    // run_heterogenous
    // Fixed fleet; each chromosome locus maps to a specific radar range.
    // ------------------------------------------------------------------
    void run_heterogenous(const std::vector<RadarConfig>& fleetConfig, int generations) {
        // Flatten fleet into a per-locus range array (metres)
        std::vector<double> flattenedRangesM;
        for (const auto& config : fleetConfig)
            for (int i = 0; i < config.count; ++i)
                flattenedRangesM.push_back(config.rangeKm * 1000.0);

        int pathLength = (int)flattenedRangesM.size();

        population.clear();
        population.reserve(popSize);

        std::cout << "[ga] Extracting baseline Area of Operation boundaries...\n";
        auto f = Bin::getFeasibleGrid(g.bin, WRITE);

        for (int i = 0; i < popSize; i++)
            population.push_back(generate_random_chromosome(f, pathLength));

        std::cout << "[ga] Dispatching parallel coverage evaluation for Generation 0...\n";
        parallelEvaluatePopulation(population, g, f, flattenedRangesM, false);

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

            parallelEvaluatePopulation(nextGen, g, f, flattenedRangesM, false);
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
                                         5.0, "final_optimized_coverage.png",
                                         flattenedRangesM);
        Bin::writeRadarScatterPng(finalBest->radar, "final_optimized_scatter.png");
    }

    // ------------------------------------------------------------------
    // run_auto
    // Variable-length, variable-count, heterogeneous fleet.
    // Each gene carries its own typeIndex; ranges are looked up from configs.
    // Uses repair, add/remove/move mutations, and 3 crossover strategies.
    // ------------------------------------------------------------------
    void run_auto(const std::vector<FlexRadarConfig>& configs, int generations) {
        // Print fleet spec
        int totalMin = 0, totalMax = 0;
        std::cout << "[ga] Fleet specification:\n";
        for (int t = 0; t < (int)configs.size(); t++) {
            std::cout << "  Type " << t << " : " << configs[t].rangeKm
                      << " km  [" << configs[t].minCount
                      << " – "    << configs[t].maxCount << " units]\n";
            totalMin += configs[t].minCount;
            totalMax += configs[t].maxCount;
        }
        std::cout << "  Total radar range : [" << totalMin
                  << ", " << totalMax << "]\n";
        std::cout << "  Lambda (cost/radar): " << lambda << "\n\n";

        std::cout << "[ga] Extracting feasible grid...\n";
        auto f = Bin::getFeasibleGrid(g.bin, WRITE);

        // --- Seed population (Generation 0) ---
        std::vector<AutoChromosome> autoPop;
        autoPop.reserve(popSize);
        for (int i = 0; i < popSize; i++)
            autoPop.push_back(generateAutoChromosome(f, configs));

        std::cout << "[ga] Evaluating Generation 0...\n";
        parallelEvaluateAutoPopulation(autoPop, g, f, configs, lambda, false);

        // --- Evolution ---
        for (int gen = 1; gen <= generations; ++gen) {
            // Sort ascending; best sits at back (easy elitism)
            std::sort(autoPop.begin(), autoPop.end(),
                [](const AutoChromosome& a, const AutoChromosome& b) {
                    return a.fitness < b.fitness; });

            std::vector<AutoChromosome> nextGen;
            nextGen.reserve(popSize);
            nextGen.push_back(autoPop.back()); // elitism

            bool useSUS = (autoSelection == SelectionMethod::SUS);
            std::vector<AutoChromosome> matingPool;
            size_t poolIdx = 0;
            if (useSUS)
                matingPool = stochasticUniversalSamplingAuto(autoPop, popSize - 1);

            while ((int)nextGen.size() < popSize) {
                AutoChromosome p1, p2;
                if (useSUS) {
                    p1 = matingPool[poolIdx % matingPool.size()]; poolIdx++;
                    p2 = matingPool[poolIdx % matingPool.size()]; poolIdx++;
                } else {
                    p1 = executeSelectionAuto(autoPop);
                    p2 = executeSelectionAuto(autoPop);
                }

                // Crossover → repair → mutate → repair
                auto [c1, c2] = executeCrossoverAuto(p1, p2);

                repairChromosome(c1, f, configs);
                repairChromosome(c2, f, configs);

                executeMutationAuto(c1, f, configs);
                executeMutationAuto(c2, f, configs);

                repairChromosome(c1, f, configs);
                repairChromosome(c2, f, configs);

                nextGen.push_back(c1);
                if ((int)nextGen.size() < popSize) nextGen.push_back(c2);
            }

            parallelEvaluateAutoPopulation(nextGen, g, f, configs, lambda, false);
            autoPop = std::move(nextGen);

            auto best = std::max_element(autoPop.begin(), autoPop.end(),
                [](const AutoChromosome& a, const AutoChromosome& b) {
                    return a.fitness < b.fitness; });

            // Per-type breakdown for the best chromosome this generation
            std::string breakdown;
            for (int t = 0; t < (int)configs.size(); t++)
                breakdown += " [T" + std::to_string(t) + ":"
                           + std::to_string(best->countType(t)) + "]";

            std::cout << "  Gen " << gen
                      << " | Fitness: " << best->fitness
                      << " | N=" << best->totalCount()
                      << breakdown << "\n";
        }

        // --- Final export ---
        auto finalBest = std::max_element(autoPop.begin(), autoPop.end(),
            [](const AutoChromosome& a, const AutoChromosome& b) {
                return a.fitness < b.fitness; });

        std::cout << "\n[ga] Evolution complete.\n";
        std::cout << "  Final fitness  : " << finalBest->fitness << "\n";
        std::cout << "  Total radars   : " << finalBest->totalCount() << "\n";
        for (int t = 0; t < (int)configs.size(); t++)
            std::cout << "  Type " << t << " (" << configs[t].rangeKm << " km) : "
                      << finalBest->countType(t) << " units\n";

        auto positions = extractPositions(finalBest->radars);
        auto ranges    = extractRangesM(finalBest->radars, configs);
        calculate_radar_coverage_overlap(g, f, positions,
                                         5.0, "final_optimized_coverage.png", ranges);
        Bin::writeRadarScatterPng(positions, "final_optimized_scatter.png");
    }
};

// ============================================================
// ENTRY POINT
// ============================================================
int main() {
    ThreadPool pool;
    GaEngine g(pool);

    // --- Example: run_homogenous ---
    g.setOperators(SelectionMethod::TOURNAMENT, CrossoverMethod::UNIFORM,
                   MutationMethod::SPATIAL_NUDGE, 0.08, 0.85);
    g.run_homogenous(2, 50, WRITE, 15);

    // --- Example: run_heterogenous ---
    // g.setOperators(SelectionMethod::TOURNAMENT, CrossoverMethod::UNIFORM,
    //                MutationMethod::SPATIAL_NUDGE, 0.08, 0.85);
    // std::vector<RadarConfig> fleet = { {50.0, 2}, {150.0, 1}, {300.0, 1} };
    // g.run_heterogenous(fleet, 50);

    // --- Example: run_auto ---
    // g.setLambda(0.0);
    // g.setAutoOperators(SelectionMethod::TOURNAMENT,
    //                    AutoCrossoverMethod::RANDOM_SUBSET,
    //                    0.10,  // mutation rate
    //                    0.85); // crossover rate
    // std::vector<FlexRadarConfig> flexFleet = {
    //     {  50.0, 2, 5 },   // short-range: 2–5 units
    //     { 300.0, 1, 3 },   // long-range:  1–3 units
    // };
    // g.run_auto(flexFleet, 5);

    return 0;
}
