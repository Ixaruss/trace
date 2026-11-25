
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <vector>
#include <deque>
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
#include <filesystem>
#include <fstream>
#include <json/json.h>


#include "../common.h"
#include "h5.h"
#include <H5Cpp.h>
#include <H5Classes.h>
#include <H5DataSet.h>
#include <H5Group.h>
#include <H5File.h>
#include "radar.h"
#include "config.h"

using namespace std;
using namespace H5;
// ============================================================
// CONSTANTS / DEFAULTS
// (all overridden at runtime when a config file is present)
// ============================================================
const int NUM_CHROMOSOMES = 100;


// Default fitness weights — used only when no config is loaded
static constexpr double DEFAULT_W1 = 0.9;   // unique coverage reward
static constexpr double DEFAULT_W2 = 0.75;  // low-overlap reward
static constexpr double DEFAULT_W3 = 0.9;   // high-overlap penalty
static constexpr double DEFAULT_W4 = 0.9;   // spacing-gap penalty
static constexpr double DEFAULT_X  = 0.0;   // lambda (radar-count cost, run_auto)

// Size of the "recent best fitness" window calculateVariance()/push_top()
// use for plateau detection. FIXED (added alongside the recent-best
// variance change): a deque with fewer than 2 entries has a mathematically
// degenerate variance of EXACTLY 0 (Welford's algorithm with one sample
// has nothing to differ from) — that used to satisfy
// `currVariance < varianceThreshold` on generation 1 every single time,
// regardless of whether the GA had actually converged. Callers now require
// `top.size() >= TOP_WINDOW` (a full window) before trusting the variance
// at all, so the convergence check can't fire on a single meaningless
// sample.
static constexpr size_t TOP_WINDOW = 10;

// ============================================================
// ENUMS
// ============================================================
enum class SelectionMethod     { TOURNAMENT, ROULETTE_WHEEL, RANK_BASED, SUS };
enum class CrossoverMethod     { SINGLE_POINT, MULTI_POINT, UNIFORM };
enum class MutationMethod      { SINGLE_POINT, SPATIAL_NUDGE, GAUSSIAN_CREEP };
enum class AutoCrossoverMethod { RANDOM_SUBSET, LENGTH_NORMALIZED, POSITIONAL_PADDING };

// Which run_* the hyperparameter sweep in run_benchmark() exercises.
enum class BenchmarkCase { HOMOGENOUS, HETEROGENOUS, AUTO };

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

// Struct representing a configuration case for the hyperparameter search.
//
// CHANGED: dropped popSize/generations/mutationRate/crossoverRate. Those
// were hardcoded PER ROW below, which looked like each profile ran with
// its own population/rate settings — but setOperators() only applies
// mRate/cRate when `!confLoaded`, and run_benchmark() requires a loaded
// GaConf to make "read conf" meaningful at all. So those four numbers were
// silently ignored on every real run; only the operator enums ever took
// effect. What actually varies between sweep profiles is the operator
// combination, so that's all this struct carries now — popSize/
// generations/mutationRate/crossoverRate are read once from conf and
// shared by every profile in the sweep (see run_benchmark()).
struct GAConfig {
    std::string name;
    SelectionMethod selection;
    CrossoverMethod crossover;
    MutationMethod mutation;
};

// run_auto's own sweep row: SelectionMethod is shared with the fixed-length
// GA, but crossover is AutoCrossoverMethod (a different enum — run_auto's
// crossover works on variable-length chromosomes), and there is no
// mutation axis at all: executeMutationAuto() always internally
// randomizes add/remove/move, it isn't a strategy the caller picks.
struct AutoGAConfig {
    std::string name;
    SelectionMethod selection;
    AutoCrossoverMethod crossover;
};

// Struct to keep track of benchmark results for the final leaderboard.
//
// CHANGED: removed totalTimeMs/evalTimeMs/evoTimeMs (timing dropped from
// the benchmark entirely) and bestDistance (this was never actually
// computed — every row hardcoded `resultDistance = 0` — replaced with
// coveragePercent, a real measured result). popSize/mutationRate/
// crossoverRate dropped too: they're now uniform across every row in the
// sweep (see GAConfig), so repeating the same three numbers on every
// leaderboard line added nothing.
struct BenchmarkResult {
    std::string name;
    std::string selection;
    std::string crossover;
    std::string mutation;
    int    generations;       // ACTUAL generations the run took (early-stop aware)
    double coveragePercent;   // % of the AOI covered by the final best solution
};

// ============================================================
// CHROMOSOME TYPES
// ============================================================

struct Chromosome {
    int id;
    std::vector<GridPoint> radar;
    double fitness = 0.0;

    explicit Chromosome(int length) : radar(length), fitness(0.0) {}
    Chromosome() : fitness(0.0) {}
    int countType(int t) const {
        return (int)std::count_if(radar.begin(), radar.end(),
            [t](const GridPoint& g) { return g.typeIndex == t; });
    }
    int totalCount() const { return (int)radar.size(); }
};

enum RunCase {
    HOMO,
    HETRO,
    AUTO
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
        else{
            threads=abs(threads*0.75);
        }

        for (size_t i = 0; i < threads; ++i) {
            workers.emplace_back([this]() {
                // Population-level parallelism already fills the cores —
                // stop computeLocalMask from spawning a nested OpenMP team
                // per worker (see rc_detail::tls_serial_rays in radar.h).
                rc_detail::tls_serial_rays = true;
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

    vector<Chromosome> population;
    ThreadPool& pool;
    GaConf conf;
    bool confLoaded = false;
    double totalFitnessEvalTimeMs = 0.0;
    double totalEvolutionaryTimeMs = 0.0;


    // --- Hyperparameters (set from config or defaults) ---
    int    popSize       = NUM_CHROMOSOMES;
    double mutationRate  = 0.08;
    double crossoverRate = 0.85;
    double lambda        = DEFAULT_X;
    Weights weights      = { DEFAULT_W1, DEFAULT_W2, DEFAULT_W3, DEFAULT_W4, DEFAULT_X };
    int cellSize;
    // Output paths (overridden by config when present)
    std::string coverageOutput = "final_optimized_coverage.png";
    std::unordered_map<uint64_t, double> cache;   // key: content hash, value: fitness
    std::mutex cacheMutex;
    // Physical spec for the homogenous evaluator (all radar share these values)
    RadarConf homoParams;

    // Rolling window of the best chromosome from each of the last
    // TOP_WINDOW generations. Used by calculateVariance() for plateau
    // detection — variance of "how much has the best fitness moved
    // recently", not variance across the current population.
    //
    // FIXED: this member is shared by run_homogenous/run_heterogenous/
    // run_auto and was never cleared between calls. run_native_random_
    // simulation() constructs ONE GaEngine and calls run_homogenous() in a
    // loop with a freshly randomized config each iteration — without a
    // clear, the first several generations of run N+1 computed variance
    // over a mix of run N's leftover chromosomes (a different fitness
    // scale entirely) and run N+1's new ones. Each run_* method now clears
    // `top` at the start, same place population.clear() already runs.
    std::deque<Chromosome> top;

    // --- Fixed-length operator selectors ---
    SelectionMethod    currentSelection  = SelectionMethod::TOURNAMENT;
    CrossoverMethod    currentCrossover  = CrossoverMethod::UNIFORM;
    MutationMethod     currentMutation   = MutationMethod::SPATIAL_NUDGE;

    // --- Variable-length operator selectors ---
    SelectionMethod     autoSelection    = SelectionMethod::TOURNAMENT;
    AutoCrossoverMethod autoCrossover    = AutoCrossoverMethod::RANDOM_SUBSET;
    double autoMutationRate  = 0.10;
    double autoCrossoverRate = 0.85;
    int cachehits = 0;
    int totalcount = 0;
    std::mt19937 rng;
    std::pair<double,double> centroid;
    std::pair<int,int> midpoint;

    // Set by run_homogenous/run_heterogenous/run_auto just before they
    // return, for run_benchmark() (a different GaEngine instance) to read
    // afterward. "Generations" is the ACTUAL count the run stopped at
    // (early-stop via variance convergence, or hit the configured max) —
    // not the configured target, which the old benchmark loop reported
    // instead and which said nothing about whether a profile converged
    // early or ran the full budget.
    int    lastGenerationsRun   = 0;
    double lastCoveragePercent  = 0.0;   // % of the AOI covered by the final best solution
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
        use_heatmap = p.use_heatmap;
        seedTopPercent = p.warm_seed_percentage;

        if (!p.output.empty())        coverageOutput = p.output;

        // Physical radar params for the homogenous run
        if (conf.homo.is_case) {
            const auto& r = conf.homo.radar;
            if (r.range        > 0) homoParams.range        = r.range * 1000.0; // km → m
            if (r.antennaHeight> 0) homoParams.antennaHeight = r.antennaHeight;
            if (r.target       > 0) homoParams.target  = r.target;
            homoParams.minAngle = r.minAngle;  // float, always valid
            homoParams.maxAngle = r.maxAngle;
        }

        if (!Config::is_available()) {
            cerr << "[ga] ERROR: applyConfig requires the environment config "
                    "(grid.resolution, area_op) — none available. Aborting "
                    "rather than continuing with an empty/zero Config.\n";
            return;
        }
        auto c = Config::load();
        cellSize = c.grid.resolution;
        if (cellSize <= 0) {
            cerr << "[ga] ERROR: applyConfig — config grid.resolution is not "
                    "set (or <= 0). Aborting.\n";
            return;
        }
        if(conf.parameters.threat_sector){
            auto x = polygonPointOnSurface(c.area_op);
              if(x.has_value()) {
                centroid = x.value();
                cout << "latlon centroid: " << x.value().first <<  " " << x.value().second << endl;
             }else {
                cout << "[ga] centroid lies outside aoi..." <<endl;
             }
              midpoint = get_midpoint_from_two_closest_points_to_sector(g, c.area_op, centroid, conf.parameters.loa, conf.parameters.roa);
        }
        std::cout << "[ga] Config loaded — "
                  << "pop=" << popSize
                  << "  mutRate=" << mutationRate
                  << "  cxRate="  << crossoverRate
                  << "  λ="       << lambda
                  << "  output=\"" << coverageOutput << "\"\n";
    }

    bool writeResultsToJson( const std::vector<GridPoint>& points, const std::string& filePath = "./output/ga/result.json") {
        std::filesystem::path targetPath(filePath);

        if (targetPath.has_parent_path()) {
            std::filesystem::create_directories(targetPath.parent_path());
        }

        Json::Value jsonArray(Json::arrayValue);
        for (const auto& point : points) {
            double x,y;
            g.bin.getLatLngFromGrid(point.r, point.c, x, y);
            jsonArray.append(x);
            jsonArray.append(y);
        }

        std::ofstream file(targetPath);
        if (!file.is_open()) {
            std::cerr << "Failed to create file at: " << filePath << "\n";
            return false;
        }

        Json::StreamWriterBuilder builder;
        builder["commentStyle"] = "None";
        builder["indentation"] = ""; // Keeps array flat on a single line

        std::unique_ptr<Json::StreamWriter> writer(builder.newStreamWriter());
        writer->write(jsonArray, &file);
        file << std::endl;

        std::cout << "[ga] Successfully created file and wrote coordinates to: " << filePath << "\n";
        return true;
    }
    static uint64_t hashChromosome(const Chromosome& cn) {
        std::vector<uint64_t> geneHashes;
        geneHashes.reserve(cn.radar.size());
        for (const auto& g : cn.radar) {
            uint64_t h = std::hash<int>()(g.r);
            h ^= std::hash<int>()(g.c)         + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
            h ^= std::hash<int>()(g.typeIndex) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
            geneHashes.push_back(h);
        }
        std::sort(geneHashes.begin(), geneHashes.end()); // order-independent
        uint64_t result = geneHashes.size();
        for (auto h : geneHashes)
            result ^= h + 0x9e3779b97f4a7c15ULL + (result << 6) + (result >> 2);
        return result;
    }
    // =========================================================
    // SECTION 1 — FIXED-LENGTH GA (run_homogenous / run_heterogenous)
    // =========================================================

    // Variance of the RECENT-BEST-FITNESS trend (the last up-to-TOP_WINDOW
    // generations' winners), not variance across the current population.
    //
    // NOTE the degenerate case this does NOT itself guard against: with
    // `candidates.size() == 1`, Welford's update sets mean=fitness on the
    // only iteration, so delta2 == 0 and the returned variance is
    // IDENTICALLY 0 regardless of what that fitness value actually is —
    // there's nothing to compute a spread over yet. That's mathematically
    // correct, not a bug in this function; the bug would be a CALLER
    // trusting that 0 as "converged" before the window has enough samples
    // to mean anything. Every call site below now requires
    // `top.size() >= TOP_WINDOW` before checking the threshold, so a
    // still-filling window can't trigger it.
    double calculateVariance(const std::deque<Chromosome>& candidates) {
          if (candidates.empty()) return 0.0;

          double mean = 0.0;
          double M2 = 0.0;
          int count = 0;

          for (const Chromosome& cn : candidates) {
              count++;
              double delta = cn.fitness - mean;
              mean += delta / count;
              double delta2 = cn.fitness - mean;
              M2 += delta * delta2;
          }

          return M2 / count; // population variance of recent best-fitness values
      }
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
        const Weights& w, const RadarConf& rp, bool writePng, vector<RadarConf> radar = {},
        const FeasibleGrid& aoi_grid = {}, bool wf = 0, double overlapping_factor = 0.25,OutputCrs opt = OutputCrs::Wgs84,
        BandLayout band = BandLayout::Rgba, bool use_curv= true, bool threat_sector=false, double loa = 0.5,
        double roa = 60, bool use_radar_model = false,int peak_power = 1, int antenna_gain=34, double wavelength = 0.1, int crossSection = 1, double detachable_threshold = 0.00000000001
        ,bool useAMSL = false, std::pair<int,int> midpoint = {}, int cellSize = 100
    )   // full AOI (allCells=WITHALL); fitness now scored against this
    {
        std::string png_path = writePng
                   ? "./output/ga/chromosome_" + std::to_string(index) + "_coverage.png" : "";

         std::vector<int> res;
        if(radar.empty()) {
                res = calculate_radar_coverage_overlap(
                   g, base_grid, cn.radar,
                   5.0,
                   png_path,
                   opt, band,
                   rp.antennaHeight,
                   rp.minAngle,
                   rp.target,
                   rp.maxAngle,
                   rp.range,
                   {},
                   aoi_grid,
                   false,
                   0,
                   use_curv,
                   threat_sector, loa, roa,
                   use_radar_model, peak_power,
                   antenna_gain, wavelength,
                   crossSection, detachable_threshold,useAMSL,midpoint);
        }else
        {
            res = calculate_radar_coverage_overlap(
                g, base_grid, cn.radar, 5.0, png_path, opt,band,0,0,0,0,0, radar, aoi_grid,
                false,
                0,
                use_curv,
                threat_sector, loa, roa,
                use_radar_model,0,0,0,0,0,useAMSL,midpoint);
        }

               // total denominator matches whatever calculate_radar_coverage_overlap scored
               // result[] against: full AOI when provided, feasible-only otherwise.
               int    total        = !aoi_grid.cells.empty()
                                    ? (int)aoi_grid.cells.size()
                                    : (int)base_grid.cells.size();
               double unique       = (double)res[1] / total;
               double low_overlap  = (double)res[2] / total;
               double high_overlap = (double)res[3] / total;
               double sector;
               if(threat_sector){
                   int digits = res[4] == 0 ? 1 : (int)log10(abs(res[4])) ;
                   int delta = 1;
                   for(int i =0; i<digits; i++) delta = delta * 10;
                   delta = delta / cellSize;
                   sector = 1.00 / ((double)res[4] / (cellSize * delta));
                   cout << "ditance  sum:" <<res[4] << " dist sum norm: "  << res[4] / (cellSize * delta);
               }

               double spacing_gaps = (double)res[0] / total;
               double overlap_factor = (low_overlap + high_overlap) / unique;
               if (wf == true){
                   // Coverage = cells seen by 1 or 2 radars. Low (2-radar) overlap is
                   // not penalised; the closeness term steers it toward the target.
                   // High (3+ radar) overlap is undesired -> explicit w4 penalty.
                   double covered = unique + low_overlap;
                   double target_error = abs(overlapping_factor - overlap_factor);
                   double closeness_reward = std::max<double>(0.0, 1.0 - target_error);
                   cn.fitness = (w.w1 * covered) + (w.w2 * closeness_reward)
                              + (threat_sector ? (w.w3 * sector) : 0);
                   (void)spacing_gaps;
               }
               else{
                   double target_error = abs(overlapping_factor - overlap_factor);
                   double closeness_reward = std::max<double>(0.0, 1.0 - target_error);
                   cn.fitness = (w.w1 * unique) + (w.w2 * closeness_reward); //+ (threat_sector ? (w.w3 * sector) : 0);
                   // double target_error = abs(0.25 - overlap_factor);
                   // double closeness_reward = std::max<double>(0.0, 1.0 - target_error);
                   // cn.fitness = (0.70 * unique) + (0.30 * closeness_reward);
               }
               if(threat_sector) std::cout << "sector :" << sector <<std::endl;
    }

    // Heterogenous fitness — weighted cell-count formula (case 2).
    // Weights struct kept for interface consistency; case-2 formula uses
    // fixed multipliers by design, but w4 is reused for the proximity penalty.


    // --- Parallel evaluation (homogenous) ---
    void parallelEvaluatePopulation(
        std::vector<Chromosome>& pop, Global& globalRef,
        const FeasibleGrid& base_grid, bool writePng = false,
        const FeasibleGrid& aoi_grid = {},OutputCrs opt = OutputCrs::Wgs84, BandLayout band = BandLayout::Rgba)
    {
        auto startTime = std::chrono::high_resolution_clock::now();

        std::vector<std::future<void>> futures;
        std::vector<size_t>   pendingIdx;
        std::vector<uint64_t> pendingKeys;
        futures.reserve(pop.size());
        pendingIdx.reserve(pop.size());
        pendingKeys.reserve(pop.size());

        std::vector<RadarConf> empty(0);

        for (size_t i = 0; i < pop.size(); i++) {
            uint64_t key = hashChromosome(pop[i]);
            totalcount++;
            bool hit = false;
            double cachedFitness = 0.0;
            {
                std::lock_guard<std::mutex> lock(cacheMutex);
                auto it = cache.find(key);
                if (it != cache.end()) { hit = true; cachedFitness = it->second; }
            }

            if (hit) {
                pop[i].fitness = cachedFitness;
                cachehits++;
                continue; // skip enqueueing — no coverage recompute needed
            }

            pendingIdx.push_back(i);
            pendingKeys.push_back(key);
            if (conf.homo.is_case) {
                futures.push_back(pool.enqueue(
                    evaluate_fitness,
                    std::ref(globalRef), std::ref(base_grid),
                    std::ref(pop[i]), (int)i,
                    weights, homoParams, writePng, empty, std::ref(aoi_grid), wf, overlapping_factor,opt,band,
                    conf.parameters.use_curv, conf.parameters.threat_sector, conf.parameters.loa, conf.parameters.roa,
                    conf.parameters.use_radar_model, conf.homo.radar.peak_power,conf.homo.radar.antenna_gain,conf.homo.radar.wavelength,
                    conf.homo.radar.crossSection, conf.homo.radar.detachable_threshold, conf.parameters.useAMSL,midpoint,cellSize));
            } else {
                futures.push_back(pool.enqueue(
                    evaluate_fitness,
                    std::ref(globalRef), std::ref(base_grid),
                    std::ref(pop[i]), (int)i,
                    weights, homoParams, writePng, std::ref(conf.hetero.radars), std::ref(aoi_grid), wf, overlapping_factor, opt,
                    band,conf.parameters.use_curv, conf.parameters.threat_sector, conf.parameters.loa, conf.parameters.roa,
                    conf.parameters.use_radar_model,0,0,0,0,0,conf.parameters.useAMSL, midpoint,cellSize));
            }
        }

        for (auto& f : futures) f.wait();

        // Populate cache with freshly-computed results
        {
            std::lock_guard<std::mutex> lock(cacheMutex);
            for (size_t k = 0; k < pendingIdx.size(); k++)
                cache[pendingKeys[k]] = pop[pendingIdx[k]].fitness;
        }

        auto endTime = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> duration = endTime - startTime;
        totalFitnessEvalTimeMs += duration.count();
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

       static std::vector<GridPoint> extractPositions(const std::vector<GridPoint>& genes) {
           std::vector<GridPoint> pts;
           pts.reserve(genes.size());
           for (const auto& gene : genes)
               if (gene.typeIndex >= 0) pts.push_back({gene.r, gene.c});
           return pts;
       }

       static std::vector<double> extractRangesM(
           const std::vector<GridPoint>& genes,
           const std::vector<FlexRadarConfig>& configs)
       {
           std::vector<double> ranges;
           ranges.reserve(genes.size());
           for (const auto& gene : genes)
               if (gene.typeIndex >= 0 && gene.typeIndex < (int)configs.size())
                   ranges.push_back(configs[gene.typeIndex].rangeKm * 1000.0);
           return ranges;
       }

       Chromosome generateChromosome(
           const FeasibleGrid& f, const std::vector<FlexRadarConfig>& configs)
       {
           Chromosome cn;
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
                           cn.radar.push_back({cell.r, cell.c, t});
                           occupied.insert(key);
                           break;
                       }
                   }
               }
           }
           return cn;
       }


    // Builds one RadarConf (quantity=1) per active gene, aligned 1:1 with
    // extractPositions()'s output order — feeds calculate_radar_coverage_overlap's
    // radarRanges param directly, so each placement gets its own type's real spec
    // (antennaHeight/minAngle/maxAngle/target/range) instead of just a range number.
    static std::vector<RadarConf> extractRadarConfigs(
        const std::vector<GridPoint>& genes,
        const std::vector<RadarConf>& typeSpecs)
    {
        std::vector<RadarConf> out;
        out.reserve(genes.size());
        for (const auto& gene : genes) {
            if (gene.typeIndex >= 0 && gene.typeIndex < (int)typeSpecs.size()) {
                RadarConf rc = typeSpecs[gene.typeIndex];
                rc.quantity = 1; // one specific unit per placement, not a pool to expand
                out.push_back(rc);
            }
        }
        return out;
    }

    // Auto fitness — same layered-count formula as heterogenous, plus λ·n penalty.
    // Weights::w4 scales the proximity penalty; weights::x is λ (radar-count cost).
    // typeSpecs = conf.autoG.radar (physical spec per type); aoi_grid = full AOI,
    // used for the result[]/fitness denominator same as evaluate_fitness.
    static void evaluate_auto_fitness(
        Global& g, const FeasibleGrid& base_grid,
        Chromosome& cn, int index,
        const std::vector<FlexRadarConfig>& configs,
        const std::vector<RadarConf>& typeSpecs,
        const Weights& w, bool writePng,
        const FeasibleGrid& aoi_grid = {}, bool wf = false,OutputCrs opt = OutputCrs::Wgs84, BandLayout band = BandLayout::Rgba
        ,bool use_curv= true, bool threat_sector=false, double loa = 0.5,
            double roa = 60, bool use_radar_model = false, double overlapping_factor = 0.25, bool useAMSL = false, std::pair<int,int> midpoint = {},int cellsize = 100)
    {
        auto positions    = extractPositions(cn.radar);
        auto radarConfigs = extractRadarConfigs(cn.radar, typeSpecs);

        std::string png_path = writePng
            ? "chromosome_" + std::to_string(index) + "_coverage.png" : "";

        auto res = calculate_radar_coverage_overlap(
            g, base_grid, positions, 5.0, png_path, opt, band,
            0, 0, 0, 0, 0, radarConfigs, aoi_grid, false, 0,  use_curv,  threat_sector,  loa,
             roa ,  use_radar_model, 0,0,0,0,0,useAMSL,midpoint, cellsize);
        // int    total = !aoi_grid.cells.empty()
        //                       ? (int)aoi_grid.cells.size()
        //                       : (int)base_grid.cells.size();
        //         const double penaltyScale = (w.w4 > 0.0) ? w.w4 * 5000.0 / DEFAULT_W4 : 5000.0;
        //         int n = cn.totalCount();
        //         cn.fitness = (res[1] * 1.0) + (res[2] * 2.0) + (res[3] * 3.0)
        //                    - (res[0] * penaltyScale)
        //                    - w.x * n;           // λ sourced from weights.x
        //         (void)total; // kept for parity with evaluate_fitness's shape; formula here is count-based, not ratio-based
        int    total        = !aoi_grid.cells.empty()
                             ? (int)aoi_grid.cells.size()
                             : (int)base_grid.cells.size();
        double unique       = (double)res[1] / total;
        double low_overlap  = (double)res[2] / total;
        double high_overlap = (double)res[3] / total;
        double spacing_gaps = (double)res[0] / total;
        double sector =  1.00 / res[4];
        double overlap_factor = (low_overlap + high_overlap) / unique;
        if (wf == true){
            // Coverage = cells seen by 1 or 2 radars. Low (2-radar) overlap is
            // not penalised; the closeness term steers it toward the target.
            // High (3+ radar) overlap is undesired -> explicit w4 penalty.
            double covered = unique + low_overlap;
            double target_error = abs(overlapping_factor - overlap_factor);
            double closeness_reward = std::max<double>(0.0, 1.0 - target_error);
            cn.fitness = (w.w1 * covered) + (w.w2 * closeness_reward)
                       + (threat_sector ? (w.w3 * sector) : 0) - w.x * cn.totalCount();
            (void)spacing_gaps;
        }
        else{
            double target_error = abs(overlapping_factor - overlap_factor);
            double closeness_reward = std::max<double>(0.0, 1.0 - target_error);
            cn.fitness = (w.w1* unique) + (w.w2 * closeness_reward) - w.x * cn.totalCount();
            // cn.fitness = (w.w1* unique) + (w.w2 * closeness_reward) +(threat_sector ? (w.w3 * sector) : 0) - w.x * cn.totalCount();
            // double target_error = abs(0.25 - overlap_factor);
            // double closeness_reward = std::max<double>(0.0, 1.0 - target_error);
            // cn.fitness = (0.70 * unique) + (0.30 * closeness_reward) - w.x * cn.totalCount();
        }
    }

    void parallelEvaluateAutoPopulation(
        std::vector<Chromosome>& pop, Global& globalRef,
        const FeasibleGrid& base_grid,
        const std::vector<FlexRadarConfig>& configs,
        const std::vector<RadarConf>& typeSpecs,
        bool writePng = false,
        const FeasibleGrid& aoi_grid = {},OutputCrs opt = OutputCrs::Wgs84, BandLayout band = BandLayout::Rgba)
    {
        std::vector<std::future<void>> futures;
        std::vector<size_t>   pendingIdx;
        std::vector<uint64_t> pendingKeys;
        futures.reserve(pop.size());
        pendingIdx.reserve(pop.size());
        pendingKeys.reserve(pop.size());

        for (size_t i = 0; i < pop.size(); i++) {
            uint64_t key = hashChromosome(pop[i]);
            totalcount++;
            bool hit = false;
            double cachedFitness = 0.0;
            {
                std::lock_guard<std::mutex> lock(cacheMutex);
                auto it = cache.find(key);
                if (it != cache.end()) { hit = true; cachedFitness = it->second; }
            }

            if (hit) {
                pop[i].fitness = cachedFitness;
                cachehits++;
                continue;
            }

            pendingIdx.push_back(i);
            pendingKeys.push_back(key);

            futures.push_back(pool.enqueue(
                evaluate_auto_fitness,
                std::ref(globalRef), std::ref(base_grid),
                std::ref(pop[i]), (int)i,
                std::ref(configs), std::ref(typeSpecs), weights, writePng,
                std::ref(aoi_grid), wf, opt,band,conf.parameters.use_curv, conf.parameters.threat_sector, conf.parameters.loa, conf.parameters.roa,
                conf.parameters.use_radar_model, overlapping_factor, conf.parameters.useAMSL, midpoint, cellSize));
        }

        for (auto& f : futures) f.wait();

        {
            std::lock_guard<std::mutex> lock(cacheMutex);
            for (size_t k = 0; k < pendingIdx.size(); k++)
                cache[pendingKeys[k]] = pop[pendingIdx[k]].fitness;
        }
    }
    // --- REPAIR ---
    void repairChromosome(
        Chromosome& ind, const FeasibleGrid& f,
        const std::vector<FlexRadarConfig>& configs)
    {
        // 1. Remove NULL sentinels
        ind.radar.erase(
            std::remove_if(ind.radar.begin(), ind.radar.end(),
                [](const GridPoint& g) { return g.typeIndex < 0; }),
            ind.radar.end());

        // 2. Deduplicate by position
        {
            std::unordered_set<int64_t> seen;
            ind.radar.erase(
                std::remove_if(ind.radar.begin(), ind.radar.end(),
                    [&](const GridPoint& g) {
                        return !seen.insert(posKey(g.r, g.c)).second;
                    }),
                ind.radar.end());
        }

        // 3. Trim types that exceed maxCount
        for (int t = 0; t < (int)configs.size(); t++) {
            int excess = ind.countType(t) - configs[t].maxCount;
            while (excess > 0) {
                std::vector<int> byType;
                for (int i = 0; i < (int)ind.radar.size(); i++)
                    if (ind.radar[i].typeIndex == t) byType.push_back(i);
                if (byType.empty()) break;
                int rem = byType[
                    std::uniform_int_distribution<int>(0, (int)byType.size()-1)(rng)];
                ind.radar.erase(ind.radar.begin() + rem);
                excess--;
            }
        }

        // 4. Fill types below minCount
        {
            std::unordered_set<int64_t> occupied;
            for (const auto& gene : ind.radar)
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
                            ind.radar.push_back({cell.r, cell.c, t});
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

    void mutate_add(Chromosome& ind, const FeasibleGrid& f,
                    const std::vector<FlexRadarConfig>& configs)
    {
        std::vector<int> eligible;
        for (int t = 0; t < (int)configs.size(); t++)
            if (ind.countType(t) < configs[t].maxCount) eligible.push_back(t);
        if (eligible.empty()) return;

        int typeIdx = eligible[
            std::uniform_int_distribution<int>(0, (int)eligible.size()-1)(rng)];

        std::unordered_set<int64_t> occupied;
        for (const auto& gene : ind.radar)
            if (gene.typeIndex >= 0) occupied.insert(posKey(gene.r, gene.c));

        std::uniform_int_distribution<int> distF(0, (int)f.cells.size()-1);
        for (int attempt = 0; attempt < 300; attempt++) {
            const auto& cell = f.cells[distF(rng)];
            if (!occupied.count(posKey(cell.r, cell.c))) {
                ind.radar.push_back({cell.r, cell.c, typeIdx});
                return;
            }
        }
    }

    void mutate_remove(Chromosome& ind,
                       const std::vector<FlexRadarConfig>& configs)
    {
        std::vector<int> removable;
        for (int i = 0; i < (int)ind.radar.size(); i++) {
            int t = ind.radar[i].typeIndex;
            if (t >= 0 && ind.countType(t) > configs[t].minCount)
                removable.push_back(i);
        }
        if (removable.empty()) return;

        int idx = removable[
            std::uniform_int_distribution<int>(0, (int)removable.size()-1)(rng)];
        ind.radar.erase(ind.radar.begin() + idx);
    }

    void mutate_move(Chromosome& ind, const FeasibleGrid& f) {
        if (ind.radar.empty()) return;

        auto& gene = ind.radar[
            std::uniform_int_distribution<int>(0, (int)ind.radar.size()-1)(rng)];
        if (gene.typeIndex < 0) return;

        // Multi-scale step: log-uniform sigma from 5 cells (~0.5 km) up to
        // 600 cells (~54 km @ 90 m). Small steps polish a layout, large steps
        // let a radar jump into a different gap. The old fixed +/-10 cell
        // nudge (~0.9 km) could never move a 85-100 km radar meaningfully.
        double logLo = std::log(5.0), logHi = std::log(600.0);
        double sigma = std::exp(std::uniform_real_distribution<double>(logLo, logHi)(rng));
        std::normal_distribution<double> step(0.0, sigma);
        GridPoint snapped = snapToFeasible(
            gene.r + (int)std::lround(step(rng)),
            gene.c + (int)std::lround(step(rng)), f);
        gene.r = snapped.r;
        gene.c = snapped.c;
    }

    void executeMutationAuto(Chromosome& ind, const FeasibleGrid& f,
                             const std::vector<FlexRadarConfig>& configs)
    {
        if (std::uniform_real_distribution<double>(0.0, 1.0)(rng) > autoMutationRate)
            return;
        // Weighted: mostly relocate (that is what refines a layout); add /
        // remove only occasionally change the fleet size.
        double u = std::uniform_real_distribution<double>(0.0, 1.0)(rng);
        if      (u < 0.15) mutate_add(ind, f, configs);
        else if (u < 0.30) mutate_remove(ind, configs);
        else               mutate_move(ind, f);
    }

    // --- CROSSOVER (variable-length) ---

    std::pair<Chromosome, Chromosome> crossover_randomSubset(
        const Chromosome& p1, const Chromosome& p2)
    {
        auto splitParent = [&](const Chromosome& parent)
            -> std::pair<std::vector<GridPoint>, std::vector<GridPoint>>
        {
            int n = (int)parent.radar.size();
            int k = (n > 0) ? std::uniform_int_distribution<int>(0, n)(rng) : 0;
            std::vector<int> idx(n);
            std::iota(idx.begin(), idx.end(), 0);
            std::shuffle(idx.begin(), idx.end(), rng);
            std::vector<GridPoint> chosen, rest;
            for (int i = 0; i < n; i++)
                (i < k ? chosen : rest).push_back(parent.radar[idx[i]]);
            return {chosen, rest};
        };

        auto [a_chosen, a_rest] = splitParent(p1);
        auto [b_chosen, b_rest] = splitParent(p2);

        Chromosome c1, c2;
        c1.radar.insert(c1.radar.end(), a_chosen.begin(), a_chosen.end());
        c1.radar.insert(c1.radar.end(), b_chosen.begin(), b_chosen.end());
        c2.radar.insert(c2.radar.end(), a_rest.begin(),   a_rest.end());
        c2.radar.insert(c2.radar.end(), b_rest.begin(),   b_rest.end());
        return {c1, c2};
    }

    std::pair<Chromosome, Chromosome> crossover_lengthNormalized(
        const Chromosome& p1, const Chromosome& p2)
    {
        int avgLen = ((int)p1.radar.size() + (int)p2.radar.size() + 1) / 2;
        std::vector<GridPoint> pool;
        pool.insert(pool.end(), p1.radar.begin(), p1.radar.end());
        pool.insert(pool.end(), p2.radar.begin(), p2.radar.end());
        int take = std::min(avgLen, (int)pool.size());

        auto buildChild = [&]() {
            std::shuffle(pool.begin(), pool.end(), rng);
            Chromosome child;
            child.radar.assign(pool.begin(), pool.begin() + take);
            return child;
        };
        return {buildChild(), buildChild()};
    }

    std::pair<Chromosome, Chromosome> crossover_positionalPadding(
        const Chromosome& p1, const Chromosome& p2)
    {
        int maxLen = (int)std::max(p1.radar.size(), p2.radar.size());
        if (maxLen == 0) return {p1, p2};

        const GridPoint null_gene = {0, 0, -1};
        std::vector<GridPoint> pad1(p1.radar), pad2(p2.radar);
        while ((int)pad1.size() < maxLen) pad1.push_back(null_gene);
        while ((int)pad2.size() < maxLen) pad2.push_back(null_gene);

        int cut = std::uniform_int_distribution<int>(1, maxLen)(rng);
        Chromosome c1, c2;
        for (int i = 0; i < maxLen; i++) {
            if (i < cut) { c1.radar.push_back(pad1[i]); c2.radar.push_back(pad2[i]); }
            else          { c1.radar.push_back(pad2[i]); c2.radar.push_back(pad1[i]); }
        }

        auto stripNull = [](std::vector<GridPoint>& v) {
            v.erase(std::remove_if(v.begin(), v.end(),
                [](const GridPoint& g) { return g.typeIndex < 0; }), v.end());
        };
        stripNull(c1.radar);
        stripNull(c2.radar);
        return {c1, c2};
    }

    std::pair<Chromosome, Chromosome> executeCrossoverAuto(
        const Chromosome& p1, const Chromosome& p2)
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

    Chromosome tournamentSelectionAuto(const std::vector<Chromosome>& pop) {
        std::uniform_int_distribution<int> dist(0, (int)pop.size()-1);
        Chromosome best = pop[dist(rng)];
        for (int i = 0; i < 2; ++i) {
            const Chromosome& comp = pop[dist(rng)];
            if (comp.fitness > best.fitness) best = comp;
        }
        return best;
    }

    Chromosome rouletteWheelSelectionAuto(const std::vector<Chromosome>& pop) {
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

    Chromosome rankBasedSelectionAuto(const std::vector<Chromosome>& sorted_pop) {
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

    std::vector<Chromosome> stochasticUniversalSamplingAuto(
        const std::vector<Chromosome>& pop, int numToSelect)
    {
        double minFit = pop[0].fitness;
        for (const auto& ind : pop) if (ind.fitness < minFit) minFit = ind.fitness;
        double offset = (minFit < 0) ? std::abs(minFit) + 1.0 : 0.0;

        double total = 0.0;
        for (const auto& ind : pop) total += (ind.fitness + offset);

        double ptr_dist  = total / numToSelect;
        double start_ptr = std::uniform_real_distribution<double>(0.0, ptr_dist)(rng);

        std::vector<Chromosome> selected;
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

    Chromosome executeSelectionAuto(const std::vector<Chromosome>& pop) {
        if (autoSelection == SelectionMethod::TOURNAMENT)     return tournamentSelectionAuto(pop);
        if (autoSelection == SelectionMethod::ROULETTE_WHEEL) return rouletteWheelSelectionAuto(pop);
        return rankBasedSelectionAuto(pop);
    }

    // Pushes ch as the newest "best of this generation" sample, evicting
    // the oldest once the window exceeds TOP_WINDOW entries.
    void push_top(std::deque<Chromosome>& pop, const Chromosome& ch) {
        pop.push_back(ch);
        if (pop.size() > TOP_WINDOW) {
            pop.pop_front();
        }
    }

    // =========================================================
    // SECTION 3 — EVOLUTION ENTRY POINTS
    // =========================================================

public:
    bool wf = false;   // true = coverage-based overlap-target fitness (testing)
    double overlapping_factor = conf.parameters.overlapFactor;
    Global g;
    bool use_curv = true;
    bool use_heatmap = false;
    double seedTopPercent = 80;

    void reloadConfig() {
        if (GaConf::is_available()) {
            conf = GaConf::load();
            confLoaded = true;
            applyConfig();
            if (!Config::is_available()) {
                cerr << "[ga] ERROR: reloadConfig requires the environment "
                        "config (area_op, fcase flags) — none available. "
                        "g.bin's flags/area_op were NOT refreshed.\n";
                return;
            }
            {
                // CHANGED: a separate `g.bin.conf.area_op = Config::load().area_op;`
                // (its own Config::load() call) used to run right before the
                // line below and was immediately overwritten by it —
                // g.bin.conf = Config::load() already brings area_op along
                // as part of the whole struct. Removed the redundant call;
                // behavior is identical, one fewer file parse per reload.
                g.bin.conf = Config::load();
                g.bin.flag_ib = g.bin.conf.fcase.ib;
                g.bin.flag_rail = g.bin.conf.fcase.rails;
                g.bin.flag_road = g.bin.conf.fcase.roads;
                g.bin.flag_water_area = g.bin.conf.fcase.waterbodies;
                g.bin.flag_water_line = g.bin.conf.fcase.waterlines;
                g.bin.flag_slope = g.bin.conf.fcase.slope;
            }
        }
    }
    // A GA over zero feasible cells can't do anything useful:
    // generate_random_chromosome() builds uniform_int_distribution(0, -1),
    // which is undefined behaviour, so the run hangs or crashes with no
    // explanation. Call this right after getFeasibleGrid() and bail out with
    // the reason instead. Works the same whether the constraints matrix was
    // loaded from map_tiles/ (block_size > 0) or a single map.bin (== 0).
    bool feasibleGridUsable(const FeasibleGrid& f, const char* who) const {
        if (!f.cells.empty()) return true;
        std::cerr << "[ga] " << who << ": feasible grid is EMPTY — aborting this run.\n"
                  << "  constraints loaded : " << (g.bin.vector_loaded    ? "yes" : "NO") << "\n"
                  << "  elevation loaded   : " << (g.bin.elevation_loaded ? "yes" : "NO") << "\n"
                  << "  window valid       : " << (g.bin.window_valid     ? "yes" : "NO") << "\n"
                  << "  Look above for the [bin] / [tiles] / [feasible] lines that say why.\n";
        return false;
    }
    explicit GaEngine(ThreadPool& p) : pool(p) {
        rng.seed((unsigned int)
            std::chrono::system_clock::now().time_since_epoch().count());

        // Load config if available; otherwise the defaults declared above are used
        if (GaConf::is_available()) {
            conf = GaConf::load();
            use_curv = conf.parameters.use_curv;
            confLoaded = true;
            applyConfig();
        } else {
            std::cout << "[ga] No config file found — exiting\n";
            return;
        }

        std::cout << "[ga] Synchronously pre-loading all .bin matrices to memory...\n";
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


        std::cout << "--> Starting " << num_experiments << " Native Randomized Homogeneous Simulations..." << std::endl;
        ThreadPool pool;
        GaEngine ga(pool);

        ga.setOperators(SelectionMethod::TOURNAMENT, CrossoverMethod::UNIFORM,
                        MutationMethod::SPATIAL_NUDGE,
                        /*mRate=*/0.08, /*cRate=*/0.85);

        for (int i = 1; i <= num_experiments; ++i) {
            std::cout << "\n==================================================" << std::endl;
            std::cout << " RUN SIMULATION [" << i << " / " << num_experiments << "]" << std::endl;
            std::cout << "==================================================" << std::endl;

            ga.reloadConfig();
            ga.run_homogenous();
            std::cout << ">> Run [" << i << "] completed and saved successfully." << std::endl;
        }
    }
    int seed_population_from_heatmap(std::vector<Chromosome>& pop, const FeasibleGrid& f, int length) {

               string dataFilePath = "./output/heatmap/data.h5";
               bool matchRadarParams = false;
               if (seedTopPercent <= 0.0 || length <= 0 || popSize <= 0 || f.cells.empty()) return 0;

               std::vector<h5_lib::RunData> runs;
               try { runs = read_h5_runs(dataFilePath); }
               catch (const std::exception& e) {
                   std::cout << "[ga] could not read '" << dataFilePath << "': " << e.what() << "\n";
                   return 0;
               }
               if (runs.empty()) {
                   std::cout << "[ga] no runs found in '" << dataFilePath << "' - starting fully random\n";
                   return 0;
               }

               const size_t K = static_cast<size_t>(
                   std::ceil(popSize * std::min(seedTopPercent, 100.0) / 100.0));

               auto sites = collect_seed_sites(runs, matchRadarParams,
                                               homoParams.antennaHeight, homoParams.target,
                                               homoParams.minAngle, homoParams.maxAngle,
                                               homoParams.range);
               // Only the best few times K*length sites can end up used; trim before
               // the feasibility scan.
               const size_t cap = std::max<size_t>(K * static_cast<size_t>(length) * 10, 1000);
               if (sites.size() > cap) sites.resize(cap);
               std::cout<< "site in ga1: " << sites.size() <<std::endl;
               // Keep only sites that are feasible cells of THIS run's grid (the file
               // may come from another AOI or constraint set). One pass over f.cells.
               std::unordered_set<int64_t> wanted, feasible;
               wanted.reserve(sites.size() * 2);
               for (const auto& s : sites) wanted.insert(posKey(s.r, s.c));
               for (const auto& cell : f.cells) {
                   const int64_t key = posKey(cell.r, cell.c);
                   if (wanted.count(key)) {
                       feasible.insert(key);
                       if (feasible.size() == wanted.size()) break;
                   }
               }
               // sites.erase(std::remove_if(sites.begin(), sites.end(),
               //                 [&](const h5_lib::SeedSite& s) { return !feasible.count(posKey(s.r, s.c)); }),
               //             sites.end());
               std::cout<< "site in ga2: " << sites.size() <<std::endl;
               auto groups = deal_seed_chromosomes(sites, K, static_cast<size_t>(length),0.3,rng);


               for (const auto& grp : groups) {
                   Chromosome cn(length);
                   for (int i = 0; i < length; ++i) cn.radar[i] = {grp[i].r, grp[i].c, 0};
                   pop.push_back(std::move(cn));
               }
               std::cout << "[ga]" << groups.size() << " of " << K
                         << " seed chromosomes from '" << dataFilePath << "' ("
                         << sites.size() << " usable sites in "
                         << runs.size() << " run(s))\n";
               return static_cast<int>(groups.size());
           }
    // ------------------------------------------------------------------
    // run_homogenous
    // All radar share a single fixed range (set in radar.h defaults).
    // chromosomeLength and generations sourced from config when available.
    // ------------------------------------------------------------------
    void run_homogenous(int reduction = 50, bool write_png = false, OutputCrs opts = OutputCrs::Wgs84, BandLayout band = BandLayout::Rgba, int pathLength   = 2,
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

        const double marginKm = homoParams.range / 1000.0;   // range is metres here
        g.load({ CONSTRAINTS, ELEVATION }, AOI, 0.0, 0.0, marginKm);
        // if (threat_sector_flag) {
        //     auto ax = computeBearing(g, centroid.first, centroid.second, cell.r, cell.c);
        //     // std::cout << "bearing angle for cell :" << ax << " loa: " << loa << " roa: " << roa << std::endl;
        //     if(ax > loa && ax < roa) {
        //         cell.overlap = 9999;
        //         // std::cout << "true \n";
        //     }// }else std::cout << "false \n";

        // }
        population.clear();
        population.reserve(popSize);

        top.clear();

        lastGenerationsRun = 0;
        lastCoveragePercent = 0.0;

        std::cout << "[ga] [homogenous] chromosomes=" << pathLength
                  << "  generations=" << generations
                  << "  popSize=" << popSize << "\n";
        std::cout << "[ga] Extracting baseline Area of Operation boundaries...\n";
        auto f = Bin::getFeasibleGrid(g.bin, 0);
        auto aoiGrid = Bin::getFeasibleGrid(g.bin, 0, 0, WITHALL); // unfiltered AOI, same AABB, for PNG background
        if (!feasibleGridUsable(f, "run_homogenous")) return;
         auto initStart = std::chrono::high_resolution_clock::now();

         if(use_heatmap)
            seed_population_from_heatmap(population, f, pathLength);

        for (int i = 0; i < popSize; i++)
            population.push_back(generate_random_chromosome(f, pathLength));

        auto initEnd = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> initDuration = initEnd - initStart;

        std::cout << "[ga] Dispatching parallel coverage evaluation for Generation 0...\n";
        parallelEvaluatePopulation(population, g, f, false, aoiGrid);

        std::normal_distribution<double> gaussDist(0.0, 10.0);
        int gen = 1;
        int count = 0;
        while (true) {

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
            push_top(top, *best);
            double currVariance = calculateVariance(top);
            std::cout << "  -> Gen " << gen
                      << " Best Fitness: " << best->fitness << "\n";
            std::cout << "  -> Variance: " << currVariance << "\n";
            if (opt) {
                vector<GridPoint> top_k;
                for (int i = 0; i < population.size() /** conf.parameters.saveSolutions */ && i < (int)population.size(); i++)
                    for (const auto& pt : population[i].radar)
                        top_k.push_back(pt);
                write_h5_indiviual(g, f, "./output/heatmap/data.h5", top_k,homoParams.antennaHeight,
                homoParams.target,
                homoParams.minAngle,
                homoParams.maxAngle,
                homoParams.range);
            }

            if (top.size() >= TOP_WINDOW && currVariance < conf.parameters.varianceThreshold) {
                if (++count >= conf.parameters.forgiveFactor) break;
            } else {
                count = 0;
            }
            if(gen >= generations) break;
            ++gen;
        }

        std::sort(population.begin(), population.end(),
            [](const Chromosome& a, const Chromosome& b) {
                return a.fitness > b.fitness; });



        auto& finalBest = population[0];
        std::cout << "\n[ga] Evolution Complete. Best Lifetime Score: "
                  << finalBest.fitness << "\n";

        std::cout << "[ga] Cache hits: " << cachehits << " / " << totalcount << "\n";
        cachehits = 0;
        totalcount = 0;

        lastGenerationsRun = gen;
        auto finalRes = calculate_radar_coverage_overlap(
            g, f, finalBest.radar,
            5.0, coverageOutput, opts, band,
            homoParams.antennaHeight,
            homoParams.minAngle,
            homoParams.target,
            homoParams.maxAngle,
            homoParams.range,
            {},
            aoiGrid, write_png, reduction,

            conf.parameters.use_curv, conf.parameters.threat_sector,
            conf.parameters.loa, conf.parameters.roa,
            conf.parameters.use_radar_model, conf.homo.radar.peak_power,
            conf.homo.radar.antenna_gain, conf.homo.radar.wavelength,
            conf.homo.radar.crossSection, conf.homo.radar.detachable_threshold,
            conf.parameters.useAMSL, midpoint, cellSize);
        {
            const double total = !aoiGrid.cells.empty()
                                ? (double)aoiGrid.cells.size()
                                : (double)f.cells.size();
            lastCoveragePercent = total > 0
                ? 100.0 * (finalRes[1] + finalRes[2] + finalRes[3]) / total
                : 0.0;
        }
         std::cout << "[ga] total percentage covered: " << ((finalRes[1] + finalRes[2]) / aoiGrid.cells.size()) *100 << "% \n";

        writeResultsToJson(finalBest.radar);
    }

    // ------------------------------------------------------------------
    // run_heterogenous
    // Fixed fleet; each chromosome locus maps to a specific radar range.
    // Fleet sourced from config (hetero.radar) when available.
    // ------------------------------------------------------------------
    void run_heterogenous(int reduction = 50, bool write_png = false,OutputCrs opts = OutputCrs::Wgs84, BandLayout band = BandLayout::Rgba, std::vector<RadarConfig> fleetConfig = {},
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
            std::cerr << "[ga] run_heterogenous: no fleet config supplied and"
                         "config unavailable — aborting.\n";
            return ;
        }

        // Window margin = the LARGEST range in the fleet, so even the
        // longest-reaching radar has terrain out to its full radius. See
        // run_homogenous for why. Done after the fleet is built so the max
        // is actually known.
        double marginKm = 0.0;
        for (const auto& config : fleetConfig)
            marginKm = std::max(marginKm, config.rangeKm);
        g.load({ CONSTRAINTS, ELEVATION }, AOI, 0.0, 0.0, marginKm);
        // FIXED: same missing-clear issue as run_homogenous — see the
        // comment on the `top` member declaration.
        top.clear();
        lastGenerationsRun = 0;
        lastCoveragePercent = 0.0;

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
        auto f = Bin::getFeasibleGrid(g.bin, 0);
        auto aoiGrid = Bin::getFeasibleGrid(g.bin, 0, 0, WITHALL); // unfiltered AOI, same AABB, for PNG background
        if (!feasibleGridUsable(f, "run_heterogenous")) return;

        for (int i = 0; i < popSize; i++)
            population.push_back(generate_random_chromosome(f, pathLength));

        std::cout << "[ga] Dispatching parallel coverage evaluation for Generation 0...\n";
        parallelEvaluatePopulation(population, g, f, false, aoiGrid);

        std::normal_distribution<double> gaussDist(0.0, 10.0);

        int gen = 1;
        int count = 0;
        while (true) {
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
            push_top(top, *best);
            double currVariance = calculateVariance(top);
            std::cout << "  -> Gen " << gen
                      << " Best Fitness: " << best->fitness << "\n";
                      std::cout << "  -> Variance: " << currVariance << "\n";
            if (top.size() >= TOP_WINDOW && currVariance < conf.parameters.varianceThreshold) {
                if (++count >= conf.parameters.forgiveFactor) break;
            } else {
                count = 0;
            }
            if(gen >= generations) break;
            ++gen;
        }

        auto finalBest = std::max_element(population.begin(), population.end(),
            [](const Chromosome& a, const Chromosome& b) {
                return a.fitness < b.fitness; });
        std::cout << "\n[ga] Evolution Complete. Best Lifetime Score: "
                  << finalBest->fitness << "\n";
        std::cout << "[ga] Cache hits: " << cachehits << " / " << totalcount << "\n";
        cachehits = 0;
        totalcount = 0;

        lastGenerationsRun = gen;
        auto finalRes = calculate_radar_coverage_overlap(g, f, finalBest->radar,
                                         5.0, coverageOutput,opts, band,0,0,0,0,0, conf.hetero.radars,
                                         aoiGrid,write_png,reduction,
                                         // FIX: same model as the fitness evaluations
                                         conf.parameters.use_curv, conf.parameters.threat_sector,
                                         conf.parameters.loa, conf.parameters.roa,
                                         conf.parameters.use_radar_model, 0,0,0,0,0,
                                         conf.parameters.useAMSL, midpoint, cellSize);
        {
            const double total = !aoiGrid.cells.empty()
                                ? (double)aoiGrid.cells.size()
                                : (double)f.cells.size();
            lastCoveragePercent = total > 0
                ? 100.0 * (finalRes[1] + finalRes[2] + finalRes[3]) / total
                : 0.0;
        }

        writeResultsToJson(finalBest->radar);
    }

    // ------------------------------------------------------------------
    // run_auto
    // Variable-length, variable-count, heterogeneous fleet.
    // Fleet sourced from config (autoG.radar) when available.
    // ------------------------------------------------------------------
    void run_auto(int reduction = 50, bool write_png = false,OutputCrs opts = OutputCrs::Wgs84, BandLayout band = BandLayout::Rgba, std::vector<FlexRadarConfig> configs = {},
                      int generations = 50)
        {
            if (confLoaded && conf.autoG.is_case) {
                if (conf.parameters.numGenerations > 0)
                    generations = conf.parameters.numGenerations;

                // Build flex fleet from config if none supplied by caller.
                // Each RadarConf.quantity is used as maxCount; minCount is
                // clamped to at least 0 and at most quantity.
                if (configs.empty()) {
                    int numTypes = (int)conf.autoG.radars.size();
                    for (const auto& rc : conf.autoG.radars) {
                        int minC = (numTypes > 0)
                            ? std::max(0, conf.autoG.minRadars / numTypes)
                            : 0;
                        int maxC = (numTypes > 0)
                            ? std::max(minC, conf.autoG.maxRadars / numTypes)
                            : minC;
                        configs.push_back({ (double)rc.range, minC, maxC });
                    }
                }
            }

            if (configs.empty()) {
                std::cerr << "[ga] run_auto: no fleet config supplied and "
                             "config unavailable — aborting.\n";
                return;
            }
            // Window margin = the LARGEST type range in the fleet. See
            // run_homogenous for why the AOI bbox alone is not enough.
            double marginKm = 0.0;
            for (const auto& c : configs)
                marginKm = std::max(marginKm, c.rangeKm);
            g.load({ CONSTRAINTS, ELEVATION }, AOI, 0.0, 0.0, marginKm);
            // FIXED: same missing-clear issue as run_homogenous — see the
            // comment on the `top` member declaration.
            top.clear();
            lastGenerationsRun = 0;
            lastCoveragePercent = 0.0;
            // typeSpecs carries the FULL physical spec per type (antennaHeight, angles,
            // target, range) — configs (FlexRadarConfig) only carries min/max counts for
            // chromosome generation; typeSpecs is what feeds actual coverage math.
            const std::vector<RadarConf>& typeSpecs = conf.autoG.radars;

            // Print fleet spec
            int totalMin = 0, totalMax = 0;
            std::cout << "[ga] [auto] Fleet specification:\n";
            for (int t = 0; t < (int)configs.size(); t++) {
                std::cout << "  Type " << t << " : " << configs[t].rangeKm
                          << " km  [" << configs[t].minCount
                          << " – "    << configs[t].maxCount << " units]\n";
                totalMin += configs[t].minCount;
                totalMax += configs[t].maxCount;
            }
            std::cout << "  Total radar range : [" << totalMin << ", " << totalMax << "]\n";
            std::cout << "  Lambda (cost/radar): " << weights.x << "\n";
            std::cout << "  generations=" << generations
                      << "  popSize=" << popSize << "\n\n";

            std::cout << "[ga] Extracting feasible grid...\n";
            auto f = Bin::getFeasibleGrid(g.bin, 0);
            auto aoiGrid = Bin::getFeasibleGrid(g.bin, 0, 0, WITHALL); // unfiltered AOI, same AABB, for coverage + PNG background
            if (!feasibleGridUsable(f, "run_auto")) return;

            // Seed population
            std::vector<Chromosome> autoPop;
            autoPop.reserve(popSize);
            for (int i = 0; i < popSize; i++)
                autoPop.push_back(generateChromosome(f, configs));

            std::cout << "[ga] Evaluating Generation 0...\n";
            parallelEvaluateAutoPopulation(autoPop, g, f, configs, typeSpecs, false, aoiGrid);
            int gen = 1;
            int count = 0;
            while(true) {
                std::sort(autoPop.begin(), autoPop.end(),
                    [](const Chromosome& a, const Chromosome& b) {
                        return a.fitness < b.fitness; });

                std::vector<Chromosome> nextGen;
                nextGen.reserve(popSize);
                nextGen.push_back(autoPop.back()); // elitism

                bool useSUS = (autoSelection == SelectionMethod::SUS);
                std::vector<Chromosome> matingPool;
                size_t poolIdx = 0;
                if (useSUS)
                    matingPool = stochasticUniversalSamplingAuto(autoPop, popSize - 1);

                while ((int)nextGen.size() < popSize) {
                    Chromosome p1, p2;
                    if (useSUS) {
                        p1 = matingPool[poolIdx % matingPool.size()]; poolIdx++;
                        p2 = matingPool[poolIdx % matingPool.size()]; poolIdx++;
                    } else {
                        p1 = executeSelectionAuto(autoPop);
                        p2 = executeSelectionAuto(autoPop);
                    }

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

                parallelEvaluateAutoPopulation(nextGen, g, f, configs, typeSpecs, false, aoiGrid);
                autoPop = std::move(nextGen);

                auto best = std::max_element(autoPop.begin(), autoPop.end(),
                    [](const Chromosome& a, const Chromosome& b) {
                        return a.fitness < b.fitness; });

                std::string breakdown;
                for (int t = 0; t < (int)configs.size(); t++)
                    breakdown += " [T" + std::to_string(t) + ":"
                               + std::to_string(best->countType(t)) + "]";

                push_top(top, *best);
                double currVariance = calculateVariance(top);

                std::cout << "  Gen " << gen
                          << " | Fitness: " << best->fitness
                          << " | N=" << best->totalCount()
                          << breakdown << "\n";
                          std::cout << "  -> Variance: " << currVariance << "\n";

                if (top.size() >= TOP_WINDOW && currVariance < conf.parameters.varianceThreshold) {
                    if (++count >= conf.parameters.forgiveFactor) break;
                } else {
                    count = 0;
                }
                if(gen >= generations) break;
                ++gen;
            }

            auto finalBest = std::max_element(autoPop.begin(), autoPop.end(),
                [](const Chromosome& a, const Chromosome& b) {
                    return a.fitness < b.fitness; });

            std::cout << "\n[ga] Evolution complete.\n";
            std::cout << "  Final fitness  : " << finalBest->fitness << "\n";
            std::cout << "  Total radar   : " << finalBest->totalCount() << "\n";
            for (int t = 0; t < (int)configs.size(); t++)
                std::cout << "  Type " << t << " (" << configs[t].rangeKm << " km) : "
                          << finalBest->countType(t) << " units\n";
            std::cout << "  Cache hits: " << cachehits << " / " << totalcount << "\n";
            cachehits = 0;
            totalcount = 0;
            auto positions    = extractPositions(finalBest->radar);
            auto radarConfigs = extractRadarConfigs(finalBest->radar, typeSpecs);

            lastGenerationsRun = gen;
            auto finalRes = calculate_radar_coverage_overlap(g, f, positions,
                                             5.0, coverageOutput, opts, band,
                                             0, 0, 0, 0, 0, radarConfigs,
                                             aoiGrid,write_png,reduction,
                                             // FIX: same model as evaluate_auto_fitness
                                             conf.parameters.use_curv, conf.parameters.threat_sector,
                                             conf.parameters.loa, conf.parameters.roa,
                                             conf.parameters.use_radar_model, 0,0,0,0,0,
                                             conf.parameters.useAMSL, midpoint, cellSize);
            {
                const double total = !aoiGrid.cells.empty()
                                    ? (double)aoiGrid.cells.size()
                                    : (double)f.cells.size();
                lastCoveragePercent = total > 0
                    ? 100.0 * (finalRes[1] + finalRes[2] + finalRes[3]) / total
                    : 0.0;
            }

            writeResultsToJson(finalBest->radar);
        }
    // which: BenchmarkCase::HOMOGENOUS / HETEROGENOUS / AUTO — selects which
    // of the three run_* methods every profile in the sweep is run through.
    //
    // CHANGED (full rewrite):
    //  - Takes `which` instead of only ever calling run_homogenous().
    //  - sweepConfigs no longer carries its own popSize/generations/
    //    mutationRate/crossoverRate per row (see GAConfig) — those four
    //    numbers are read ONCE from conf, here, and shared by every
    //    profile. They used to be hardcoded per row and looked
    //    configurable, but setOperators()/setAutoOperators() only apply
    //    mRate/cRate when `!confLoaded`, which is never true once a
    //    GaConf has loaded (the normal case for a benchmark run) — so
    //    those per-row numbers were dead weight, always overridden by
    //    conf regardless of what the row said.
    //  - Timing removed entirely (no start/end/elapsed, no
    //    evalTimeMs/evoTimeMs/totalTimeMs).
    //  - Leaderboard reports ACTUAL generations run (early-stop aware,
    //    via lastGenerationsRun) and real coverage% (via
    //    lastCoveragePercent), sorted by coverage% descending — replacing
    //    bestDistance, which was never actually computed (every row
    //    hardcoded `resultDistance = 0`).
    // which: BenchmarkCase::HOMOGENOUS / HETEROGENOUS / AUTO — selects which
    // of the three run_* methods the sweep is run through. HOMOGENOUS and
    // HETEROGENOUS share one 16-row operator table (they use the same
    // fixed-length SelectionMethod/CrossoverMethod/MutationMethod operator
    // set); AUTO uses its own table (SelectionMethod x AutoCrossoverMethod
    // — see AutoGAConfig) since it has a different crossover enum and no
    // selectable mutation strategy at all.
    //
    // popSize/generations/mutationRate/crossoverRate are read ONCE from
    // conf and shared by every row in whichever sweep runs — they used to
    // be hardcoded per row, which looked configurable but was dead
    // weight: setOperators()/setAutoOperators() only apply mRate/cRate
    // when `!confLoaded`, never true once a GaConf has loaded.
    //
    // No timing (removed entirely). Leaderboard reports ACTUAL generations
    // run (early-stop aware, via lastGenerationsRun) and real coverage%
    // (via lastCoveragePercent), sorted by coverage% descending.
    int run_benchmark(BenchmarkCase which = BenchmarkCase::HOMOGENOUS) {
        // Fixed-length sweep: HOMOGENOUS and HETEROGENOUS.
        std::vector<GAConfig> sweepConfigs = {
            // Control & Baseline
            {"Baseline Control",     SelectionMethod::TOURNAMENT,     CrossoverMethod::SINGLE_POINT, MutationMethod::SPATIAL_NUDGE},

            // Lightweight / Fast Profiles
            {"Lightweight (Fast)",   SelectionMethod::TOURNAMENT,     CrossoverMethod::SINGLE_POINT, MutationMethod::SPATIAL_NUDGE},
            {"Micro-Sweep (SUS)",    SelectionMethod::SUS,            CrossoverMethod::MULTI_POINT,  MutationMethod::GAUSSIAN_CREEP},

            // Intensive Searches
            {"Intensive Search",     SelectionMethod::TOURNAMENT,     CrossoverMethod::MULTI_POINT,  MutationMethod::SPATIAL_NUDGE},
            {"Intensive SUS",        SelectionMethod::SUS,            CrossoverMethod::UNIFORM,      MutationMethod::GAUSSIAN_CREEP},

            // Exploration-Heavy
            {"Uniform Exploration",  SelectionMethod::RANK_BASED,     CrossoverMethod::UNIFORM,      MutationMethod::GAUSSIAN_CREEP},
            {"SUS High-Mutation",    SelectionMethod::SUS,            CrossoverMethod::UNIFORM,      MutationMethod::GAUSSIAN_CREEP},

            // Exploitation-Heavy
            {"Exploitation Focus",   SelectionMethod::TOURNAMENT,     CrossoverMethod::MULTI_POINT,  MutationMethod::SPATIAL_NUDGE},

            // Roulette Wheel Variations
            {"Roulette & Uniform",   SelectionMethod::ROULETTE_WHEEL, CrossoverMethod::UNIFORM,      MutationMethod::GAUSSIAN_CREEP},
            {"Roulette & Nudge",     SelectionMethod::ROULETTE_WHEEL, CrossoverMethod::MULTI_POINT,  MutationMethod::SPATIAL_NUDGE},

            // Rank-Based Variations
            {"Rank & Gaussian",      SelectionMethod::RANK_BASED,     CrossoverMethod::MULTI_POINT,  MutationMethod::GAUSSIAN_CREEP},
            {"Rank & Nudge",         SelectionMethod::RANK_BASED,     CrossoverMethod::SINGLE_POINT, MutationMethod::SPATIAL_NUDGE},

            // Single Point Reset Mutation Variations
            {"Tournament & Reset",   SelectionMethod::TOURNAMENT,     CrossoverMethod::SINGLE_POINT, MutationMethod::SINGLE_POINT},
            {"SUS & Reset",          SelectionMethod::SUS,            CrossoverMethod::SINGLE_POINT, MutationMethod::SINGLE_POINT},
            {"Rank & Reset",         SelectionMethod::RANK_BASED,     CrossoverMethod::SINGLE_POINT, MutationMethod::SINGLE_POINT},

            // High Crossover Exploration
            {"High-Cross SUS Nudge", SelectionMethod::SUS,            CrossoverMethod::SINGLE_POINT, MutationMethod::SPATIAL_NUDGE}
        };

        // Variable-length sweep: AUTO. Full factorial over AUTO's actual
        // parameter space (4 SelectionMethod x 3 AutoCrossoverMethod = 12
        // rows) rather than a curated subset, since there's no natural
        // basis to prefer some combinations over others here the way the
        // fixed-length table's hand-picked profiles do.
        std::vector<AutoGAConfig> sweepConfigsAuto = {
            {"Tournament + RandomSubset",   SelectionMethod::TOURNAMENT,     AutoCrossoverMethod::RANDOM_SUBSET},
            {"Tournament + LengthNorm",     SelectionMethod::TOURNAMENT,     AutoCrossoverMethod::LENGTH_NORMALIZED},
            {"Tournament + PositionalPad",  SelectionMethod::TOURNAMENT,     AutoCrossoverMethod::POSITIONAL_PADDING},

            {"Roulette + RandomSubset",     SelectionMethod::ROULETTE_WHEEL, AutoCrossoverMethod::RANDOM_SUBSET},
            {"Roulette + LengthNorm",       SelectionMethod::ROULETTE_WHEEL, AutoCrossoverMethod::LENGTH_NORMALIZED},
            {"Roulette + PositionalPad",    SelectionMethod::ROULETTE_WHEEL, AutoCrossoverMethod::POSITIONAL_PADDING},

            {"RankBased + RandomSubset",    SelectionMethod::RANK_BASED,     AutoCrossoverMethod::RANDOM_SUBSET},
            {"RankBased + LengthNorm",      SelectionMethod::RANK_BASED,     AutoCrossoverMethod::LENGTH_NORMALIZED},
            {"RankBased + PositionalPad",   SelectionMethod::RANK_BASED,     AutoCrossoverMethod::POSITIONAL_PADDING},

            {"SUS + RandomSubset",          SelectionMethod::SUS,            AutoCrossoverMethod::RANDOM_SUBSET},
            {"SUS + LengthNorm",            SelectionMethod::SUS,            AutoCrossoverMethod::LENGTH_NORMALIZED},
            {"SUS + PositionalPad",         SelectionMethod::SUS,            AutoCrossoverMethod::POSITIONAL_PADDING}
        };

        const bool isAuto  = (which == BenchmarkCase::AUTO);
        const size_t nRows = isAuto ? sweepConfigsAuto.size() : sweepConfigs.size();

        ThreadPool sharedPool;
        GaEngine gaEngine(sharedPool);

        if (!gaEngine.confLoaded) {
            std::cerr << "[ga] run_benchmark: no GaConf loaded — every profile's "
                         "popSize/generations/mutationRate/crossoverRate comes "
                         "from config, so there is nothing to sweep with. Aborting.\n";
            return -1;
        }

        const char* caseName = which == BenchmarkCase::HOMOGENOUS   ? "HOMOGENOUS"
                              : which == BenchmarkCase::HETEROGENOUS ? "HETEROGENOUS"
                                                                      : "AUTO";

        std::cout << "========================================================================\n";
        std::cout << "   STARTING HYPERPARAMETER OPTIMIZATION SWEEP (" << nRows
                  << " profiles, " << caseName << ")\n";
        std::cout << "========================================================================\n";
        std::cout << "  popSize="        << gaEngine.popSize
                  << "  generations="    << gaEngine.conf.parameters.numGenerations
                  << "  mutationRate="   << gaEngine.mutationRate
                  << "  crossoverRate="  << gaEngine.crossoverRate
                  << "  (uniform across every profile, from conf)\n";
        std::cout << "Running evaluation sweep over C++ ThreadPool...\n\n";

        std::vector<BenchmarkResult> leaderboard;
        leaderboard.reserve(nRows);

        for (size_t idx = 0; idx < nRows; ++idx) {
            std::string rowName, selectionLabel, crossoverLabel, mutationLabel;

            if (isAuto) {
                const auto& c = sweepConfigsAuto[idx];
                rowName        = c.name;
                selectionLabel = toString(c.selection);
                crossoverLabel = toString(c.crossover);
                mutationLabel  = "Auto";   // no selectable mutation strategy in run_auto

                std::cout << " [*] Running Test " << std::setw(2) << (idx + 1) << " / " << nRows
                          << ": " << std::left << std::setw(28) << rowName << " ... " << std::flush;

                gaEngine.setAutoOperators(c.selection, c.crossover, 0.0, 0.0);
                gaEngine.run_auto();

            } else {
                const auto& c = sweepConfigs[idx];
                rowName        = c.name;
                selectionLabel = toString(c.selection);
                crossoverLabel = toString(c.crossover);
                mutationLabel  = toString(c.mutation);

                std::cout << " [*] Running Test " << std::setw(2) << (idx + 1) << " / " << nRows
                          << ": " << std::left << std::setw(28) << rowName << " ... " << std::flush;

                // HOMOGENOUS and HETEROGENOUS share the same fixed-length
                // operator set, so no mapping is needed between them.
                gaEngine.setOperators(c.selection, c.crossover, c.mutation, 0.0, 0.0);
                if (which == BenchmarkCase::HOMOGENOUS) gaEngine.run_homogenous();
                else                                    gaEngine.run_heterogenous();
            }

            leaderboard.push_back({
                rowName,
                selectionLabel,
                crossoverLabel,
                mutationLabel,
                gaEngine.lastGenerationsRun,
                gaEngine.lastCoveragePercent
            });

            std::cout << "Done. (Gen: " << gaEngine.lastGenerationsRun
                      << ", Coverage: " << std::fixed << std::setprecision(2)
                      << gaEngine.lastCoveragePercent << "%)\n";
        }

        // Best (highest) coverage first.
        std::sort(leaderboard.begin(), leaderboard.end(), [](const BenchmarkResult& a, const BenchmarkResult& b) {
            return a.coveragePercent > b.coveragePercent;
            });

        std::cout << "\n\n";
        std::cout << "=================================================================================\n";
        std::cout << "                      HYPERPARAMETER OPTIMIZATION LEADERBOARD (" << caseName << ")\n";
        std::cout << "=================================================================================\n";
        std::cout << " Rank | " << std::left << std::setw(28) << "Config Profile Name"
                  << " | " << std::setw(11) << "Selection"
                  << " | " << std::setw(13) << "Crossover"
                  << " | " << std::setw(9)  << "Mutation"
                  << " | " << std::right << std::setw(4) << "Gen"
                  << " | " << std::setw(12) << "Coverage %"
                  << " |\n";
        std::cout << "------|------------------------------|-------------|---------------|-----------|------|--------------|\n";

        for (size_t i = 0; i < leaderboard.size(); ++i) {
            const auto& r = leaderboard[i];
            std::cout << "  " << std::right << std::setw(2) << (i + 1) << "  | "
                      << std::left << std::setw(28) << r.name << " | "
                      << std::setw(11) << r.selection << " | "
                      << std::setw(13) << r.crossover << " | "
                      << std::setw(9)  << r.mutation << " | "
                      << std::right << std::setw(4) << r.generations << " | "
                      << std::setw(11) << std::fixed << std::setprecision(2) << r.coveragePercent << "%"
                      << " |\n";
        }
        std::cout << "=================================================================================\n";
        std::cout << "Note: sorted by coverage % (descending, highest first).\n\n";

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

    std::string toString(AutoCrossoverMethod cm) {
        switch(cm) {
            case AutoCrossoverMethod::RANDOM_SUBSET:      return "RandomSubset";
            case AutoCrossoverMethod::LENGTH_NORMALIZED:  return "LengthNorm";
            case AutoCrossoverMethod::POSITIONAL_PADDING: return "PositionalPad";
        }
        return "Unknown";
    }

    double getFitnessEvalTime() const { return totalFitnessEvalTimeMs; }
    double getEvolutionaryTime() const { return totalEvolutionaryTimeMs; }

};
