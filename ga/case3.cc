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
#include <numeric>       // std::iota

// Only include common.h for Global, Bin, Config, and FeasibleGrid structures.
#include "../common.h"

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
enum class SelectionMethod { TOURNAMENT, ROULETTE_WHEEL, RANK_BASED, SUS };

// Three crossover strategies designed for variable-length chromosomes.
// One is chosen at random each time crossover fires.
enum class CrossoverMethod { RANDOM_SUBSET, LENGTH_NORMALIZED, POSITIONAL_PADDING };

// ============================================================
// RADAR TYPE SPECIFICATION
// Each entry describes one class of radar:
//   rangeKm  — effective detection range in km
//   minCount — fewest radars of this class the solution must deploy
//   maxCount — most radars of this class the solution may deploy
//
// Example: { 50.0, 2, 5 }  means "2 to 5 short-range radars"
// ============================================================
struct RadarConfig {
    double rangeKm;
    int    minCount;
    int    maxCount;
};

// ============================================================
// GENE — a single radar deployment
//   r, c       : grid cell (row, column)
//   typeIndex  : which RadarConfig entry this radar belongs to
//                -1 is a NULL sentinel used only during
//                positional-padding crossover; always stripped before evaluation.
// ============================================================
struct RadarGene {
    int r, c;
    int typeIndex;
};

// ============================================================
// CHROMOSOME — variable-length array of radar genes
// Length can vary between sum(minCount) and sum(maxCount)
// across all types.
// ============================================================
struct Chromosome {
    std::vector<RadarGene> radars;
    double fitness = 0.0;

    // Count genes belonging to a specific radar type
    int countType(int t) const {
        return (int)std::count_if(radars.begin(), radars.end(),
            [t](const RadarGene& g) { return g.typeIndex == t; });
    }

    int totalCount() const { return (int)radars.size(); }
};

// ============================================================
// LOCAL COVERAGE CONTAINER (unchanged)
// ============================================================
struct LocalCoverage {
    std::vector<uint8_t> mask;
    int rowMin, rowMax, colMin, colMax;
    int localCols, localRows;
};

// ============================================================
// THREAD POOL (unchanged)
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
            if (stop) throw std::runtime_error("Enqueue on stopped ThreadPool");
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
        for (std::thread& w : workers)
            if (w.joinable()) w.join();
    }
};

// ============================================================
// GENETIC ALGORITHM
// ============================================================
class ga {
    Global g;
    std::vector<Chromosome> population;
    ThreadPool& pool;

    // --- Hyperparameters ---
    int    popSize       = NUM_CHROMOSOMES;
    double mutationRate  = 0.10;   // P(a chromosome gets any mutation)
    double crossoverRate = 0.85;

    // Radar cost coefficient in F(k) = coverage_score - λ * n
    //   λ = 0.0  → radars are free; pure coverage maximisation
    //   λ > 0    → penalise each additional radar; algorithm learns
    //              to cover with fewer units
    double lambda = 0.0;

    SelectionMethod currentSelection = SelectionMethod::TOURNAMENT;
    CrossoverMethod currentCrossover = CrossoverMethod::RANDOM_SUBSET;

    std::mt19937 rng;

    // ============================================================
    // HELPERS
    // ============================================================

    // Encode (r,c) as a single 64-bit key for O(1) set lookups
    static int64_t posKey(int r, int c) {
        return static_cast<int64_t>(r) * 1000000LL + c;
    }

    // Extract plain GridPoints from a gene vector (for PNG helpers that
    // still expect a vector<GridPoint>)
    static std::vector<GridPoint> extractPositions(const std::vector<RadarGene>& genes) {
        std::vector<GridPoint> pts;
        pts.reserve(genes.size());
        for (const auto& gene : genes)
            if (gene.typeIndex >= 0) pts.push_back({gene.r, gene.c});
        return pts;
    }

    // ============================================================
    // RADAR COMPUTATION (identical to original)
    // ============================================================

    static std::vector<GridPoint> get_circle_boundary(int cr, int cc, int radius) {
        std::vector<GridPoint> pts;
        pts.reserve(8 * radius);
        int x = 0, y = radius, d = 3 - 2 * radius;
        auto push8 = [&](int px, int py) {
            pts.push_back({cr + px, cc + py}); pts.push_back({cr - px, cc + py});
            pts.push_back({cr + px, cc - py}); pts.push_back({cr - px, cc - py});
            pts.push_back({cr + py, cc + px}); pts.push_back({cr - py, cc + px});
            pts.push_back({cr + py, cc - px}); pts.push_back({cr - py, cc - px});
        };
        while (y >= x) {
            push8(x, y); x++;
            if (d > 0) { y--; d += 4 * (x - y) + 10; } else { d += 4 * x + 6; }
        }
        std::sort(pts.begin(), pts.end(), [](const GridPoint& a, const GridPoint& b) {
            return a.r < b.r || (a.r == b.r && a.c < b.c);
        });
        pts.erase(std::unique(pts.begin(), pts.end(),
            [](const GridPoint& a, const GridPoint& b) { return a.r == b.r && a.c == b.c; }),
            pts.end());
        return pts;
    }

    static void cast_radar_ray_local(
        int srcRow, int srcCol, int endRow, int endCol,
        double srcElev, double targetHeight,
        double minSlope, double maxSlope, double horizonM, int maxRadiusCells,
        const short* __restrict__ dem, uint8_t* __restrict__ stat,
        int gridCols, int gridRows,
        int rowMin, int colMin, int colMax)
    {
        int adr = std::abs(endRow - srcRow), adc = std::abs(endCol - srcCol);
        int sr  = (srcRow < endRow) ? 1 : -1;
        int sc  = (srcCol < endCol) ? 1 : -1;
        int err = adr - adc;

        const std::ptrdiff_t demStrideR  = static_cast<std::ptrdiff_t>(sr) * gridCols;
        const std::ptrdiff_t demStrideC  = sc;
        std::ptrdiff_t demOff = static_cast<std::ptrdiff_t>(srcRow) * gridCols + srcCol;
        const short* demCell = dem + demOff;

        int localCols = colMax - colMin;
        const std::ptrdiff_t statStrideR = static_cast<std::ptrdiff_t>(sr) * localCols;
        const std::ptrdiff_t statStrideC = sc;
        std::ptrdiff_t statOff = static_cast<std::ptrdiff_t>(srcRow - rowMin) * localCols
                               + (srcCol - colMin);
        uint8_t* statCell = stat + statOff;

        int    row = srcRow, col = srcCol;
        double dist = 0.0, curvDrop = 0.0;
        double maxTerrainSlope = -1e18;
        int    cellCount = 0;

        while (true) {
            if (row == endRow && col == endCol) break;

            int e2 = 2 * err;
            bool mr = (e2 > -adc), mc = (e2 < adr);
            if (mr) { err -= adc; row += sr; demCell  += demStrideR; statCell += statStrideR; }
            if (mc) { err += adr; col += sc; demCell  += demStrideC; statCell += statStrideC; }

            double step = (mr && mc) ? RESOLUTION * 1.41421356237 : RESOLUTION;
            curvDrop += (2.0 * dist * step + step * step) * INV_2_KEFFR;
            dist     += step;
            ++cellCount;

            if (cellCount > maxRadiusCells) break;
            if (dist > horizonM) break;

            short raw = *demCell;
            if (raw < -9000) break;

            double terrainElev  = static_cast<double>(raw) - curvDrop;
            double terrainSlope = (terrainElev - srcElev) / dist;
            double targetSlope  = terrainSlope + targetHeight / dist;

            if (terrainSlope > maxTerrainSlope) {
                maxTerrainSlope = terrainSlope;
                *statCell = 1;
            } else {
                *statCell = 0;
            }
        }
    }

    static LocalCoverage compute_local_coverage(
        Global& g, int srcRow, int srcCol,
        double antennaHeight, double targetHeight,
        double minElevDeg, double maxElevDeg, double maxRangeM)
    {
        int gridCols = g.conf.grid.cols;
        int gridRows = g.conf.grid.rows;

        double horizonM = std::sqrt(2.0 * K_EFF_R * antennaHeight)
                        + std::sqrt(2.0 * K_EFF_R * targetHeight);
        if (horizonM > maxRangeM) horizonM = maxRangeM;

        double minSlopeRatio = std::tan(minElevDeg * DEG2RAD);
        double maxSlopeRatio = std::tan(maxElevDeg * DEG2RAD);

        int rawRadius   = (int)std::ceil(maxRangeM / RESOLUTION);
        int maxPossible = std::min({srcRow, gridRows - srcRow - 1,
                                    srcCol, gridCols - srcCol - 1});
        int maxRadiusCells = std::min(rawRadius, std::max(0, maxPossible));

        int rowMin = std::max(0,        srcRow - maxRadiusCells - 1);
        int rowMax = std::min(gridRows, srcRow + maxRadiusCells + 2);
        int colMin = std::max(0,        srcCol - maxRadiusCells - 1);
        int colMax = std::min(gridCols, srcCol + maxRadiusCells + 2);

        LocalCoverage cov;
        cov.rowMin   = rowMin; cov.rowMax   = rowMax;
        cov.colMin   = colMin; cov.colMax   = colMax;
        cov.localCols = colMax - colMin;
        cov.localRows = rowMax - rowMin;
        cov.mask.assign((size_t)cov.localRows * cov.localCols, 0);

        if (srcRow >= 0 && srcRow < gridRows && srcCol >= 0 && srcCol < gridCols) {
            short  srcTerrain = g.bin.elev_matrix[(size_t)srcRow * gridCols + srcCol];
            double srcElev    = static_cast<double>(srcTerrain) + antennaHeight;

            auto boundary = get_circle_boundary(srcRow, srcCol, maxRadiusCells);
            cov.mask[(size_t)(srcRow - rowMin) * cov.localCols + (srcCol - colMin)] = 1;

            const short* demPtr  = g.bin.elev_matrix.data();
            uint8_t*     statPtr = cov.mask.data();

            for (const auto& bp : boundary) {
                int er = bp.r, ec = bp.c;
                if (er >= 0 && ec >= 0 && er < gridRows && ec < gridCols) {
                    cast_radar_ray_local(
                        srcRow, srcCol, er, ec,
                        srcElev, targetHeight, minSlopeRatio, maxSlopeRatio,
                        horizonM, maxRadiusCells, demPtr, statPtr,
                        gridCols, gridRows, rowMin, colMin, colMax
                    );
                }
            }
        }
        return cov;
    }

    // ============================================================
    // COVERAGE EVALUATION
    // Now takes vector<RadarGene> — each gene carries its own typeIndex
    // so the correct range is looked up from configs[typeIndex].
    //
    // result[0] = proximity violation count
    // result[1] = cells with single coverage
    // result[2] = cells with dual coverage
    // result[3] = cells with triple+ coverage
    // ============================================================
    static std::vector<int> calculate_radar_coverage_overlap(
        Global& g, FeasibleGrid grid,
        const std::vector<RadarGene>& placements,
        const std::vector<RadarConfig>& configs,
        double minSpacingCells = 5.0,
        const std::string& pngPath = "")
    {
        std::vector<int> result(4, 0);
        if (placements.empty() || grid.cells.empty()) return result;

        // --- Proximity violations ---
        for (size_t i = 0; i < placements.size(); i++) {
            for (size_t j = i + 1; j < placements.size(); j++) {
                double dr = placements[i].r - placements[j].r;
                double dc = placements[i].c - placements[j].c;
                if (std::sqrt(dr*dr + dc*dc) < minSpacingCells) result[0]++;
            }
        }

        // --- Per-radar LOS coverage using type-specific range ---
        for (const auto& p : placements) {
            if (p.typeIndex < 0 || p.typeIndex >= (int)configs.size()) continue;
            double maxRangeM = configs[p.typeIndex].rangeKm * 1000.0;

            LocalCoverage cov = compute_local_coverage(
                g, p.r, p.c, 50.0, 1000.0, 10.0, 30.0, maxRangeM);

            for (auto& cell : grid.cells) {
                if (cell.r >= cov.rowMin && cell.r < cov.rowMax &&
                    cell.c >= cov.colMin && cell.c < cov.colMax) {
                    size_t flat = (size_t)(cell.r - cov.rowMin) * cov.localCols
                                + (cell.c - cov.colMin);
                    if (flat < cov.mask.size() && cov.mask[flat] == 1)
                        cell.overlap++;
                }
            }
        }

        for (const auto& cell : grid.cells) {
            if      (cell.overlap == 1) result[1]++;
            else if (cell.overlap == 2) result[2]++;
            else if (cell.overlap >= 3) result[3]++;
        }

        if (!pngPath.empty()) {
            auto pos = extractPositions(placements);
            Bin::writeFeasibleGridOverlapPng(grid, pngPath, pos);
        }
        return result;
    }

    // ============================================================
    // FITNESS
    //
    //   F(k) = single_cov×1 + dual_cov×2 + triple_cov×3
    //         - proximity_violations × 5000
    //         - λ × n
    //
    //   n = total number of radars in this chromosome
    //   λ = radar cost coefficient (set via setLambda)
    // ============================================================
    static void evaluate_fitness(
        Global& g, const FeasibleGrid& base_grid, Chromosome& cn, int index,
        const std::vector<RadarConfig>& configs, double lambda, bool writePng)
    {
        std::string png_path = writePng
            ? "chromosome_" + std::to_string(index) + "_coverage.png" : "";

        auto res = calculate_radar_coverage_overlap(
            g, base_grid, cn.radars, configs, 5.0, png_path);

        int n = cn.totalCount();
        cn.fitness = (res[1] * 1.0) + (res[2] * 2.0) + (res[3] * 3.0)
                   - (res[0] * 5000.0)
                   - lambda * n;
    }

    void parallelEvaluatePopulation(
        std::vector<Chromosome>& pop, Global& globalRef,
        const FeasibleGrid& base_grid,
        const std::vector<RadarConfig>& configs,
        double lambda, bool writePng = false)
    {
        std::vector<std::future<void>> futures;
        futures.reserve(pop.size());
        for (size_t i = 0; i < pop.size(); i++) {
            futures.push_back(pool.enqueue(
                evaluate_fitness,
                std::ref(globalRef), std::ref(base_grid), std::ref(pop[i]),
                (int)i, std::ref(configs), lambda, writePng));
        }
        for (auto& f : futures) f.wait();
    }

    // ============================================================
    // SNAP TO NEAREST FEASIBLE CELL
    // Used by MOVE mutation and REPAIR to clamp out-of-bounds positions.
    // ============================================================
    GridPoint snapToFeasible(int r, int c, const FeasibleGrid& f) {
        int minDistSq = INT_MAX;
        GridPoint best = {f.cells.front().r, f.cells.front().c};
        for (const auto& cell : f.cells) {
            int dr = cell.r - r, dc = cell.c - c;
            int d  = dr*dr + dc*dc;
            if (d < minDistSq) {
                minDistSq = d;
                best = {cell.r, cell.c};;
                if (d == 0) break; // exact hit
            }
        }
        return best;
    }

    // ============================================================
    // MUTATION  — three types, one chosen uniformly at random
    // ============================================================

    // ADD: insert one new radar whose type is still under its maxCount.
    // Picks a random eligible type, then a random unoccupied feasible cell.
    void mutate_add(Chromosome& ind, const FeasibleGrid& f,
                    const std::vector<RadarConfig>& configs)
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
        // Silent fail: all cells occupied (extremely unlikely in practice)
    }

    // REMOVE: drop one radar whose type count is still above its minCount.
    void mutate_remove(Chromosome& ind, const std::vector<RadarConfig>& configs) {
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

    // MOVE: nudge an existing radar to a nearby feasible cell.
    // The count of every radar type stays unchanged.
    void mutate_move(Chromosome& ind, const FeasibleGrid& f) {
        if (ind.radars.empty()) return;

        std::uniform_int_distribution<int> distGene(0, (int)ind.radars.size()-1);
        std::uniform_int_distribution<int> distNudge(-10, 10);

        auto& gene = ind.radars[distGene(rng)];
        if (gene.typeIndex < 0) return; // skip NULL sentinel

        GridPoint snapped = snapToFeasible(
            gene.r + distNudge(rng),
            gene.c + distNudge(rng), f);
        gene.r = snapped.r;
        gene.c = snapped.c;
    }

    // Execute one randomly chosen mutation type (if mutation roll passes)
    void executeMutation(Chromosome& ind, const FeasibleGrid& f,
                         const std::vector<RadarConfig>& configs)
    {
        std::uniform_real_distribution<double> distR(0.0, 1.0);
        if (distR(rng) > mutationRate) return;

        switch (std::uniform_int_distribution<int>(0, 2)(rng)) {
            case 0: mutate_add(ind, f, configs);    break;
            case 1: mutate_remove(ind, configs);    break;
            case 2: mutate_move(ind, f);            break;
        }
    }

    // ============================================================
    // REPAIR — called after crossover AND after mutation
    //
    // Steps:
    //  1. Strip NULL sentinel genes (from positional-padding crossover)
    //  2. Deduplicate coincident positions (keep first occurrence)
    //  3. Trim types that exceed their maxCount (random removal)
    //  4. Fill types that are below their minCount (random insertion)
    // ============================================================
    void repairChromosome(Chromosome& ind, const FeasibleGrid& f,
                          const std::vector<RadarConfig>& configs)
    {
        // 1. Remove NULL sentinels
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

        // 3. Trim excess: remove random genes of any type that exceeds maxCount
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

        // 4. Fill deficit: add random genes for types below minCount
        {
            std::unordered_set<int64_t> occupied;
            for (const auto& gene : ind.radars)
                occupied.insert(posKey(gene.r, gene.c));

            std::uniform_int_distribution<int> distF(0, (int)f.cells.size()-1);

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
                    if (!placed) break; // grid too crowded — extremely rare
                    deficit--;
                }
            }
        }
    }

    // ============================================================
    // CROSSOVER — three variable-length methods, chosen at random
    // ============================================================

    // Method 1: RANDOM_SUBSET
    // Partition each parent randomly into a "chosen" and "rest" half.
    // Child 1 = chosen(A) ∪ chosen(B)
    // Child 2 = rest(A)   ∪ rest(B)
    // Pro: inherits actual positions from both parents.
    // Con: child length is non-deterministic; repairChromosome fixes it.
    std::pair<Chromosome, Chromosome> crossover_randomSubset(
        const Chromosome& p1, const Chromosome& p2)
    {
        auto splitParent = [&](const Chromosome& parent)
            -> std::pair<std::vector<RadarGene>, std::vector<RadarGene>>
        {
            int n = (int)parent.radars.size();
            int k = (n > 0) ? std::uniform_int_distribution<int>(0, n)(rng) : 0;
            std::vector<int> idx(n);
            std::iota(idx.begin(), idx.end(), 0);
            std::shuffle(idx.begin(), idx.end(), rng);
            std::vector<RadarGene> chosen, rest;
            for (int i = 0; i < n; i++) {
                (i < k ? chosen : rest).push_back(parent.radars[idx[i]]);
            }
            return {chosen, rest};
        };

        auto [a_chosen, a_rest] = splitParent(p1);
        auto [b_chosen, b_rest] = splitParent(p2);

        Chromosome c1, c2;
        c1.radars.insert(c1.radars.end(), a_chosen.begin(), a_chosen.end());
        c1.radars.insert(c1.radars.end(), b_chosen.begin(), b_chosen.end());
        c2.radars.insert(c2.radars.end(), a_rest.begin(),   a_rest.end());
        c2.radars.insert(c2.radars.end(), b_rest.begin(),   b_rest.end());
        return {c1, c2};
    }

    // Method 2: LENGTH_NORMALIZED
    // Decide child length = ⌈(|A| + |B|) / 2⌉, then sample from the
    // combined gene pool of both parents.
    // Pro: keeps length stable across generations.
    // Con: loses organic length drift; use with REMOVE/ADD mutations.
    std::pair<Chromosome, Chromosome> crossover_lengthNormalized(
        const Chromosome& p1, const Chromosome& p2)
    {
        int avgLen = ((int)p1.radars.size() + (int)p2.radars.size() + 1) / 2;

        std::vector<RadarGene> pool;
        pool.insert(pool.end(), p1.radars.begin(), p1.radars.end());
        pool.insert(pool.end(), p2.radars.begin(), p2.radars.end());
        int take = std::min(avgLen, (int)pool.size());

        auto buildChild = [&]() {
            std::shuffle(pool.begin(), pool.end(), rng);
            Chromosome child;
            child.radars.assign(pool.begin(), pool.begin() + take);
            return child;
        };

        return {buildChild(), buildChild()};
    }

    // Method 3: POSITIONAL_PADDING
    // Pad the shorter parent with NULL (-1) sentinels to match the longer,
    // apply a single-point cut, then strip all NULLs from children.
    // Pro: mirrors standard crossover; preserves positional structure.
    // Con: NULL placement is arbitrary and can bias which genes are inherited.
    std::pair<Chromosome, Chromosome> crossover_positionalPadding(
        const Chromosome& p1, const Chromosome& p2)
    {
        int len1 = (int)p1.radars.size(), len2 = (int)p2.radars.size();
        int maxLen = std::max(len1, len2);
        if (maxLen == 0) return {p1, p2};

        const RadarGene null_gene = {0, 0, -1};
        std::vector<RadarGene> pad1(p1.radars), pad2(p2.radars);
        while ((int)pad1.size() < maxLen) pad1.push_back(null_gene);
        while ((int)pad2.size() < maxLen) pad2.push_back(null_gene);

        int cut = std::uniform_int_distribution<int>(1, maxLen)(rng);

        Chromosome c1, c2;
        for (int i = 0; i < maxLen; i++) {
            if (i < cut) {
                c1.radars.push_back(pad1[i]);
                c2.radars.push_back(pad2[i]);
            } else {
                c1.radars.push_back(pad2[i]);
                c2.radars.push_back(pad1[i]);
            }
        }

        // Strip NULL sentinels (repairChromosome does this too, but be explicit)
        auto stripNull = [](std::vector<RadarGene>& v) {
            v.erase(std::remove_if(v.begin(), v.end(),
                [](const RadarGene& g) { return g.typeIndex < 0; }), v.end());
        };
        stripNull(c1.radars);
        stripNull(c2.radars);

        return {c1, c2};
    }

    // Execute a randomly chosen crossover method (unless roll misses)
    std::pair<Chromosome, Chromosome> executeCrossover(
        const Chromosome& p1, const Chromosome& p2)
    {
        if (std::uniform_real_distribution<double>(0.0, 1.0)(rng) > crossoverRate)
            return {p1, p2};

        switch (std::uniform_int_distribution<int>(0, 2)(rng)) {
            case 0: return crossover_randomSubset(p1, p2);
            case 1: return crossover_lengthNormalized(p1, p2);
            default: return crossover_positionalPadding(p1, p2);
        }
    }

    // ============================================================
    // SELECTION (unchanged logic, adapted for new Chromosome type)
    // ============================================================

    Chromosome tournamentSelection(const std::vector<Chromosome>& pop) {
        std::uniform_int_distribution<int> dist(0, (int)pop.size()-1);
        Chromosome best = pop[dist(rng)];
        for (int i = 0; i < 2; ++i) {
            const Chromosome& comp = pop[dist(rng)];
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

        double slice = std::uniform_real_distribution<double>(0.0, total)(rng);
        double accum = 0.0;
        for (const auto& ind : pop) {
            accum += (ind.fitness + offset);
            if (accum >= slice) return ind;
        }
        return pop.back();
    }

    Chromosome rankBasedSelection(const std::vector<Chromosome>& sorted_pop) {
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

    std::vector<Chromosome> stochasticUniversalSampling(
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

    Chromosome executeSelection(const std::vector<Chromosome>& pop) {
        if (currentSelection == SelectionMethod::TOURNAMENT)     return tournamentSelection(pop);
        if (currentSelection == SelectionMethod::ROULETTE_WHEEL) return rouletteWheelSelection(pop);
        return rankBasedSelection(pop); // RANK_BASED (SUS handled in run_evolution)
    }

    // ============================================================
    // INITIAL CHROMOSOME GENERATION
    // Picks a random count in [minCount, maxCount] for each type,
    // then places genes at unique random feasible cells.
    // ============================================================
    Chromosome generateRandomChromosome(
        const FeasibleGrid& f, const std::vector<RadarConfig>& configs)
    {
        Chromosome cn;
        std::unordered_set<int64_t> occupied;
        std::uniform_int_distribution<int> distF(0, (int)f.cells.size()-1);

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

public:
    explicit ga(ThreadPool& p) : pool(p) {
        rng.seed((unsigned int)
            std::chrono::system_clock::now().time_since_epoch().count());
        std::cout << "[ga] Pre-loading all .bin matrices...\n";
        g.bin.init();
    }

    // λ = 0    → radars are free; algorithm maximises coverage only
    // λ = 10   → mild preference for fewer radars
    // λ = 1000 → strong pressure to minimise radar count
    void setLambda(double l) { lambda = l; }

    void setOperators(SelectionMethod sm, CrossoverMethod cm,
                      double mRate, double cRate) {
        currentSelection = sm;
        currentCrossover = cm;
        mutationRate     = mRate;
        crossoverRate    = cRate;
    }

    // ============================================================
    // MAIN EVOLUTION LOOP
    //
    // configs  — heterogeneous fleet; each entry has rangeKm + [min,max] count
    // generations — number of breeding cycles after the seed population
    // ============================================================
    void run_evolution(const std::vector<RadarConfig>& configs, int generations) {
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
        std::cout << "  Total radar range : [" << totalMin << ", " << totalMax << "]\n";
        std::cout << "  Lambda (cost/radar): " << lambda << "\n\n";

        std::cout << "[ga] Extracting feasible grid...\n";
        auto f = Bin::getFeasibleGrid(g.bin, WRITE);

        // --- Seed population (Generation 0) ---
        population.clear();
        population.reserve(popSize);
        for (int i = 0; i < popSize; i++)
            population.push_back(generateRandomChromosome(f, configs));

        std::cout << "[ga] Evaluating Generation 0...\n";
        parallelEvaluatePopulation(population, g, f, configs, lambda, false);

        // --- Evolution ---
        for (int gen = 1; gen <= generations; ++gen) {
            // Sort ascending; best individual sits at back (easy elitism)
            std::sort(population.begin(), population.end(),
                [](const Chromosome& a, const Chromosome& b) {
                    return a.fitness < b.fitness;
                });

            std::vector<Chromosome> nextGen;
            nextGen.reserve(popSize);

            // Elitism: copy single best unchanged
            nextGen.push_back(population.back());

            // SUS fills a mating pool in one pass; others select on demand
            bool useSUS = (currentSelection == SelectionMethod::SUS);
            std::vector<Chromosome> matingPool;
            size_t poolIdx = 0;
            if (useSUS)
                matingPool = stochasticUniversalSampling(population, popSize - 1);

            while ((int)nextGen.size() < popSize) {
                Chromosome p1, p2;
                if (useSUS) {
                    p1 = matingPool[poolIdx % matingPool.size()]; poolIdx++;
                    p2 = matingPool[poolIdx % matingPool.size()]; poolIdx++;
                } else {
                    p1 = executeSelection(population);
                    p2 = executeSelection(population);
                }

                // Crossover → repair → mutate → repair
                auto [c1, c2] = executeCrossover(p1, p2);

                repairChromosome(c1, f, configs);
                repairChromosome(c2, f, configs);

                executeMutation(c1, f, configs);
                executeMutation(c2, f, configs);

                repairChromosome(c1, f, configs);
                repairChromosome(c2, f, configs);

                nextGen.push_back(c1);
                if ((int)nextGen.size() < popSize) nextGen.push_back(c2);
            }

            parallelEvaluatePopulation(nextGen, g, f, configs, lambda, false);
            population = std::move(nextGen);

            auto best = std::max_element(population.begin(), population.end(),
                [](const Chromosome& a, const Chromosome& b) {
                    return a.fitness < b.fitness;
                });

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
        auto finalBest = std::max_element(population.begin(), population.end(),
            [](const Chromosome& a, const Chromosome& b) {
                return a.fitness < b.fitness;
            });

        std::cout << "\n[ga] Evolution complete.\n";
        std::cout << "  Final fitness  : " << finalBest->fitness << "\n";
        std::cout << "  Total radars   : " << finalBest->totalCount() << "\n";
        for (int t = 0; t < (int)configs.size(); t++) {
            std::cout << "  Type " << t << " (" << configs[t].rangeKm << " km) : "
                      << finalBest->countType(t) << " units\n";
        }

        // Export coverage overlap map and radar scatter PNG
        auto positions = extractPositions(finalBest->radars);
        calculate_radar_coverage_overlap(g, f, finalBest->radars, configs,
                                         5.0, "final_optimized_coverage.png");
        Bin::writeRadarScatterPng(positions, "final_optimized_scatter.png");
    }
};

// ============================================================
// ENTRY POINT
// ============================================================
int main() {
    ThreadPool pool;
    ga g(pool);

    // λ = 0   → free radars, pure coverage maximisation
    // λ = 50  → moderate cost; algorithm balances coverage vs. count
    // λ = 500 → heavy cost; algorithm strongly prefers fewer radars
    g.setLambda(0.0);

    g.setOperators(
        SelectionMethod::TOURNAMENT,
        CrossoverMethod::RANDOM_SUBSET,
        0.10,   // mutationRate  — 10 % of chromosomes get mutated per gen
        0.85    // crossoverRate — 85 % of pairs undergo crossover
    );

    // Fleet: each entry = { rangeKm, minCount, maxCount }
    //   "I want 2-5 short-range (50 km) and 1-3 long-range (300 km) radars"
    std::vector<RadarConfig> fleet = {
        {  50.0, 2, 5 },   // short-range
        { 300.0, 1, 3 },   // long-range
    };

    g.run_evolution(fleet, 5);
    return 0;
}
