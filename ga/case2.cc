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

// Only include common.h for Global, Bin, Config, and FeasibleGrid structures.
#include "../common.h"

using namespace std;

// --- System & Physical Constants ---
const int NUM_CHROMOSOMES = 100;
static constexpr double EARTH_R     = 6371000.0;
static constexpr double K_ATMOS     = 4.0 / 3.0;
static constexpr double K_EFF_R     = K_ATMOS * EARTH_R;
static constexpr double INV_2_KEFFR = 1.0 / (2.0 * K_EFF_R);
static constexpr double DEG2RAD     = 3.14159265358979 / 180.0;
static constexpr double RAD2DEG     = 180.0 / 3.14159265358979;

enum class SelectionMethod { TOURNAMENT, ROULETTE_WHEEL, RANK_BASED, SUS };
enum class CrossoverMethod { SINGLE_POINT, MULTI_POINT, UNIFORM };
enum class MutationMethod  { SINGLE_POINT, SPATIAL_NUDGE, GAUSSIAN_CREEP };


// Represents a distinct radar specification requested by the user
struct RadarConfig {
    double rangeKm;
    int count;
};

struct Chromosome {
    std::vector<GridPoint> radar;
    double fitness;

    public: Chromosome(int length) : radar(length), fitness(0.0) {}
};

// --- Local Radar Coverage Container ---
struct LocalCoverage {
    std::vector<uint8_t> mask;
    int rowMin, rowMax, colMin, colMax;
    int localCols, localRows;
};

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
                        this->cv.wait(lock, [this]() { return this->stop || !this->tasks.empty(); });
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
        -> std::future<typename std::invoke_result<F, Args...>::type> {
        using return_type = typename std::invoke_result<F, Args...>::type;

        auto task = std::make_shared<std::packaged_task<return_type()>>(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...)
        );

        std::future<return_type> res = task->get_future();
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            if (stop) throw std::runtime_error("Enqueue requested on stopped context");
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
        for (std::thread &worker : workers) {
            if (worker.joinable()) worker.join();
        }
    }
};

class ga {
    Global g;
    vector<Chromosome> population;
    ThreadPool& pool;

    // --- Evolutionary Hyperparameters ---
    int popSize = NUM_CHROMOSOMES;
    double mutationRate = 0.08;
    double crossoverRate = 0.85;

    SelectionMethod currentSelection = SelectionMethod::TOURNAMENT;
    CrossoverMethod currentCrossover = CrossoverMethod::UNIFORM;
    MutationMethod  currentMutation  = MutationMethod::SPATIAL_NUDGE;

    std::mt19937 rng;

    // =========================================================================
    // EMBEDDED RADAR LOGIC (Memory-Optimized Local Scope)
    // =========================================================================

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
            [](const GridPoint& a, const GridPoint& b) { return a.r == b.r && a.c == b.c; }), pts.end());
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
        int sr = (srcRow < endRow) ? 1 : -1, sc = (srcCol < endCol) ? 1 : -1;
        int err = adr - adc;

        const std::ptrdiff_t demStrideR = static_cast<std::ptrdiff_t>(sr) * gridCols;
        const std::ptrdiff_t demStrideC = sc;
        std::ptrdiff_t demOff = static_cast<std::ptrdiff_t>(srcRow) * gridCols + srcCol;
        const short* demCell = dem + demOff;

        int localCols = colMax - colMin;
        const std::ptrdiff_t statStrideR = static_cast<std::ptrdiff_t>(sr) * localCols;
        const std::ptrdiff_t statStrideC = sc;
        std::ptrdiff_t statOff = static_cast<std::ptrdiff_t>(srcRow - rowMin) * localCols + (srcCol - colMin);
        uint8_t* statCell = stat + statOff;

        int row = srcRow, col = srcCol;
        double dist = 0.0, curvDrop = 0.0;
        double maxTerrainSlope = -1e18;
        int cellCount = 0;

        while (true) {
            if (row == endRow && col == endCol) break;

            int e2 = 2 * err;
            bool mr = (e2 > -adc), mc = (e2 < adr);
            if (mr) { err -= adc; row += sr; demCell += demStrideR; statCell += statStrideR; }
            if (mc) { err += adr; col += sc; demCell += demStrideC; statCell += statStrideC; }

            double step = (mr && mc) ? RESOLUTION * 1.41421356237 : RESOLUTION;
            curvDrop += (2.0 * dist * step + step * step) * INV_2_KEFFR;
            dist += step;
            ++cellCount;

            if (cellCount > maxRadiusCells) break;
            if (dist > horizonM) break;

            short raw = *demCell;
            if (raw < -9000) break;

            double terrainElev = static_cast<double>(raw) - curvDrop;
            double terrainSlope = (terrainElev - srcElev) / dist;
            double targetSlope = terrainSlope + targetHeight / dist;

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

        double horizonM = std::sqrt(2.0 * K_EFF_R * antennaHeight) + std::sqrt(2.0 * K_EFF_R * targetHeight);
        if (horizonM > maxRangeM) horizonM = maxRangeM;

        double minSlopeRatio = std::tan(minElevDeg * DEG2RAD);
        double maxSlopeRatio = std::tan(maxElevDeg * DEG2RAD);

        int rawRadius = (int)std::ceil(maxRangeM / RESOLUTION);
        int maxPossible = std::min({srcRow, gridRows - srcRow - 1, srcCol, gridCols - srcCol - 1});
        int maxRadiusCells = std::min(rawRadius, std::max(0, maxPossible));

        int rowMin = std::max(0, srcRow - maxRadiusCells - 1);
        int rowMax = std::min(gridRows, srcRow + maxRadiusCells + 2);
        int colMin = std::max(0, srcCol - maxRadiusCells - 1);
        int colMax = std::min(gridCols, srcCol + maxRadiusCells + 2);

        int localCols = colMax - colMin;
        int localRows = rowMax - rowMin;

        LocalCoverage cov;
        cov.rowMin = rowMin; cov.rowMax = rowMax;
        cov.colMin = colMin; cov.colMax = colMax;
        cov.localCols = localCols; cov.localRows = localRows;

        cov.mask.assign((size_t)localRows * localCols, 0);

        if (srcRow >= 0 && srcRow < gridRows && srcCol >= 0 && srcCol < gridCols) {
            short srcTerrain = g.bin.elev_matrix[(size_t)srcRow * gridCols + srcCol];
            double srcElev = static_cast<double>(srcTerrain) + antennaHeight;

            auto boundary = get_circle_boundary(srcRow, srcCol, maxRadiusCells);
            cov.mask[(size_t)(srcRow - rowMin) * localCols + (srcCol - colMin)] = 1;

            const short* demPtr = g.bin.elev_matrix.data();
            uint8_t* statPtr = cov.mask.data();

            const size_t N = boundary.size();
            for (size_t i = 0; i < N; ++i) {
                int er = boundary[i].r, ec = boundary[i].c;
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

    static std::vector<int> calculate_radar_coverage_overlap(
        Global& g, FeasibleGrid grid, const std::vector<GridPoint>& placements,
        const std::vector<double>& radarRanges,
        double minSpacingCells = 5.0, const std::string& pngPath = "")
    {
        std::vector<int> result(4, 0);
        if (placements.empty() || grid.cells.empty()) return result;

        for (size_t i = 0; i < placements.size(); i++) {
            for (size_t j = i + 1; j < placements.size(); j++) {
                double dr = placements[i].r - placements[j].r;
                double dc = placements[i].c - placements[j].c;
                if (std::sqrt(dr*dr + dc*dc) < minSpacingCells) result[0]++;
            }
        }

        // Evaluate coverage dynamically based on the range mapped to this specific index locus
        for (size_t i = 0; i < placements.size(); i++) {
            const auto& p = placements[i];
            double maxRangeM = radarRanges[i];

            LocalCoverage cov = compute_local_coverage(
                g, p.r, p.c, 50.0, 1000.0, 10.0, 30.0, maxRangeM
            );

            for (auto& cell : grid.cells) {
                if (cell.r >= cov.rowMin && cell.r < cov.rowMax &&
                    cell.c >= cov.colMin && cell.c < cov.colMax) {
                    size_t local_r = cell.r - cov.rowMin;
                    size_t local_c = cell.c - cov.colMin;
                    size_t flat = local_r * cov.localCols + local_c;

                    if (flat < cov.mask.size() && cov.mask[flat] == 1) {
                        cell.overlap++;
                    }
                }
            }
        }

        for (const auto& cell : grid.cells) {
            if      (cell.overlap == 1) result[1]++;
            else if (cell.overlap == 2) result[2]++;
            else if (cell.overlap >= 3) result[3]++;
        }

        if (!pngPath.empty()) Bin::writeFeasibleGridOverlapPng(grid, pngPath,placements);

        return result;
    }

    static void evaluate_fitness(Global& g, const FeasibleGrid& base_grid, Chromosome& cn, int index, const std::vector<double>& radarRanges, bool writePng) {
        std::string png_path = writePng ? "chromosome_" + std::to_string(index) + "_coverage.png" : "";
        auto res = calculate_radar_coverage_overlap(g, base_grid, cn.radar, radarRanges, 5.0, png_path);

        // Target: Maximize layered tracking landscape coverage, heavily fine proximity violations
        cn.fitness = (res[1] * 1.0) + (res[2] * 2.0) + (res[3] * 3.0) - (res[0] * 5000.0);
    }

    void parallelEvaluatePopulation(std::vector<Chromosome>& pop, Global& globalRef, const FeasibleGrid& base_grid, const std::vector<double>& radarRanges, bool writePng = false) {
        std::vector<std::future<void>> futures;
        futures.reserve(pop.size());

        for (size_t i = 0; i < pop.size(); i++) {
            futures.push_back(pool.enqueue(evaluate_fitness, std::ref(globalRef), std::ref(base_grid), std::ref(pop[i]), static_cast<int>(i), std::ref(radarRanges), writePng));
        }

        for (auto& f : futures) f.wait();
    }

    // =========================================================================
    // EVOLUTIONARY OPERATORS (SELECTION, CROSSOVER, MUTATION)
    // =========================================================================

    // Snaps an invalid generated point back to the closest valid cell in the grid
    GridPoint snapToFeasible(int r, int c, const FeasibleGrid& f) {
        int minDistSq = 2e9;
        GridPoint best = {f.cells.front().r, f.cells.front().c};
        for (const auto& cell : f.cells) {
            int dr = cell.r - r;
            int dc = cell.c - c;
            int distSq = dr*dr + dc*dc;
            if (distSq < minDistSq) {
                minDistSq = distSq;
                best = {cell.r, cell.c};
                if (distSq == 0) break; // Exact match found
            }
        }
        return best;
    }

    // --- SELECTION ---
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

        double totalFitness = 0.0;
        for (const auto& ind : pop) totalFitness += (ind.fitness + offset);

        std::uniform_real_distribution<double> dist(0.0, totalFitness);
        double slice = dist(rng);
        double currentSum = 0.0;

        for (const auto& ind : pop) {
            currentSum += (ind.fitness + offset);
            if (currentSum >= slice) return ind;
        }
        return pop.back();
    }

    Chromosome rankBasedSelection(const std::vector<Chromosome>& sorted_pop) {
        int n = sorted_pop.size();
        int totalRankSum = (n * (n + 1)) / 2;
        std::uniform_int_distribution<int> dist(1, totalRankSum);
        int target = dist(rng);
        int rankSum = 0;

        for (int i = 0; i < n; ++i) {
            rankSum += (i + 1);
            if (rankSum >= target) return sorted_pop[i];
        }
        return sorted_pop.back();
    }

    std::vector<Chromosome> stochasticUniversalSampling(const std::vector<Chromosome>& pop, int numToSelect) {
        std::vector<Chromosome> selected;

        double minFit = pop[0].fitness;
        for (const auto& ind : pop) if (ind.fitness < minFit) minFit = ind.fitness;
        double offset = (minFit < 0) ? std::abs(minFit) + 1.0 : 0.0;

        double totalFitness = 0.0;
        for (const auto& ind : pop) totalFitness += (ind.fitness + offset);

        double pointerDistance = totalFitness / numToSelect;
        std::uniform_real_distribution<double> dist(0.0, pointerDistance);
        double startPointer = dist(rng);

        double currentSum = 0.0;
        size_t idx = 0;

        for (int i = 0; i < numToSelect; ++i) {
            double currentPointer = startPointer + (i * pointerDistance);
            while (currentSum + (pop[idx].fitness + offset) < currentPointer && idx < pop.size() - 1) {
                currentSum += (pop[idx].fitness + offset);
                idx++;
            }
            selected.push_back(pop[idx]);
        }
        return selected;
    }

    Chromosome executeSelection(const std::vector<Chromosome>& pop) {
        if (currentSelection == SelectionMethod::TOURNAMENT) return tournamentSelection(pop);
        if (currentSelection == SelectionMethod::ROULETTE_WHEEL) return rouletteWheelSelection(pop);
        return rankBasedSelection(pop);
    }

    // --- CROSSOVER ---
    std::pair<Chromosome, Chromosome> executeCrossover(const Chromosome& p1, const Chromosome& p2) {
        std::uniform_real_distribution<double> distZeroOne(0.0, 1.0);
        if (distZeroOne(rng) > crossoverRate) return {p1, p2};

        int len = p1.radar.size();
        Chromosome c1(0), c2(0);

        if (currentCrossover == CrossoverMethod::SINGLE_POINT) {
            std::uniform_int_distribution<int> distPoint(1, len - 2);
            int cp = distPoint(rng);
            c1.radar.insert(c1.radar.end(), p1.radar.begin(), p1.radar.begin() + cp);
            c1.radar.insert(c1.radar.end(), p2.radar.begin() + cp, p2.radar.end());
            c2.radar.insert(c2.radar.end(), p2.radar.begin(), p2.radar.begin() + cp);
            c2.radar.insert(c2.radar.end(), p1.radar.begin() + cp, p1.radar.end());
        }
        else if (currentCrossover == CrossoverMethod::MULTI_POINT) {
            std::uniform_int_distribution<int> distPoint(1, len - 2);
            int cp1 = distPoint(rng), cp2 = distPoint(rng);
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
                if (distZeroOne(rng) < 0.5) {
                    c1.radar.push_back(p1.radar[i]); c2.radar.push_back(p2.radar[i]);
                } else {
                    c1.radar.push_back(p2.radar[i]); c2.radar.push_back(p1.radar[i]);
                }
            }
        }
        return {c1, c2};
    }

    // --- MUTATION ---
    void executeMutation(Chromosome& ind, const FeasibleGrid& f, std::normal_distribution<double>& gauss) {
        std::uniform_real_distribution<double> distZeroOne(0.0, 1.0);
        std::uniform_int_distribution<int> distF(0, f.cells.size() - 1);
        std::uniform_int_distribution<int> distNudge(-5, 5);

        for (auto& gene : ind.radar) {
            if (distZeroOne(rng) < mutationRate) {
                if (currentMutation == MutationMethod::SINGLE_POINT) {
                    gene = {f.cells[distF(rng)].r, f.cells[distF(rng)].c};
                }
                else if (currentMutation == MutationMethod::SPATIAL_NUDGE) {
                    gene = snapToFeasible(gene.r + distNudge(rng), gene.c + distNudge(rng), f);
                }
                else { // GAUSSIAN_CREEP
                    int rOff = static_cast<int>(std::round(gauss(rng)));
                    int cOff = static_cast<int>(std::round(gauss(rng)));
                    gene = snapToFeasible(gene.r + rOff, gene.c + cOff, f);
                }
            }
        }
    }

public:
    ga(ThreadPool& p) : pool(p) {
        rng.seed(static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count()));
        std::cout << "[ga] Synchronously pre-loading all .bin matrices to memory...\n";
        g.bin.init();
    }

    void setOperators(SelectionMethod sm, CrossoverMethod cm, MutationMethod mm, double mRate, double cRate) {
        currentSelection = sm;
        currentCrossover = cm;
        currentMutation = mm;
        mutationRate = mRate;
        crossoverRate = cRate;
    }

    Chromosome generate_random_chromosome(const FeasibleGrid& f, int length) {
        int max_cells = f.cells.size();
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

    void run_evolution(const std::vector<RadarConfig>& fleetConfig, int generations) {
        // Flatten the configured fleet into a strict index array mapped directly to meters
        std::vector<double> flattenedRangesM;
        for (const auto& config : fleetConfig) {
            for (int i = 0; i < config.count; ++i) {
                flattenedRangesM.push_back(config.rangeKm * 1000.0);
            }
        }

        int pathLength = flattenedRangesM.size();

        population.clear();
        population.reserve(popSize);

        std::cout << "[ga] Extracting baseline Area of Operation boundaries...\n";
        auto f = Bin::getFeasibleGrid(g.bin, WRITE);

        // --- Init Generation 0 ---
        for (int i = 0; i < popSize; i++) {
            population.push_back(generate_random_chromosome(f, pathLength));
        }

        std::cout << "[ga] Dispatching parallel coverage evaluation for Generation 0...\n";
        parallelEvaluatePopulation(population, g, f, flattenedRangesM, false);

        std::normal_distribution<double> gaussDist(0.0, 10.0);

        for (int gen = 1; gen <= generations; ++gen) {

            // Sort population (vital for Rank-Based Selection and Elitism extraction)
            std::sort(population.begin(), population.end(), [](const Chromosome& a, const Chromosome& b) {
                return a.fitness < b.fitness;
            });

            std::vector<Chromosome> nextGen;
            nextGen.reserve(popSize);

            // Carry over the absolute best to prevent regression
            nextGen.push_back(population.back());

            bool useSUS = (currentSelection == SelectionMethod::SUS);
            std::vector<Chromosome> matingPool;
            if (useSUS) {
                matingPool = stochasticUniversalSampling(population, popSize - 1);
            }

            size_t poolIdx = 0;
            while (nextGen.size() < (size_t)popSize) {
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
                if (nextGen.size() < (size_t)popSize) nextGen.push_back(c2);
            }

            // Parallel Evaluate the new breed using the strict locus mapping array
            parallelEvaluatePopulation(nextGen, g, f, flattenedRangesM, false);
            population = std::move(nextGen);

            auto bestInd = std::max_element(population.begin(), population.end(),
                [](const Chromosome& a, const Chromosome& b) { return a.fitness < b.fitness; });

            std::cout << "  -> Gen " << gen << " Best Local Fitness Score: " << bestInd->fitness << "\n";
        }

        // --- Post-Evolution Cleanup & Export ---
        auto finalBest = std::max_element(population.begin(), population.end(),
            [](const Chromosome& a, const Chromosome& b) { return a.fitness < b.fitness; });

        std::cout << "\n[ga] Evolution Complete. Best Lifetime Score: " << finalBest->fitness << "\n";

        // Dump visual matrices for the final absolute best setup
        calculate_radar_coverage_overlap(g, f, finalBest->radar, flattenedRangesM, 5.0, "final_optimized_coverage.png");
        Bin::writeRadarScatterPng(finalBest->radar, "final_optimized_scatter.png");
    }
};

int main () {
    ThreadPool pool;
    ga g(pool);

    // Test parameters: Selection, Crossover, Mutation, M-Rate, C-Rate
    g.setOperators(SelectionMethod::TOURNAMENT, CrossoverMethod::UNIFORM, MutationMethod::SPATIAL_NUDGE, 0.08, 0.85);


    std::vector<RadarConfig> fleet = {
        {50.0, 2},
        {150.0, 1},
        {300.0, 1}
    };
    g.run_evolution(fleet, 50);

    return 0;
}
