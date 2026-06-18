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

// Only include common.h for Global, Bin, Config, and FeasibleGrid structures.
// operations.h is intentionally omitted to keep radar calculations localized and memory-safe.
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

struct GridPoint {
    int r, c;
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
        for (std::thread &worker : workers) {
            if (worker.joinable()) worker.join();
        }
    }
};

class ga {
    Global g;
    vector<Chromosome> population;
    ThreadPool& pool;

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

        // Global DEM traversal pointers
        const std::ptrdiff_t demStrideR = static_cast<std::ptrdiff_t>(sr) * gridCols;
        const std::ptrdiff_t demStrideC = sc;
        std::ptrdiff_t demOff = static_cast<std::ptrdiff_t>(srcRow) * gridCols + srcCol;
        const short* demCell = dem + demOff;

        // Local bounding-box mask traversal pointers
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

                    // CRITICAL FIX:
                    // If current terrain slope is HIGHER than our established maxSlope,
                    // we have hit a ridge/mountain. This cell is visible, but EVERYTHING
                    // behind this cell on this ray is now MASKED.
                    if (terrainSlope > maxTerrainSlope) {
                        maxTerrainSlope = terrainSlope;
                        *statCell = 1; // This cell is visible
                    } else {
                        // We are in the shadow of a previous higher terrain point
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

        // Clamp radius tightly to valid grid extents to prevent buffer escapes
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

        // This instantly assigns ONLY the bounding box memory (~4MB - ~72MB depending on range),
        // completely resolving the massive 10+ GB memory spike from the old global allocation.
        cov.mask.assign((size_t)localRows * localCols, 0);

        if (srcRow >= 0 && srcRow < gridRows && srcCol >= 0 && srcCol < gridCols) {
            short srcTerrain = g.bin.elev_matrix[(size_t)srcRow * gridCols + srcCol];
            double srcElev = static_cast<double>(srcTerrain) + antennaHeight;

            auto boundary = get_circle_boundary(srcRow, srcCol, maxRadiusCells);

            // Self-cell is automatically visible
            cov.mask[(size_t)(srcRow - rowMin) * localCols + (srcCol - colMin)] = 1;

            const short* demPtr = g.bin.elev_matrix.data();
            uint8_t* statPtr = cov.mask.data();

            // Fire standard raycaster sequentially (ThreadPool handles the macro-parallelism)
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

    // =========================================================================

public:
    ga(ThreadPool& p) : pool(p) {
        std::cout << "[ga] Synchronously pre-loading all .bin matrices to memory...\n";
        g.bin.init();
    }

    static Chromosome generate_random_chromosome(const FeasibleGrid& f, pair<int,int> range, int length) {
        if ((size_t)length > f.cells.size()) {
            std::cerr << "[ga] WARN: requested " << length << " radars but only "
                      << f.cells.size() << " feasible cells exist — clamping\n";
            length = (int)f.cells.size();
        }

        random_device rd;
        mt19937 code(rd());
        uniform_int_distribution<int> gen(range.first, range.second);

        Chromosome cn(length);
        std::unordered_set<int> used;

        for (int i = 0; i < length; i++) {
            int x;
            do { x = gen(code); } while (!used.insert(x).second);

            cn.radar[i].r = f.cells[x].r;
            cn.radar[i].c = f.cells[x].c;
        }

        cn.fitness = 0.0;
        return cn;
    }

    static std::vector<int> calculate_radar_coverage_overlap(
        Global&                       g,
        FeasibleGrid                  grid,              // BY VALUE — this thread's private copy
        const std::vector<GridPoint>& placements,
        double                        minSpacingCells = 5.0,
        const std::string&            pngPath = "")
    {
        std::vector<int> result(4, 0);

        if (placements.empty() || grid.cells.empty()) return result;

        // --- [0] spacing gaps: pairwise under-separation within this chromosome ---
        for (size_t i = 0; i < placements.size(); i++) {
            for (size_t j = i + 1; j < placements.size(); j++) {
                double dr = placements[i].r - placements[j].r;
                double dc = placements[i].c - placements[j].c;
                if (std::sqrt(dr*dr + dc*dc) < minSpacingCells) result[0]++;
            }
        }

        // --- coverage overlap per feasible cell ---
        for (const auto& p : placements) {
            // Radar specifications: 50m AGL, 1000m tgt, 10-30 deg Elev, 300km Range
            LocalCoverage cov = compute_local_coverage(
                g, p.r, p.c,
                50.0, 1000.0, 10.0, 30.0, 300000.0
            );

            // Shift grid cell coordinates into local radar bounds
            for (auto& cell : grid.cells) {
                if (cell.r >= cov.rowMin && cell.r < cov.rowMax &&
                    cell.c >= cov.colMin && cell.c < cov.colMax) {

                    size_t local_r = cell.r - cov.rowMin;
                    size_t local_c = cell.c - cov.colMin;
                    size_t flat = local_r * cov.localCols + local_c;

                    if (flat < cov.mask.size()) {
                        // Only add to overlap if the radar actually covers this specific cell
                        if (cov.mask[flat] == 1) {
                            cell.overlap++;
                        }
                    }
                }
            }
        }

        // --- bucket final per-cell overlap counts ---
        for (const auto& cell : grid.cells) {
            if      (cell.overlap == 1) result[1]++;
            else if (cell.overlap == 2) result[2]++;
            else if (cell.overlap >= 3) result[3]++;
        }

        // --- Render Visual Integrity Map (Uses overlapToColor from bin.h natively) ---
        if (!pngPath.empty()) {
            Bin::writeFeasibleGridOverlapPng(grid, pngPath);
        }

        return result;
    }

    static void evaluate_fitness(Global& g, const FeasibleGrid& base_grid, Chromosome& cn, int index) {
        std::string png_path = "chromosome_" + std::to_string(index) + "_coverage.png";

        auto res = calculate_radar_coverage_overlap(g, base_grid, cn.radar, 5.0, png_path);

        // Target: Maximize layered tracking landscape coverage, heavily fine proximity violations
        double fitness = (res[1] * 1.0) + (res[2] * 2.0) + (res[3] * 3.0) - (res[0] * 5000.0);

        cn.fitness = fitness;
    }

    void parallelEvaluatePopulation(std::vector<Chromosome>& pop, Global& g, const FeasibleGrid& base_grid) {
        std::vector<std::future<void>> futures;
        futures.reserve(pop.size());

        for (size_t i = 0; i < pop.size(); i++) {
            futures.push_back(pool.enqueue(evaluate_fitness, std::ref(g), std::ref(base_grid), std::ref(pop[i]), static_cast<int>(i)));
        }

        for (auto& f : futures) f.wait();
    }

    void init_population(int length) {
        population.reserve(NUM_CHROMOSOMES);

        std::cout << "[ga] Extracting baseline Area of Operation boundaries...\n";

        // Pass the loaded g.bin instance directly to prevent ghost matrix allocations
        auto f = Bin::getFeasibleGrid(g.bin, WRITE);

        for (int i = 0; i < NUM_CHROMOSOMES; i++){
            auto x = generate_random_chromosome(f, {0, (int)f.cells.size()-1}, length);
            population.push_back(x);
        }

        std::cout << "[ga] Dispatching parallel coverage evaluation and map exports...\n";
        parallelEvaluatePopulation(population, g, f);

        // Export full Generation 0 layout summary (Just red scatter dots)
        std::vector<GridPoint> all_radars;
        for (auto& cn : population)
            all_radars.insert(all_radars.end(), cn.radar.begin(), cn.radar.end());

        Bin::writeRadarScatterPng(all_radars, "population_combined_scatter.png");

        auto best = std::max_element(population.begin(), population.end(),
            [](const Chromosome& a, const Chromosome& b) { return a.fitness < b.fitness; });

        std::cout << "[ga] Done. Generation 0 Best Fitness Score: " << best->fitness << "\n";
    }

    const vector<Chromosome>& get_population() const { return population; }
};

int main () {
    ThreadPool pool;
    ga g(pool);
    g.init_population(4); // Generates 4 radars per chromosome

    return 0;
}
