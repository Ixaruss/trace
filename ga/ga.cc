#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <chrono>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <functional>
#include <iomanip>
#include <string>
#include <fstream>

#include "../common.h"
#include "../operations.h"

using namespace std;

const int NUM_CHROMOSOMES = 100;
const double EARTH_RADIUS_KM = 6371.0;
const double PI = 3.14159265358979323846;

enum class SelectionMethod { TOURNAMENT, ROULETTE_WHEEL, RANK_BASED, SUS };
enum class CrossoverMethod { SINGLE_POINT, MULTI_POINT, UNIFORM };
enum class MutationMethod  { SINGLE_POINT, SPATIAL_NUDGE, GAUSSIAN_CREEP };


struct GridPoint {
    int r, c;
};

struct Chromosome {
    std::vector<GridPoint> radar;
    double fitness;

    public: Chromosome(int length){
        std::vector<GridPoint> r(length);
        radar = r;
    }
};

class ThreadPool {
private:
    std::vector<std::thread> workers;       // Active system threads waiting for work
    std::queue<std::function<void()>> tasks;// Thread-safe FIFO queue containing pending tasks
    std::mutex queueMutex;                  // Mutex locks ensuring unique access to the queue
    std::condition_variable cv;             // Blocks or wakes worker threads dynamically
    bool stop;                              // Execution termination flag

public:
    /**
     * Initializes the worker pool based on native hardware concurrency capabilities.
     */
    ThreadPool() : stop(false) {
        // Detect available processor cores
        size_t threads = std::thread::hardware_concurrency();
        if (threads == 0) threads = 4; // Defensive fallback configuration if detection fails

        for (size_t i = 0; i < threads; ++i) {
            // Emplace active lambdas acting as continuous internal processing engines
            workers.emplace_back([this]() {
                while (true) {
                    std::function<void()> task;
                    {
                        // Safely lock the queue to verify state changes
                        std::unique_lock<std::mutex> lock(this->queueMutex);

                        // Sleep thread until new tasks are queued or system shuts down
                        this->cv.wait(lock, [this]() { return this->stop || !this->tasks.empty(); });

                        // Exit worker execution loop if termination is requested and queue is drained
                        if (this->stop && this->tasks.empty()) return;

                        // Extract task from front of queue via efficient move semantics
                        task = std::move(this->tasks.front());
                        this->tasks.pop();
                    } // Mutex is implicitly unlocked here to allow other threads to read the queue

                    task(); // Execute the extracted task
                }
            });
        }
    }

    /**
     * Schedules a task asynchronously. Accepts any callable function and arguments.
     * Returns a std::future containing the eventual execution results.
     */
    template<class F, class... Args>
    auto enqueue(F&& f, Args&&... args)
        -> std::future<typename std::invoke_result<F, Args...>::type> {
        using return_type = typename std::invoke_result<F, Args...>::type;

        // Package the function pointer alongside its arguments into a shared reference block
        auto task = std::make_shared<std::packaged_task<return_type()>>(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...)
        );

        std::future<return_type> res = task->get_future();
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            if (stop) throw std::runtime_error("Enqueue requested on stopped ThreadPool context");

            // Wrap the packaged task execution inside a generalized void wrapper
            tasks.emplace([task]() { (*task)(); });
        }
        cv.notify_one(); // Wake up exactly one sleeping thread to execute the new task
        return res;      // Return future context to caller
    }

    /**
     * Destructor: Ensures graceful shutdowns, joining all threads safely.
     */
    ~ThreadPool() {
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            stop = true;
        }
        cv.notify_all(); // Wake up all worker threads to finalize execution loops

        for (std::thread &worker : workers) {
            if (worker.joinable()) worker.join(); // Block main execution loop until thread exits
        }
    }
};

class ga {
    Global g;
    vector<Chromosome> population;
    ThreadPool& pool;

public:
    // FIX: no constructor existed, so the `ThreadPool& pool` reference member
    // could never be initialized — this class wasn't constructible at all.
    ga(ThreadPool& p) : pool(p) {}

    static Chromosome generate_random_chromosome(FeasibleGrid f,pair<int,int> range, int length){
        random_device rd;
        Chromosome cn(length);
        mt19937 code(rd());
        uniform_int_distribution<int> gen(range.first, range.second);
        for (int i=0; i<length; i++) {
            // int x = gen(code);
            // cn.radar[i].r = f.cells[x].r;
            // cn.radar[i].c = f.cells[x].c;
            // FIX: this called gen(code) twice more here, so it printed
            // different random numbers than what actually got stored above.
            cout << "{" << cn.radar[i].r << "," << cn.radar[i].c << "}" << endl;
        }
        cout<<endl;
        cn.fitness = 0.0;
        return cn;
    }

    static std::vector<int> calculate_radar_coverage_overlap(
        Global&                 g,
        FeasibleGrid             grid,             // BY VALUE — this chromosome's private copy
        const std::vector<int>&  radarIndices,      // genes: indices into grid.cells
        double                   minSpacingCells = 5.0,
        const std::string&       pngPath = "")
    {
        std::vector<int> result(4, 0);

        if (radarIndices.empty() || grid.cells.empty()) return result;

        std::vector<GridPoint> placements;
        placements.reserve(radarIndices.size());
        for (int idx : radarIndices) {
            if (idx < 0 || (size_t)idx >= grid.cells.size()) continue; // guard a stale/out-of-range gene
            placements.push_back({grid.cells[idx].r, grid.cells[idx].c});
        }
        if (placements.empty()) return result;

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
            // construction lat/lon is a throwaway — position is overwritten
            // below using the gene's actual (r,c)
            RadarEngine radar(g, 33.9871, 74.7742, 50.0, 1000.0, 10.0, 30.0, 300.0);

            radar.srcRow = p.r;
            radar.srcCol = p.c;

            // maxRadiusCells is also position-clamped in the real constructor —
            // re-derive it here too, or a radar near an edge keeps whatever
            // radius the throwaway construction position happened to allow.
            int rawRadius   = (int)std::ceil(radar.maxRangeM / RESOLUTION);
            int maxPossible = std::min(std::min(p.r, radar.gridRows - p.r - 1),
                                        std::min(p.c, radar.gridCols - p.c - 1));
            radar.maxRadiusCells = std::min(rawRadius, std::max(0, maxPossible));

            radar.rowMin = std::max(0, p.r - radar.maxRadiusCells - 1);
            radar.rowMax = std::min(radar.gridRows, p.r + radar.maxRadiusCells + 2);
            radar.colMin = std::max(0, p.c - radar.maxRadiusCells - 1);
            radar.colMax = std::min(radar.gridCols, p.c + radar.maxRadiusCells + 2);
            short srcTerrain = g.bin.elev_matrix[(size_t)p.r * radar.gridCols + p.c];
            radar.srcElev = (double)srcTerrain + radar.antennaHeight;
            radar.boundary = RadarEngine::circleBoundary(p.r, p.c, radar.maxRadiusCells);

            std::vector<uint8_t> mask = radar.computeFast();

            for (auto& cell : grid.cells) {
                size_t flat = (size_t)cell.r * radar.gridCols + cell.c;
                if (flat < mask.size() && mask[flat])
                    cell.overlap++;
            }
        }

        // --- bucket final per-cell overlap counts ---
        for (const auto& cell : grid.cells) {
            if      (cell.overlap == 1) result[1]++;
            else if (cell.overlap == 2) result[2]++;
            else if (cell.overlap >= 3) result[3]++;
        }

        if (!pngPath.empty())
            Bin::writeFeasibleGridOverlapPng(grid, pngPath);

        return result;
    }
    // FIX: Global was passed by value. Global/Bin owns a raw-pointer handle
    // with a custom destructor and no copy ctor, so every call silently
    // shallow-copied it and double-destroyed the handle on return (confirmed
    // above). const& avoids the copy entirely. If RadarEngine's constructor
    // needs non-const access, change this to `Global&` and use std::ref(g)
    // instead of cref(g) at the call site below.
    static void evaluate_fitness(const Global& g, Chromosome &cn){
        vector<int> overlaps(4);
        RadarEngine radar(
                g,
                33.9871, 74.7742, // Srinagar airport
                50.0,    // antenna height AGL (m)
                1000.0,  // target height AGL (m) — e.g. aircraft at 1km
                10.0,     // min elevation angle (degrees)
                30.0,    // max elevation angle (degrees)
                300.0    // max range (km)
        );
        for (size_t k = 0; k < cn.radar.size(); k++){
            calculate_radar_coverage_overlap(&g,f, const std::vector<int> &radarIndices)
        }
        // FIX: w1..w4 had no initializer and were used immediately — reading
        // uninitialized floats is undefined behavior. Placeholder 1.0s so
        // this is at least deterministic; actual weights are a tuning call
        // for you, not something I'm guessing at here.
        float w1 = 1.0f, w2 = 1.0f, w3 = 1.0f, w4 = 1.0f;
        double fitness = w1 * overlaps[1] + w2 * overlaps[2] - w3 * overlaps[3] - w4 * overlaps[0];
        cn.fitness = fitness;
    }

    // FIX: matches the const& change above.
    void parallelEvaluatePopulation(std::vector<Chromosome>& population, const Global& g) {
        std::vector<std::future<void>> futures;
        futures.reserve(population.size());

        for (auto& ind : population) {
            // FIX: args were swapped here (ind then g) against
            // evaluate_fitness(Global, Chromosome&) — wouldn't compile.
            futures.push_back(pool.enqueue(evaluate_fitness, cref(g), std::ref(ind)));
        }

        for (auto& f : futures) f.wait();
    }

    // FIX: was declared to return vector<Chromosome> but never returned
    // anything (UB), and nothing used the return value anyway — population
    // is already a member, so this is void now.
    void init_population(int length) {
        population.reserve(NUM_CHROMOSOMES);
        auto f = Bin::getFeasibleGrid(WRITE);
        for (int i=0; i < NUM_CHROMOSOMES; i++){
            // FIX: was `population[i] = x` after only reserve() (no resize) —
            // writing into unconstructed memory (confirmed segfault above).
            // push_back is the correct way to grow into reserved capacity.
            auto x = generate_random_chromosome(f,{0,(int)f.cells.size()-1}, length);
            population.push_back(x);
        }
        parallelEvaluatePopulation(population, g);

        // --- visualization (new) ---
        // one PNG per chromosome: just its radar placements, red dots, no
        // coverage/overlap rendering
        for (size_t i = 0; i < population.size(); i++) {
            Bin::writeRadarScatterPng(population[i].radar,
                                       "chromosome_" + std::to_string(i) + ".png");
        }

        // combined view: every chromosome's radar positions on one shared
        // canvas — still just presence (red), not coverage
        std::vector<GridPoint> all_radars;
        for (auto& cn : population)
            all_radars.insert(all_radars.end(), cn.radar.begin(), cn.radar.end());
        Bin::writeRadarScatterPng(all_radars, "population_combined.png");
    }

    const vector<Chromosome>& get_population() const { return population; }
};

int main () {
    ThreadPool pool;
    ga g(pool);
    g.init_population(/*length=*/4);   // radars per chromosome
}
