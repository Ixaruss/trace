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
#include "radar.h"

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

    // --- Evolutionary Hyperparameters ---
    int popSize = NUM_CHROMOSOMES;
    double mutationRate = 0.08;
    double crossoverRate = 0.85;

    SelectionMethod currentSelection = SelectionMethod::TOURNAMENT;
    CrossoverMethod currentCrossover = CrossoverMethod::UNIFORM;
    MutationMethod  currentMutation  = MutationMethod::SPATIAL_NUDGE;

    std::mt19937 rng;

public:
    ga(ThreadPool& p) : pool(p) {
        std::cout << "[ga] Synchronously pre-loading all .bin matrices to memory...\n";
        g.bin.init();
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
                        // FIX 4: Explicit assignment
                        cn.radar[i] = {f.cells[x].r, f.cells[x].c};
                    }
            return cn;
        }


    static void evaluate_fitness(Global& g, const FeasibleGrid& base_grid, Chromosome& cn, int index, bool writePng) {
        std::string png_path = writePng ? "chromosome_" + std::to_string(index) + "_coverage.png" : "";
        auto res = calculate_radar_coverage_overlap(g, base_grid, cn.radar, 5.0, png_path);

        // Target: Maximize layered tracking landscape coverage, heavily fine proximity violations
        double w1 = 0.9, w2 = 0.75, w3 = 0.9, w4 = 0.9;
        double fitness = (res[1] * w1) + (res[2] * w2) - (res[3] * w3) - (res[0] * w4);
        cn.fitness = fitness;
    }

    void parallelEvaluatePopulation(std::vector<Chromosome>& pop, Global& globalRef, const FeasibleGrid& base_grid, bool writePng = false) {
        std::vector<std::future<void>> futures;
        futures.reserve(pop.size());

        for (size_t i = 0; i < pop.size(); i++) {
            futures.push_back(pool.enqueue(evaluate_fitness, std::ref(globalRef), std::ref(base_grid), std::ref(pop[i]), static_cast<int>(i), writePng));
        }

        for (auto& f : futures) f.wait();
    }

    // =========================================================================
    // EVOLUTIONARY OPERATORS (SELECTION, CROSSOVER, MUTATION)
    // =========================================================================

    // Snaps an invalid generated point back to the closest valid cell in the grid
    GridPoint snapToFeasible(int r, int c, const FeasibleGrid& f) {
            int minDistSq = 2e9;
            // FIX 1: Explicitly construct GridPoint
            GridPoint best = {f.cells.front().r, f.cells.front().c};

            for (const auto& cell : f.cells) {
                int dr = cell.r - r;
                int dc = cell.c - c;
                int distSq = dr*dr + dc*dc;
                if (distSq < minDistSq) {
                    minDistSq = distSq;
                    // FIX 2: Explicit assignment
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
        for (int i = 0; i < 2; ++i) { // Tournament size 3
            Chromosome comp = pop[dist(rng)];
            if (comp.fitness > best.fitness) best = comp;
        }
        return best;
    }

    Chromosome rouletteWheelSelection(const std::vector<Chromosome>& pop) {
        // Shift up to accommodate potential negative fitness from overlap fines
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

    // Assumes the population vector is already sorted by fitness ascending
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
        Chromosome c1(0), c2(0); // init empty

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
        std::uniform_int_distribution<int> distNudge(-5, 5); // Cell offset jump

        for (auto& gene : ind.radar) {
            if (distZeroOne(rng) < mutationRate) {
                if (currentMutation == MutationMethod::SINGLE_POINT) {
                                    // FIX 3: Explicit assignment
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
    // void init_population(int length) {
    //     population.reserve(NUM_CHROMOSOMES);

    //     std::cout << "[ga] Extracting baseline Area of Operation boundaries...\n";

    //     // Pass the loaded g.bin instance directly to prevent ghost matrix allocations
    //     auto f = Bin::getFeasibleGrid(g.bin, WRITE);

    //     for (int i = 0; i < NUM_CHROMOSOMES; i++){
    //         auto x = generate_random_chromosome(f, {0, (int)f.cells.size()-1}, length);
    //         population.push_back(x);
    //     }

    //     std::cout << "[ga] Dispatching parallel coverage evaluation and map exports...\n";
    //     parallelEvaluatePopulation(population, g, f);

    //     // Export full Generation 0 layout summary (Just red scatter dots)
    //     std::vector<GridPoint> all_radars;
    //     for (auto& cn : population)
    //         all_radars.insert(all_radars.end(), cn.radar.begin(), cn.radar.end());

    //     Bin::writeRadarScatterPng(all_radars, "population_combined_scatter.png");

    //     auto best = std::max_element(population.begin(), population.end(),
    //         [](const Chromosome& a, const Chromosome& b) { return a.fitness < b.fitness; });

    //     std::cout << "[ga] Done. Generation 0 Best Fitness Score: " << best->fitness << "\n";
    // }

    const vector<Chromosome>& get_population() const { return population; }

    void setOperators(SelectionMethod sm, CrossoverMethod cm, MutationMethod mm, double mRate, double cRate) {
        currentSelection = sm;
        currentCrossover = cm;
        currentMutation = mm;
        mutationRate = mRate;
        crossoverRate = cRate;
    }

    void run_evolution(int pathLength, int generations) {
        population.clear();
        population.reserve(popSize);

        std::cout << "[ga] Extracting baseline Area of Operation boundaries...\n";
        auto f = Bin::getFeasibleGrid(g.bin, WRITE);

        // --- Init Generation 0 ---
        for (int i = 0; i < popSize; i++) {
            population.push_back(generate_random_chromosome(f, pathLength));
        }

        std::cout << "[ga] Dispatching parallel coverage evaluation for Generation 0...\n";
        parallelEvaluatePopulation(population, g, f, false); // Keep false to avoid dumping 100 images per gen

        std::normal_distribution<double> gaussDist(0.0, 10.0); // 10 cell standard deviation

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

            // Parallel Evaluate the new breed
            parallelEvaluatePopulation(nextGen, g, f, false);
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
        calculate_radar_coverage_overlap(g, f, finalBest->radar, 5.0, "final_optimized_coverage.png");
        Bin::writeRadarScatterPng(finalBest->radar, "final_optimized_scatter.png");
    }
};

int main () {
    ThreadPool pool;
    ga g(pool);

    // Test parameters: Selection, Crossover, Mutation, M-Rate, C-Rate
    g.setOperators(SelectionMethod::TOURNAMENT, CrossoverMethod::UNIFORM, MutationMethod::SPATIAL_NUDGE, 0.08, 0.85);

    // Optimize placement of 4 radars across 50 evolutionary generations
    g.run_evolution(2, 100);

    return 0;
}
