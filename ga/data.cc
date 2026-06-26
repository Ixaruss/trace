// ================================================================================
// GA RADAR PLACEMENT  →  SAVE TOP-15 SOLUTIONS INTO HDF5
// ================================================================================
//
// Build:
//   g++ -std=c++17 -O2 ga_radar_hdf5.cpp -o ga_radar_hdf5 -lhdf5_cpp -lhdf5
//
// HDF5 output layout:
//   /config/                         ← radar & GA params as attributes
//   /top_solutions/
//       solution_0000/               ← rank-1 (best)
//           attrs : rank, chromosome_fitness, generation_found, num_radars
//           positions     (N × 2)    float32  [x, y]
//           radar_fitness (N,)       float32  per-radar score
//           radar_summary (N × 4)   float32  [x, y, radar_fitness, chrom_fitness]
//       solution_0001/ … solution_0014/
// ================================================================================

#include <H5Cpp.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <random>
#include <string>
#include <thread>
#include <vector>

// ─────────────────────────────────────────────────────────────────────────────
// CONFIG  —  swap these out for a config-file parser when ready
// ─────────────────────────────────────────────────────────────────────────────

// Radar hardware parameters
constexpr double RADAR_RANGE_PX     = 3'000.0;   // effective detection radius (px)
constexpr double ANTENNA_HGT_M      =    50.0;   // antenna height (m)
constexpr double THREAT_HGT_M       = 5'000.0;   // target threat altitude (m)
constexpr double HORIZ_AZ_DEG       =   180.0;   // horizontal beam azimuth (deg)
constexpr double VERT_EL_DEG        =    30.0;   // vertical beam elevation (deg)
constexpr int    NUM_RADARS         =    10;      // radars per solution

// Deployment grid
constexpr float  GRID_MAX           = 10'000.0f;
constexpr float  GRID_BORDER        =  2'000.0f; // minimum distance from edge

// GA parameters
constexpr int    POPULATION_SIZE    = 100;
constexpr int    NUM_GENERATIONS    = 200;
constexpr double MUTATION_RATE      =   0.02;
constexpr double CROSSOVER_RATE     =   0.80;
constexpr int    TOURNAMENT_SIZE    =   5;
constexpr int    TOP_K              =  15;        // solutions to save

// Normalisation reference (global maximum range used as denominator)
constexpr double MAX_RANGE_REF      = 5'000.0;

// HDF5 output path
const std::string HDF5_OUT_FILE     = "ga_top_solutions.h5";

constexpr double PI = 3.14159265358979323846;

// ─────────────────────────────────────────────────────────────────────────────
// DATA TYPES
// ─────────────────────────────────────────────────────────────────────────────

// Chromosome: flat [x0,y0, x1,y1, … , x(N-1),y(N-1)]
using Chromosome = std::vector<float>;

struct Individual {
    Chromosome         genes;
    double             fitness         = 0.0;
    int                generation_found = 0;
    std::vector<float> radar_fitness;           // per-radar score, length = NUM_RADARS
};

// ─────────────────────────────────────────────────────────────────────────────
// PER-RADAR FITNESS
//
// Each radar i in a chromosome gets an individual score built from four
// components whose weights sum to 1.0:
//
//   (a) range_factor   — how far this radar type can reach          × 0.35
//   (b) elev_factor    — how high it can see (sin of elevation)     × 0.20
//   (c) pos_factor     — how central the placement is on the grid   × 0.25
//   (d) isolation      — penalty for overlapping with other radars  × 0.20
//
// chromosome_fitness = mean of all per-radar scores
// ─────────────────────────────────────────────────────────────────────────────
std::vector<float> compute_per_radar_fitness(const Chromosome& genes,
                                             double&           chromosome_fitness)
{
    const int    n             = NUM_RADARS;
    const double range_score   = (RADAR_RANGE_PX / MAX_RANGE_REF) * 0.35;  // same for every radar
    const double elev_score    = std::sin(VERT_EL_DEG * PI / 180.0) * 0.20;// same for every radar

    const float  cx            = GRID_MAX * 0.5f;
    const float  cy            = GRID_MAX * 0.5f;
    const float  max_dist      = std::sqrt(cx * cx + cy * cy);              // corner-to-centre

    std::vector<float> rf(n, 0.0f);

    for (int i = 0; i < n; ++i)
    {
        const float xi = genes[i * 2 + 0];
        const float yi = genes[i * 2 + 1];

        // (c) Centralisation: radar closer to grid centre scores higher
        float  d_centre   = std::sqrt((xi - cx) * (xi - cx) + (yi - cy) * (yi - cy));
        double pos_score  = (1.0 - static_cast<double>(d_centre) / max_dist) * 0.25;

        // (d) Isolation: count other radars whose coverage circle overlaps radar i
        //     overlap defined as inter-radar distance < RADAR_RANGE_PX
        int overlaps = 0;
        for (int j = 0; j < n; ++j)
        {
            if (i == j) continue;
            float dx = genes[j * 2 + 0] - xi;
            float dy = genes[j * 2 + 1] - yi;
            if (std::sqrt(dx * dx + dy * dy) < static_cast<float>(RADAR_RANGE_PX))
                ++overlaps;
        }
        // n-1 is the maximum possible overlaps for radar i
        double isolation_score = (1.0 - static_cast<double>(overlaps) / (n - 1)) * 0.20;

        // Sum all components; clamp to [0,1]
        double total = range_score + elev_score + pos_score + isolation_score;
        rf[i] = static_cast<float>(std::clamp(total, 0.0, 1.0));
    }

    // Chromosome-level fitness = mean of individual radar scores
    double sum = 0.0;
    for (float f : rf) sum += f;
    chromosome_fitness = sum / n;

    return rf;
}

// ─────────────────────────────────────────────────────────────────────────────
// GA OPERATORS
// ─────────────────────────────────────────────────────────────────────────────

Chromosome make_random_chromosome(std::mt19937& rng)
{
    std::uniform_real_distribution<float> coord(GRID_BORDER, GRID_MAX - GRID_BORDER);
    Chromosome c(NUM_RADARS * 2);
    for (auto& v : c) v = coord(rng);
    return c;
}

// Tournament selection — returns const ref into population
const Individual& tournament_select(const std::vector<Individual>& pop,
                                    std::mt19937&                   rng)
{
    std::uniform_int_distribution<int> pick(0, static_cast<int>(pop.size()) - 1);
    int best = pick(rng);
    for (int t = 1; t < TOURNAMENT_SIZE; ++t) {
        int cand = pick(rng);
        if (pop[cand].fitness > pop[best].fitness) best = cand;
    }
    return pop[best];
}

// Single-point crossover
std::pair<Chromosome, Chromosome> crossover(const Chromosome& a,
                                             const Chromosome& b,
                                             std::mt19937&     rng)
{
    std::uniform_real_distribution<double> prob(0.0, 1.0);
    if (prob(rng) > CROSSOVER_RATE) return {a, b};

    std::uniform_int_distribution<int> cut_dist(1, static_cast<int>(a.size()) - 1);
    int pt = cut_dist(rng);

    Chromosome c1 = a, c2 = b;
    for (int i = pt; i < static_cast<int>(a.size()); ++i) std::swap(c1[i], c2[i]);
    return {std::move(c1), std::move(c2)};
}

// Gaussian mutation — perturbs coordinates in-place
void mutate(Chromosome& c, std::mt19937& rng, float sigma = 200.0f)
{
    std::uniform_real_distribution<double> prob(0.0, 1.0);
    std::normal_distribution<float>        gauss(0.0f, sigma);
    const float lo = GRID_BORDER, hi = GRID_MAX - GRID_BORDER;
    for (auto& v : c)
        if (prob(rng) < MUTATION_RATE)
            v = std::clamp(v + gauss(rng), lo, hi);
}

// ─────────────────────────────────────────────────────────────────────────────
// ELITE ARCHIVE  —  maintains the TOP_K best individuals ever evaluated
// ─────────────────────────────────────────────────────────────────────────────
class EliteArchive
{
public:
    explicit EliteArchive(int capacity) : cap_(capacity) {}

    void try_insert(const Individual& ind)
    {
        if (static_cast<int>(pool_.size()) < cap_) {
            pool_.push_back(ind);
            return;
        }
        // Replace worst if new individual is better
        auto worst = std::min_element(pool_.begin(), pool_.end(),
            [](const Individual& a, const Individual& b){
                return a.fitness < b.fitness; });
        if (ind.fitness > worst->fitness) *worst = ind;
    }

    // Return archive sorted best → worst
    std::vector<Individual> sorted() const
    {
        auto copy = pool_;
        std::sort(copy.begin(), copy.end(),
            [](const Individual& a, const Individual& b){ return a.fitness > b.fitness; });
        return copy;
    }

    double best_fitness() const
    {
        if (pool_.empty()) return 0.0;
        return std::max_element(pool_.begin(), pool_.end(),
            [](const Individual& a, const Individual& b){
                return a.fitness < b.fitness; })->fitness;
    }

private:
    int                    cap_;
    std::vector<Individual> pool_;
};

// ─────────────────────────────────────────────────────────────────────────────
// HDF5 HELPERS
// ─────────────────────────────────────────────────────────────────────────────

static void write_attr(H5::Group& g, const char* name, int value)
{
    H5::DataSpace sp(H5S_SCALAR);
    auto a = g.createAttribute(name, H5::PredType::NATIVE_INT, sp);
    a.write(H5::PredType::NATIVE_INT, &value);
}
static void write_attr(H5::Group& g, const char* name, double value)
{
    H5::DataSpace sp(H5S_SCALAR);
    auto a = g.createAttribute(name, H5::PredType::NATIVE_DOUBLE, sp);
    a.write(H5::PredType::NATIVE_DOUBLE, &value);
}
static void write_1d(H5::Group& g, const char* name,
                     const float* data, hsize_t len)
{
    H5::DataSpace sp(1, &len);
    auto ds = g.createDataSet(name, H5::PredType::NATIVE_FLOAT, sp);
    ds.write(data, H5::PredType::NATIVE_FLOAT);
}
static void write_2d(H5::Group& g, const char* name,
                     const float* data, hsize_t rows, hsize_t cols)
{
    hsize_t dims[2] = {rows, cols};
    H5::DataSpace sp(2, dims);
    auto ds = g.createDataSet(name, H5::PredType::NATIVE_FLOAT, sp);
    ds.write(data, H5::PredType::NATIVE_FLOAT);
}

// ─────────────────────────────────────────────────────────────────────────────
// MAIN
// ─────────────────────────────────────────────────────────────────────────────

int main()
{
    // ── 0. Remove stale output file ───────────────────────────────────────────
    if (std::filesystem::exists(HDF5_OUT_FILE)) {
        try {
            std::filesystem::remove(HDF5_OUT_FILE);
            std::cout << "Removed stale file: " << HDF5_OUT_FILE << "\n";
        } catch (const std::filesystem::filesystem_error&) {
            std::cout << "Warning: file locked — retrying after 500 ms...\n";
            std::this_thread::sleep_for(std::chrono::milliseconds(500));
            std::filesystem::remove(HDF5_OUT_FILE);
        }
    }

    std::mt19937 rng(42);

    // ── 1. Initialise population ──────────────────────────────────────────────
    std::cout << "\nInitialising population (" << POPULATION_SIZE << " individuals)...\n";

    std::vector<Individual> population(POPULATION_SIZE);
    EliteArchive            archive(TOP_K);

    for (auto& ind : population) {
        ind.genes          = make_random_chromosome(rng);
        ind.radar_fitness  = compute_per_radar_fitness(ind.genes, ind.fitness);
        ind.generation_found = 0;
        archive.try_insert(ind);
    }

    // ── 2. GA evolution loop ──────────────────────────────────────────────────
    std::cout << "Evolving for " << NUM_GENERATIONS << " generations...\n\n";

    for (int gen = 1; gen <= NUM_GENERATIONS; ++gen)
    {
        std::vector<Individual> next_pop;
        next_pop.reserve(POPULATION_SIZE);

        while (static_cast<int>(next_pop.size()) < POPULATION_SIZE)
        {
            // Selection → crossover → mutation → evaluate
            const auto& p1    = tournament_select(population, rng);
            const auto& p2    = tournament_select(population, rng);
            auto [g1, g2]     = crossover(p1.genes, p2.genes, rng);
            mutate(g1, rng);
            mutate(g2, rng);

            Individual i1, i2;
            i1.genes = std::move(g1);   i1.generation_found = gen;
            i2.genes = std::move(g2);   i2.generation_found = gen;
            i1.radar_fitness = compute_per_radar_fitness(i1.genes, i1.fitness);
            i2.radar_fitness = compute_per_radar_fitness(i2.genes, i2.fitness);

            archive.try_insert(i1);
            archive.try_insert(i2);

            next_pop.push_back(std::move(i1));
            if (static_cast<int>(next_pop.size()) < POPULATION_SIZE)
                next_pop.push_back(std::move(i2));
        }

        population = std::move(next_pop);

        if (gen % 50 == 0 || gen == NUM_GENERATIONS)
            std::printf("  Gen %3d  |  archive best = %.5f\n",
                        gen, archive.best_fitness());
    }

    // ── 3. Serialise top-15 into HDF5 ────────────────────────────────────────
    std::cout << "\nWriting top-" << TOP_K << " solutions → " << HDF5_OUT_FILE << "\n\n";

    const auto top_solutions = archive.sorted();

    try {
        H5::H5File h5(HDF5_OUT_FILE, H5F_ACC_TRUNC);

        // ── /config  — full parameter snapshot ─────────────────────────────
        {
            H5::Group cfg = h5.createGroup("config");

            // Radar hardware
            write_attr(cfg, "radar_range_px",  RADAR_RANGE_PX);
            write_attr(cfg, "antenna_hgt_m",   ANTENNA_HGT_M);
            write_attr(cfg, "threat_hgt_m",    THREAT_HGT_M);
            write_attr(cfg, "horiz_az_deg",    HORIZ_AZ_DEG);
            write_attr(cfg, "vert_el_deg",     VERT_EL_DEG);
            write_attr(cfg, "num_radars",      NUM_RADARS);

            // Grid
            write_attr(cfg, "grid_max",        static_cast<double>(GRID_MAX));
            write_attr(cfg, "grid_border",     static_cast<double>(GRID_BORDER));

            // GA hyper-parameters
            write_attr(cfg, "population_size", POPULATION_SIZE);
            write_attr(cfg, "num_generations", NUM_GENERATIONS);
            write_attr(cfg, "mutation_rate",   MUTATION_RATE);
            write_attr(cfg, "crossover_rate",  CROSSOVER_RATE);
            write_attr(cfg, "tournament_size", TOURNAMENT_SIZE);
            write_attr(cfg, "top_k",           TOP_K);
        }

        // ── /top_solutions  — one sub-group per elite individual ────────────
        H5::Group sol_root = h5.createGroup("top_solutions");

        for (int rank = 0; rank < static_cast<int>(top_solutions.size()); ++rank)
        {
            const Individual& ind = top_solutions[rank];

            char sol_name[20];
            std::snprintf(sol_name, sizeof(sol_name), "solution_%04d", rank);
            H5::Group sol = sol_root.createGroup(sol_name);

            // Scalar attributes
            write_attr(sol, "rank",               rank + 1);          // 1-based
            write_attr(sol, "chromosome_fitness", ind.fitness);
            write_attr(sol, "generation_found",   ind.generation_found);
            write_attr(sol, "num_radars",         NUM_RADARS);

            const hsize_t N = static_cast<hsize_t>(NUM_RADARS);

            // ── Dataset: positions  (N × 2)  [x, y] ─────────────────────────
            // genes is already laid out as [x0,y0, x1,y1, ...]
            write_2d(sol, "positions",
                     ind.genes.data(), N, 2u);

            // ── Dataset: radar_fitness  (N,)  per-radar score ────────────────
            write_1d(sol, "radar_fitness",
                     ind.radar_fitness.data(), N);

            // ── Dataset: radar_summary  (N × 4) ─────────────────────────────
            //    Columns: [ x, y, radar_fitness, chromosome_fitness ]
            std::vector<float> summary(NUM_RADARS * 4);
            for (int i = 0; i < NUM_RADARS; ++i) {
                summary[i * 4 + 0] = ind.genes[i * 2 + 0];          // x
                summary[i * 4 + 1] = ind.genes[i * 2 + 1];          // y
                summary[i * 4 + 2] = ind.radar_fitness[i];           // per-radar fitness
                summary[i * 4 + 3] = static_cast<float>(ind.fitness);// chromosome fitness
            }
            write_2d(sol, "radar_summary",
                     summary.data(), N, 4u);

            // Console summary
            std::printf("  [%2d]  fitness=%.5f  gen_found=%-4d  "
                        "radar_fitness=[",
                        rank + 1, ind.fitness, ind.generation_found);
            for (int i = 0; i < NUM_RADARS; ++i)
                std::printf("%.3f%s", ind.radar_fitness[i],
                            i < NUM_RADARS - 1 ? ", " : "");
            std::printf("]\n");
        }

        std::cout << "\nHDF5 serialisation complete → " << HDF5_OUT_FILE << "\n";
    }
    catch (const H5::Exception& e) {
        std::cerr << "HDF5 Error: " << e.getCDetailMsg() << "\n";
        return 1;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
