/**
 * fitness_writer.h
 *
 * Computes per-placement individual radar fitness and writes results
 * plus radar parameters to an HDF5 file.
 *
 * HDF5 layout
 * ───────────
 *   /<run_name>/            ← unique, memorable run group (returned by write_h5_indiviual)
 *       params/
 *           antenna_height  scalar f64
 *           target_height   scalar f64
 *           min_angle_deg   scalar f64
 *           max_angle_deg   scalar f64
 *           max_range_m     scalar f64
 *       data/
 *           r               i32[N]   row    indices of placements (GridPoint::r)
 *           c               i32[N]   column indices of placements (GridPoint::c)
 *           fitness        f64[N]   actual / ideal feasible-cell ratio ∈ [0, 1]
 *
 * Build: link with  -lhdf5
 * Dependencies: radar_overlap.h  (get_individual_fitness, Global, FeasibleGrid)
 */
#pragma once

#include <random>
#include <vector>
#include <string>
#include <stdexcept>
#include <cmath>

#include <hdf5.h>

#include "radar.h"   // get_individual_fitness, Global, FeasibleGrid

#if defined(_OPENMP)
#  include <omp.h>
#endif


// ═══════════════════════════════════════════════════════════════════════════
// Internal helpers
// ═══════════════════════════════════════════════════════════════════════════
namespace h5_lib {

// ── Lightweight RAII wrappers — prevent HID leaks on exceptions ───────────
struct H5FileGuard  { hid_t id; ~H5FileGuard()  { if (id >= 0) H5Fclose(id); } };
struct H5GroupGuard { hid_t id; ~H5GroupGuard() { if (id >= 0) H5Gclose(id); } };
struct H5SpaceGuard { hid_t id; ~H5SpaceGuard() { if (id >= 0) H5Sclose(id); } };
struct H5DSetGuard  { hid_t id; ~H5DSetGuard()  { if (id >= 0) H5Dclose(id); } };

// ── Write a named scalar f64 dataset into an open group ───────────────────
inline void write_scalar(hid_t group, const char* name, double value)
{
    H5SpaceGuard sid{ H5Screate(H5S_SCALAR) };
    if (sid.id < 0)
        throw std::runtime_error(std::string("HDF5: H5Screate failed for ") + name);

    H5DSetGuard did{ H5Dcreate2(group, name, H5T_IEEE_F64LE,
                                 sid.id,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) };
    if (did.id < 0)
        throw std::runtime_error(std::string("HDF5: H5Dcreate2 failed for ") + name);

    if (H5Dwrite(did.id, H5T_NATIVE_DOUBLE,
                 H5S_ALL, H5S_ALL, H5P_DEFAULT, &value) < 0)
        throw std::runtime_error(std::string("HDF5: H5Dwrite failed for ") + name);
}

// ── Write a named 1-D i32 dataset into an open group ─────────────────────
inline void write_i32(hid_t group, const char* name,
                      const std::vector<int>& v)
{
    hsize_t dims[1] = { v.size() };
    H5SpaceGuard sid{ H5Screate_simple(1, dims, nullptr) };
    if (sid.id < 0)
        throw std::runtime_error(std::string("HDF5: H5Screate_simple failed for ") + name);

    H5DSetGuard did{ H5Dcreate2(group, name, H5T_STD_I32LE,
                                 sid.id,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) };
    if (did.id < 0)
        throw std::runtime_error(std::string("HDF5: H5Dcreate2 failed for ") + name);

    if (H5Dwrite(did.id, H5T_NATIVE_INT,
                 H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data()) < 0)
        throw std::runtime_error(std::string("HDF5: H5Dwrite failed for ") + name);
}

// ── Write a named 1-D f64 dataset into an open group ─────────────────────
inline void write_f64(hid_t group, const char* name,
                      const std::vector<double>& v)
{
    hsize_t dims[1] = { v.size() };
    H5SpaceGuard sid{ H5Screate_simple(1, dims, nullptr) };
    if (sid.id < 0)
        throw std::runtime_error(std::string("HDF5: H5Screate_simple failed for ") + name);

    H5DSetGuard did{ H5Dcreate2(group, name, H5T_IEEE_F64LE,
                                 sid.id,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) };
    if (did.id < 0)
        throw std::runtime_error(std::string("HDF5: H5Dcreate2 failed for ") + name);

    if (H5Dwrite(did.id, H5T_NATIVE_DOUBLE,
                 H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data()) < 0)
        throw std::runtime_error(std::string("HDF5: H5Dwrite failed for ") + name);
}

}  // namespace h5_lib


// ═══════════════════════════════════════════════════════════════════════════
// PUBLIC API
// ═══════════════════════════════════════════════════════════════════════════

/**
 * write_h5_indiviual
 *
 * For every placements[i], computes the individual fitness ratio via
 * get_individual_fitness() and writes all data to `h5_path` under a
 * freshly generated, memorable run group.
 *
 * The file is created if it does not exist, or appended to if it does — so
 * multiple runs from one session can accumulate in the same file.
 *
 * @param g             Global DEM / config context
 * @param grid          Feasible cells — the denominator for "ideal" fitness
 * @param h5_path       Path to output HDF5 file
 * @param placements    Radar placements as GridPoint{r, c}, length N
 * @param antennaHeight Antenna height AGL in metres   (default 50)
 * @param targetH       Target height in metres         (default 1 000)
 * @param minAngle      Min elevation angle in degrees  (default 0.5)
 * @param maxAngle      Max elevation angle in degrees  (default 30)
 * @param maxRangeM     Max radar range in metres       (default 25 000)
 *
 * @returns The run group name (e.g. "teal-ridge-143022") — use it to locate
 *          this run in the HDF5 file later.
 *
 * @throws  std::runtime_error on any HDF5 failure.
 */
inline std::string write_h5_indiviual(
    Global&                        g,
    const FeasibleGrid&            grid,
    const std::string&             h5_path,
    const std::vector<GridPoint>&  placements,
    double antennaHeight = 50.0,
    double targetH       = 1000.0,
    double minAngle      = 0.5,
    double maxAngle      = 30.0,
    double maxRangeM     = 25000.0)
{
    using namespace h5_lib;

    const size_t N = placements.size();

    // ── 1. Compute fitness for every placement ───────────────────────────
    // Each get_individual_fitness call is independent; parallelise across
    // placements. Guard against nested parallelism (e.g. caller already in
    // an OMP region) using the same pattern as radar_overlap.h.
    std::vector<double> fitness(N);

#if defined(_OPENMP)
    #pragma omp parallel for schedule(dynamic, 4) if(!omp_in_parallel())
#endif
    for (int i = 0; i < static_cast<int>(N); ++i)
        fitness[static_cast<size_t>(i)] = get_individual_coverage(
            g, grid,
            placements[i].r, placements[i].c,
            antennaHeight, targetH,
            minAngle, maxAngle, maxRangeM);

    // ── Unpack r / c into contiguous arrays for HDF5 ─────────────────────
    std::vector<int> rs(N), cs(N);
    for (size_t i = 0; i < N; ++i) { rs[i] = placements[i].r; cs[i] = placements[i].c; }

    // ── 2. Open existing file or create a new one ─────────────────────────
    // H5Fopen prints a diagnostic to stderr when the file does not exist even
    // if we handle the < 0 return — suppress HDF5 error output for this probe
    // only, then restore the original handler.
    H5E_auto2_t old_func;
    void*       old_data;
    H5Eget_auto2(H5E_DEFAULT, &old_func, &old_data);
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);          // silence
    hid_t raw_fid = H5Fopen(h5_path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    H5Eset_auto2(H5E_DEFAULT, old_func, old_data);        // restore

    if (raw_fid < 0)
        raw_fid = H5Fcreate(h5_path.c_str(),
                             H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (raw_fid < 0)
        throw std::runtime_error(
            "write_h5_indiviual: cannot open or create file '" + h5_path + "'");

    H5FileGuard fid{ raw_fid };

    // ── 3. Build run name: "run-<N>-<range>" ─────────────────────────────
    // N = number of groups already in the file root + 1, so the first run
    // in a fresh file is "run-1-...", the next is "run-2-...", etc.
    H5G_info_t root_info{};
    H5Gget_info(fid.id, &root_info);
    const int run_n = static_cast<int>(root_info.nlinks) + 1;

    // Express range as whole km when >= 1 km, otherwise as metres
    std::string range_str;
    if (maxRangeM >= 1000.0)
        range_str = std::to_string(
                        static_cast<int>(std::round(maxRangeM / 1000.0))) + "km";
    else
        range_str = std::to_string(
                        static_cast<int>(std::round(maxRangeM))) + "m";

    const std::string run_name = "run-" + std::to_string(run_n) + "-" + range_str;

    H5GroupGuard run{
        H5Gcreate2(fid.id, run_name.c_str(),
                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) };
    if (run.id < 0)
        throw std::runtime_error(
            "write_h5_indiviual: failed to create run group '" + run_name + "'");

    // ── 4. params/ — one scalar f64 dataset per radar parameter ──────────
    H5GroupGuard params{
        H5Gcreate2(run.id, "params",
                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) };
    if (params.id < 0)
        throw std::runtime_error("write_h5_indiviual: failed to create 'params' group");

    write_scalar(params.id, "antenna_height", antennaHeight);
    write_scalar(params.id, "target_height",  targetH);
    write_scalar(params.id, "min_angle_deg",  minAngle);
    write_scalar(params.id, "max_angle_deg",  maxAngle);
    write_scalar(params.id, "max_range_m",    maxRangeM);

    // ── 5. data/ — x, y (i32) and fitness (f64) arrays ─────────────────
    H5GroupGuard data{
        H5Gcreate2(run.id, "data",
                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) };
    if (data.id < 0)
        throw std::runtime_error("write_h5_indiviual: failed to create 'data' group");

    write_i32(data.id, "r",        rs);
    write_i32(data.id, "c",        cs);
    write_f64(data.id, "fitness", fitness);

    // RAII guards flush and close everything in reverse declaration order
    return run_name;
}

namespace h5_lib {

// One "run-<n>-<range>" group: the radar parameters it was scored with and
// the placements + individual coverage ratios.
struct RunData {
    std::string         name;
    double              antenna_height = 0.0;
    double              target_height  = 0.0;
    double              min_angle_deg  = 0.0;
    double              max_angle_deg  = 0.0;
    double              max_range_m    = 0.0;
    std::vector<int>    r, c;
    std::vector<double> fitness;
};

// A candidate radar site: placement + its individual coverage ratio.
struct SeedSite { int r; int c; double fitness; };

// Silences HDF5's automatic error printing for the lifetime of the object
// (probing for groups/datasets that may legitimately not exist).
struct H5ErrSilencer {
    H5E_auto2_t f{};
    void*       d{};
    H5ErrSilencer()  { H5Eget_auto2(H5E_DEFAULT, &f, &d); H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr); }
    ~H5ErrSilencer() { H5Eset_auto2(H5E_DEFAULT, f, d); }
};

inline bool read_scalar(hid_t group, const char* name, double& out)
{
    H5DSetGuard did{ H5Dopen2(group, name, H5P_DEFAULT) };
    if (did.id < 0) return false;
    return H5Dread(did.id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &out) >= 0;
}

template <class T>
inline bool read_1d(hid_t group, const char* name, hid_t memType, std::vector<T>& out)
{
    H5DSetGuard did{ H5Dopen2(group, name, H5P_DEFAULT) };
    if (did.id < 0) return false;
    H5SpaceGuard sid{ H5Dget_space(did.id) };
    if (sid.id < 0) return false;
    const hssize_t n = H5Sget_simple_extent_npoints(sid.id);
    if (n < 0) return false;
    out.resize(static_cast<size_t>(n));
    if (n == 0) return true;
    return H5Dread(did.id, memType, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data()) >= 0;
}

inline bool close_enough(double a, double b)
{
    return std::fabs(a - b) <= 1e-6 * std::max(1.0, std::fabs(b));
}

}  // namespace h5_lib


/**
 * read_h5_runs
 *
 * Reads every run group written by write_h5_indiviual. Groups that don't
 * have both params/ and data/ (for example "aoi" from write_heatmap_bin_h5)
 * and groups that fail to read are skipped. Returns an empty vector if the
 * file does not exist.
 */
inline std::vector<h5_lib::RunData> read_h5_runs(const std::string& h5_path)
{
    using namespace h5_lib;
    std::vector<RunData> runs;

    H5ErrSilencer quiet;
    H5FileGuard fid{ H5Fopen(h5_path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT) };
    if (fid.id < 0) return runs;

    H5G_info_t info{};
    if (H5Gget_info(fid.id, &info) < 0) return runs;

    for (hsize_t i = 0; i < info.nlinks; ++i) {
        const ssize_t len = H5Lget_name_by_idx(fid.id, ".", H5_INDEX_NAME, H5_ITER_INC,
                                               i, nullptr, 0, H5P_DEFAULT);
        if (len <= 0) continue;
        std::string name(static_cast<size_t>(len), '\0');
        H5Lget_name_by_idx(fid.id, ".", H5_INDEX_NAME, H5_ITER_INC,
                           i, &name[0], static_cast<size_t>(len) + 1, H5P_DEFAULT);

        H5GroupGuard run{ H5Gopen2(fid.id, name.c_str(), H5P_DEFAULT) };
        if (run.id < 0) continue;
        H5GroupGuard params{ H5Gopen2(run.id, "params", H5P_DEFAULT) };
        H5GroupGuard data  { H5Gopen2(run.id, "data",   H5P_DEFAULT) };
        if (params.id < 0 || data.id < 0) continue;

        RunData rd;
        rd.name = name;
        bool ok = read_scalar(params.id, "antenna_height", rd.antenna_height)
               && read_scalar(params.id, "target_height",  rd.target_height)
               && read_scalar(params.id, "min_angle_deg",  rd.min_angle_deg)
               && read_scalar(params.id, "max_angle_deg",  rd.max_angle_deg)
               && read_scalar(params.id, "max_range_m",    rd.max_range_m)
               && read_1d(data.id, "r",       H5T_NATIVE_INT,    rd.r)
               && read_1d(data.id, "c",       H5T_NATIVE_INT,    rd.c)
               && read_1d(data.id, "fitness", H5T_NATIVE_DOUBLE, rd.fitness);
        if (!ok || rd.r.size() != rd.c.size() || rd.r.size() != rd.fitness.size()) continue;
        runs.push_back(std::move(rd));
    }
    return runs;
}

/**
 * collect_seed_sites
 *
 * Merges all runs into one list of unique sites sorted by fitness
 * (best first). If matchParams is true, only runs scored with the same
 * radar parameters as the current homogenous radar are used, because the
 * coverage ratios of different radar models are not comparable.
 * A site that appears several times (the GA re-scores the same placement
 * every generation) is kept once, with its best value.
 */
inline std::vector<h5_lib::SeedSite> collect_seed_sites(
    const std::vector<h5_lib::RunData>& runs, bool matchParams,
    double antennaHeight, double targetH, double minAngle, double maxAngle, double maxRangeM)
{
    using namespace h5_lib;
    std::unordered_map<int64_t, SeedSite> best;
    for (const auto& run : runs) {
        // if (matchParams &&
            // !(close_enough(run.antenna_height, antennaHeight) &&
            //   close_enough(run.target_height,  targetH)       &&
            //   close_enough(run.min_angle_deg,  minAngle)      &&
            //   close_enough(run.max_angle_deg,  maxAngle)      &&
            //   close_enough(run.max_range_m,    maxRangeM)))
            // continue;
        for (size_t i = 0; i < run.r.size(); ++i) {
            if (run.r[i] < 0 || run.c[i] < 0) continue;
            const int64_t key = static_cast<int64_t>(run.r[i]) * 1'000'000LL + run.c[i];
            auto it = best.find(key);
            if (it == best.end() || run.fitness[i] > it->second.fitness)
                best[key] = SeedSite{ run.r[i], run.c[i], run.fitness[i] };
        }
    }
    std::vector<SeedSite> sites;
    std::cout<< "best: " << best.size() <<std::endl;
    sites.reserve(best.size());
    for (auto& kv : best) sites.push_back(kv.second);
    std::sort(sites.begin(), sites.end(), [](const SeedSite& a, const SeedSite& b) {
        if (a.fitness != b.fitness) return a.fitness > b.fitness;
        if (a.r != b.r) return a.r < b.r;
        return a.c < b.c;
    });
    std::cout<< "sites: " << sites.size() <<std::endl;
    return sites;
}

/**
 * deal_seed_chromosomes
 *
 * Builds up to K chromosomes of L genes from the best K*L sites by dealing
 * them round-robin: chromosome i gets sites i, i+K, i+2K, ... So every
 * chromosome mixes top-ranked and lower-ranked sites, and chromosome 0 gets
 * the best site overall. Sites must be unique and sorted best-first. If there
 * are fewer than K*L sites, fewer chromosomes are returned.
 */
inline std::vector<std::vector<h5_lib::SeedSite>> deal_seed_chromosomes(
    const std::vector<h5_lib::SeedSite>& sites, size_t K, size_t L, double topk, std::mt19937 rng)
{
    std::vector<std::vector<h5_lib::SeedSite>> out;
    int idx = std::uniform_int_distribution<int>(0, (int)(sites.size() *topk))(rng);

    if (K == 0 || L == 0) return out;
    K = std::min(K, sites.size() / L);
    if (K == 0) return out;
    out.assign(K, {});
    for (size_t j = 0; j < L; ++j)
        for (size_t i = 0; i < K; ++i)
            out[i].push_back(sites[j * K + i]);
    return out;
}
