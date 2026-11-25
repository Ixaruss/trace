// ============================================================================
// PSO single-radar siting on the project's own tiled bins.
//
// Finds the single best radar location inside the AOI by maximising UNIQUE
// coverage, scored by the project's own calculate_radar_coverage_overlap()
// (result[1] = AOI cells seen by exactly one radar; with one radar that is
// simply "cells this radar sees"). No feasible.bin / generate_bin.cc: the
// search space is Bin::getFeasibleGrid() of the loaded AOI window.
//
// Nothing in the existing project files is modified; this file only calls them.
// Everything lives in namespace pso_radar so names such as GridPoint /
// BoundingBox / ThreadPool cannot clash with project types (ODR).
// ============================================================================
#pragma once
// single-radar-deployment/pso.h - PSO single-radar siting (header only, no main()).
// Entry point: pso_radar::runAll(bool benchmark = false)
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include <numbers>
#include <random>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>
#include <json/json.h>

#include "../common.h"
#ifndef PSO_RADAR_H
#define PSO_RADAR_H "../multi-radar-deployment/radar.h"
#endif
#include PSO_RADAR_H

namespace pso_radar {

// ----------------------------------------------------------------------------
// Config
// ----------------------------------------------------------------------------
struct RadarSpec {
    double antenna = 20, minAngle = 0.5, maxAngle = 60, rangeKm = 150, target = 2000;
    int    peak_power = 660, antenna_gain = 30, crossSection = 1;
    double wavelength = 5, detachable_threshold = 3;
};
struct Hyper { double c1 = 1.494, c2 = 1.494, w = 0.729, velFrac = 0.05; };
struct Settings {
    RadarSpec radar;
    Hyper hyper;
    int swarmSize = 40, iterations = 25;
    unsigned seed = 0;
    bool use_curv = true, use_radar_model = false, useAMSL = false;
    bool write_png = true;
    int reduction = 50;
    std::string output = "output/pso/pso.tif";
    std::string topology = "global", inertia = "linear_decay", boundary = "bounce";
};

inline bool loadSettings(const std::string& path, Settings& s) {
    std::ifstream in(path);
    if (!in) return false;
    Json::Value root; Json::CharReaderBuilder b; std::string errs;
    if (!Json::parseFromStream(b, in, &root, &errs)) {
        std::cerr << "[pso] config parse error: " << errs << "\n"; return false;
    }
    auto d = [](const Json::Value& v, const char* k, double def) { return v.isMember(k) ? v[k].asDouble() : def; };
    auto i = [](const Json::Value& v, const char* k, int def)    { return v.isMember(k) ? v[k].asInt() : def; };
    auto bo= [](const Json::Value& v, const char* k, bool def)   { return v.isMember(k) ? v[k].asBool() : def; };
    auto st= [](const Json::Value& v, const char* k, const std::string& def) { return v.isMember(k) ? v[k].asString() : def; };

    const Json::Value r = root.get("radar", Json::Value(Json::objectValue));
    RadarSpec& R = s.radar;
    R.antenna = d(r, "antenna", R.antenna);       R.minAngle = d(r, "min-angle", R.minAngle);
    R.maxAngle = d(r, "max-angle", R.maxAngle);   R.rangeKm = d(r, "range", R.rangeKm);
    R.target = d(r, "target", R.target);          R.peak_power = i(r, "peak_power", R.peak_power);
    R.antenna_gain = i(r, "antenna_gain", R.antenna_gain);
    R.wavelength = d(r, "wavelength", R.wavelength);
    R.crossSection = i(r, "crossSection", R.crossSection);
    R.detachable_threshold = d(r, "detachable_threshold", R.detachable_threshold);

    s.use_curv = bo(root, "use_earth_curvature", s.use_curv);
    s.use_radar_model = bo(root, "use_radar_model", s.use_radar_model);
    s.useAMSL = bo(root, "useAMSL", s.useAMSL);
    s.output = st(root, "output", s.output);
    s.write_png = bo(root, "write_png", s.write_png);
    s.reduction = i(root, "reduction", s.reduction);

    const Json::Value h = root.get("pso_hyperparameters", Json::Value(Json::objectValue));
    s.hyper.c1 = d(h, "cognitive_weight_c1", s.hyper.c1);
    s.hyper.c2 = d(h, "social_weight_c2", s.hyper.c2);
    s.hyper.w  = d(h, "inertia_weight_w", s.hyper.w);
    s.hyper.velFrac = d(h, "velocity_clamp_fraction", s.hyper.velFrac);
    const Json::Value g = root.get("pso_global_settings", Json::Value(Json::objectValue));
    s.swarmSize = i(g, "swarm_size", s.swarmSize);
    s.iterations = i(g, "iterations", s.iterations);
    s.seed = (unsigned)i(g, "seed", 0);
    const Json::Value p = root.get("pso_strategy", Json::Value(Json::objectValue));
    s.topology = st(p, "topology", s.topology);
    s.inertia = st(p, "inertia", s.inertia);
    s.boundary = st(p, "boundary", s.boundary);
    return true;
}

// Owns the settings. Loads in the constructor from a fixed address and can be
// re-read at any time with reload() (a failed reload keeps the last good values).
class PsoConf {
public:
    static constexpr const char* kPath = "./config/psoconf.json";
    PsoConf() { reload(); }
    bool reload() {
        Settings fresh;
        if (!loadSettings(kPath, fresh)) {
            std::cerr << "[pso] cannot read " << kPath
                      << (loaded_ ? " - keeping previous settings\n" : "\n");
            return false;
        }
        s_ = fresh; loaded_ = true;
        std::cout << "[pso] config loaded from " << kPath << "\n";
        return true;
    }
    bool loaded() const { return loaded_; }
    const Settings& get() const { return s_; }
private:
    Settings s_;
    bool loaded_ = false;
};

// ----------------------------------------------------------------------------
// Geometry (grid coordinates: x = column, y = row)
// ----------------------------------------------------------------------------
struct BoundingBox { double minX, maxX, minY, maxY; };
struct Pt { double x = 0, y = 0; };
struct Particle {
    Pt pos; double vx = 0, vy = 0;
    Pt pBest; double pBestFit = -1, fit = 0;
};
struct Result {
    std::string name, topology, inertia, boundary;
    double coverage = 0;   // fraction of AOI
    Pt best; double ms = 0; long evals = 0;
};

// ----------------------------------------------------------------------------
// Strategies (kept from the original PSO)
// ----------------------------------------------------------------------------
struct IInertia { virtual ~IInertia() = default;
    virtual double weight(int it, int maxIt, double def) = 0; virtual std::string name() const = 0; };
struct IBoundary { virtual ~IBoundary() = default;
    virtual void handle(double& x, double& y, double& vx, double& vy, const BoundingBox& b) = 0;
    virtual std::string name() const = 0; };
struct ITopology { virtual ~ITopology() = default;
    virtual Pt target(int idx, int n, const Pt& gBest, const std::vector<Particle>& sw) = 0;
    virtual std::string name() const = 0; };

struct ConstantInertia : IInertia {
    double weight(int, int, double def) override { return def; }
    std::string name() const override { return "Constant"; } };
struct LinearDecayInertia : IInertia {
    double weight(int it, int maxIt, double) override { return 0.9 - 0.5 * double(it) / std::max(1, maxIt); }
    std::string name() const override { return "Linear Decay"; } };
struct RandomJitterInertia : IInertia {
    std::mt19937 rng{1337};
    double weight(int, int, double) override { return std::uniform_real_distribution<double>(0.4, 0.9)(rng); }
    std::string name() const override { return "Rand Jitter"; } };

struct ClampBoundary : IBoundary {
    void handle(double& x, double& y, double&, double&, const BoundingBox& b) override {
        x = std::clamp(x, b.minX, b.maxX); y = std::clamp(y, b.minY, b.maxY); }
    std::string name() const override { return "Clamp"; } };
struct BounceBoundary : IBoundary {
    void handle(double& x, double& y, double& vx, double& vy, const BoundingBox& b) override {
        if (x < b.minX) { x = b.minX + (b.minX - x); vx = -vx; }
        if (x > b.maxX) { x = b.maxX - (x - b.maxX); vx = -vx; }
        if (y < b.minY) { y = b.minY + (b.minY - y); vy = -vy; }
        if (y > b.maxY) { y = b.maxY - (y - b.maxY); vy = -vy; }
        x = std::clamp(x, b.minX, b.maxX); y = std::clamp(y, b.minY, b.maxY); }
    std::string name() const override { return "Bounce"; } };

struct GlobalTopology : ITopology {
    Pt target(int, int, const Pt& g, const std::vector<Particle>&) override { return g; }
    std::string name() const override { return "Global (Star)"; } };
struct RingTopology : ITopology {
    Pt target(int i, int n, const Pt&, const std::vector<Particle>& sw) override {
        const int l = (i - 1 + n) % n, r = (i + 1) % n;
        double best = sw[i].pBestFit; Pt t = sw[i].pBest;
        if (sw[l].pBestFit > best) { best = sw[l].pBestFit; t = sw[l].pBest; }
        if (sw[r].pBestFit > best) { t = sw[r].pBest; }
        return t; }
    std::string name() const override { return "Local (Ring)"; } };

inline std::shared_ptr<ITopology> makeTopology(const std::string& s) {
    if (s == "ring" || s == "local") return std::make_shared<RingTopology>();
    return std::make_shared<GlobalTopology>(); }
inline std::shared_ptr<IInertia> makeInertia(const std::string& s) {
    if (s == "constant") return std::make_shared<ConstantInertia>();
    if (s == "jitter" || s == "random_jitter") return std::make_shared<RandomJitterInertia>();
    return std::make_shared<LinearDecayInertia>(); }
inline std::shared_ptr<IBoundary> makeBoundary(const std::string& s) {
    if (s == "clamp") return std::make_shared<ClampBoundary>();
    return std::make_shared<BounceBoundary>(); }

// ----------------------------------------------------------------------------
// Engine: PSO over the project's feasible grid, fitness from the project's
// calculate_radar_coverage_overlap().
// ----------------------------------------------------------------------------
class Engine {
    Global& g;
    const Settings& S;
    FeasibleGrid f, aoi;
    std::vector<uint8_t> feasible;         // AABB-sized bitmap, O(1) membership
    BoundingBox bbox{};
    double total = 1;                      // denominator = AOI cell count
    std::mt19937 rng;

    std::unordered_map<uint64_t, double> cache;
    std::mutex cacheMx;
    std::atomic<long> evals{0};

    inline bool isFeasible(int r, int c) const {
        const int lr = r - f.origin_r, lc = c - f.origin_c;
        if (lr < 0 || lc < 0 || lr >= f.height || lc >= f.width) return false;
        return feasible[(size_t)lr * f.width + lc] != 0;
    }

    // Nearest feasible cell to (x=col, y=row): ring search, random feasible fallback.
    Pt snap(double x, double y, std::mt19937& rg) const {
        const int c0 = (int)std::lround(x), r0 = (int)std::lround(y);
        if (isFeasible(r0, c0)) return { (double)c0, (double)r0 };
        for (int rad = 1; rad <= 64; ++rad) {
            for (int dr = -rad; dr <= rad; ++dr) {
                const bool edge = (dr == -rad || dr == rad);
                const int step = edge ? 1 : 2 * rad;
                for (int dc = -rad; dc <= rad; dc += step)
                    if (isFeasible(r0 + dr, c0 + dc)) return { (double)(c0 + dc), (double)(r0 + dr) };
            }
        }
        const auto& cell = f.cells[std::uniform_int_distribution<size_t>(0, f.cells.size() - 1)(rg)];
        return { (double)cell.c, (double)cell.r };
    }

    double fitnessAt(const Pt& p) {
        const int r = (int)p.y, c = (int)p.x;
        const uint64_t key = ((uint64_t)(uint32_t)r << 32) | (uint32_t)c;
        { std::lock_guard<std::mutex> lk(cacheMx);
          auto it = cache.find(key); if (it != cache.end()) return it->second; }
        const RadarSpec& R = S.radar;
        std::vector<GridPoint> placement{ GridPoint{ r, c, 0 } };
        const auto res = calculate_radar_coverage_overlap(
            g, f, placement, 5.0, "", OutputCrs::Wgs84, BandLayout::Rgba,
            R.antenna, R.minAngle, R.target, R.maxAngle, R.rangeKm * 1000.0,
            {}, aoi, false, S.reduction,
            S.use_curv, false, 0.5, 60.0,
            S.use_radar_model, R.peak_power, R.antenna_gain, R.wavelength,
            R.crossSection, R.detachable_threshold, S.useAMSL);
        const double fit = (double)res[1] / total;     // unique coverage
        ++evals;
        std::lock_guard<std::mutex> lk(cacheMx);
        cache[key] = fit;
        return fit;
    }

    void evaluate(std::vector<Particle>& sw) {
        const unsigned nT = std::max(1u, std::min<unsigned>(std::thread::hardware_concurrency(), (unsigned)sw.size()));
        std::atomic<size_t> next{0};
        std::vector<std::thread> th;
        for (unsigned t = 0; t < nT; ++t)
            th.emplace_back([&] {
                rc_detail::tls_serial_rays = true;     // same as the GA workers
                for (size_t i; (i = next++) < sw.size();) sw[i].fit = fitnessAt(sw[i].pos);
            });
        for (auto& t : th) t.join();
    }

public:
    Engine(Global& gg, const Settings& st) : g(gg), S(st) {}

    // Load AOI window + build the search space from the project's bins.
    bool init() {
        // Radars are PLACED inside the AOI but SEE outward from it, so the
        // window must cover AOI + range (same reasoning as the GA).
        const RadarSpec& R = S.radar;
        double marginKm = R.rangeKm;
        if (S.use_radar_model) {
            const int modelKm = (int)std::pow(
                (R.peak_power * std::pow(R.antenna_gain, 2) * std::pow(R.wavelength, 2) * R.crossSection) /
                (std::pow(4 * std::numbers::pi, 2) * std::pow(10, R.detachable_threshold)), 0.25);
            marginKm = std::max(marginKm, (double)modelKm);
        }
        std::cout << "[pso] Loading AOI window (margin " << marginKm << " km)...\n";
        g.load({ CONSTRAINTS, ELEVATION }, AOI, 0.0, 0.0, marginKm);

        f   = Bin::getFeasibleGrid(g.bin, 0);
        aoi = Bin::getFeasibleGrid(g.bin, 0, 0, WITHALL);   // full AOI = scoring denominator
        if (f.cells.empty()) {
            std::cerr << "[pso] feasible grid is EMPTY (constraints loaded: "
                      << (g.bin.vector_loaded ? "yes" : "NO") << ", elevation loaded: "
                      << (g.bin.elevation_loaded ? "yes" : "NO") << ", window valid: "
                      << (g.bin.window_valid ? "yes" : "NO") << ")\n";
            return false;
        }
        total = aoi.cells.empty() ? (double)f.cells.size() : (double)aoi.cells.size();

        feasible.assign((size_t)f.height * f.width, 0);
        for (const auto& c : f.cells)
            feasible[(size_t)(c.r - f.origin_r) * f.width + (c.c - f.origin_c)] = 1;
        bbox = { (double)f.origin_c, (double)(f.origin_c + f.width - 1),
                 (double)f.origin_r, (double)(f.origin_r + f.height - 1) };
        std::cout << "[pso] feasible cells: " << f.cells.size() << "  AOI cells: " << total
                  << "  bbox cols [" << bbox.minX << "," << bbox.maxX << "] rows ["
                  << bbox.minY << "," << bbox.maxY << "]\n";
        return true;
    }

    Result run(const std::string& name, int swarmSize, int iterations,
               std::shared_ptr<ITopology> topo, std::shared_ptr<IInertia> inertia,
               std::shared_ptr<IBoundary> boundary) {
        const auto t0 = std::chrono::high_resolution_clock::now();
        evals = 0;
        rng.seed(S.seed ? S.seed : (unsigned)std::chrono::system_clock::now().time_since_epoch().count());

        const double maxVx = (bbox.maxX - bbox.minX) * S.hyper.velFrac;
        const double maxVy = (bbox.maxY - bbox.minY) * S.hyper.velFrac;
        std::uniform_real_distribution<double> vdx(-maxVx, maxVx), vdy(-maxVy, maxVy), u01(0, 1);
        std::uniform_int_distribution<size_t> pick(0, f.cells.size() - 1);

        std::vector<Particle> sw(swarmSize);
        for (auto& p : sw) {                      // start on random feasible cells
            const auto& c = f.cells[pick(rng)];
            p.pos = { (double)c.c, (double)c.r };
            p.vx = vdx(rng); p.vy = vdy(rng);
        }
        evaluate(sw);

        Pt gBest = sw[0].pos; double gFit = -1;
        for (auto& p : sw) {
            p.pBest = p.pos; p.pBestFit = p.fit;
            if (p.fit > gFit) { gFit = p.fit; gBest = p.pos; }
        }

        for (int it = 1; it <= iterations; ++it) {
            const double w = inertia->weight(it, iterations, S.hyper.w);
            for (int i = 0; i < swarmSize; ++i) {
                Particle& p = sw[i];
                const Pt tgt = topo->target(i, swarmSize, gBest, sw);
                const double r1 = u01(rng), r2 = u01(rng);
                p.vx = w * p.vx + S.hyper.c1 * r1 * (p.pBest.x - p.pos.x) + S.hyper.c2 * r2 * (tgt.x - p.pos.x);
                p.vy = w * p.vy + S.hyper.c1 * r1 * (p.pBest.y - p.pos.y) + S.hyper.c2 * r2 * (tgt.y - p.pos.y);
                p.vx = std::clamp(p.vx, -maxVx, maxVx);
                p.vy = std::clamp(p.vy, -maxVy, maxVy);
                p.pos.x += p.vx; p.pos.y += p.vy;
                boundary->handle(p.pos.x, p.pos.y, p.vx, p.vy, bbox);
                p.pos = snap(p.pos.x, p.pos.y, rng);      // keep on a feasible cell
            }
            evaluate(sw);
            for (auto& p : sw) {
                if (p.fit > p.pBestFit) { p.pBestFit = p.fit; p.pBest = p.pos; }
                if (p.fit >= gFit)      { gFit = p.fit;       gBest = p.pos; }
            }
            std::cout << "[pso] " << name << " iter " << std::setw(3) << it << "/" << iterations
                      << "  best unique coverage " << std::fixed << std::setprecision(2)
                      << gFit * 100.0 << " %\n";
        }

        Result res;
        res.name = name; res.topology = topo->name(); res.inertia = inertia->name();
        res.boundary = boundary->name(); res.coverage = gFit < 0 ? 0 : gFit; res.best = gBest;
        res.evals = evals.load();
        res.ms = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - t0).count();
        return res;
    }

    // Render the winner through the SAME coverage function (png/tif output).
    void render(const Pt& best) {
        const RadarSpec& R = S.radar;
        if (!S.output.empty()) {
            std::error_code ec;
            const auto dir = std::filesystem::path(S.output).parent_path();
            if (!dir.empty()) std::filesystem::create_directories(dir, ec);
        }
        std::vector<GridPoint> placement{ GridPoint{ (int)best.y, (int)best.x, 0 } };
        const auto res = calculate_radar_coverage_overlap(
            g, f, placement, 5.0, S.output, OutputCrs::Wgs84, BandLayout::Rgba,
            R.antenna, R.minAngle, R.target, R.maxAngle, R.rangeKm * 1000.0,
            {}, aoi, S.write_png, S.reduction,
            S.use_curv, false, 0.5, 60.0,
            S.use_radar_model, R.peak_power, R.antenna_gain, R.wavelength,
            R.crossSection, R.detachable_threshold, S.useAMSL);
        std::cout << "[pso] final render -> " << S.output << "  unique coverage "
                  << std::fixed << std::setprecision(2) << 100.0 * res[1] / total << " %\n";
    }

    void describe(const Pt& best) {
        double lat = 0, lon = 0;
        g.bin.getLatLngFromGrid((int)best.y, (int)best.x, lat, lon);
        std::cout << "[pso] best cell (row " << (int)best.y << ", col " << (int)best.x
                  << ")  lat " << std::fixed << std::setprecision(6) << lat << "  lon " << lon << "\n";
    }
};

// benchmark=false: one run with the strategy from pso_strategy ("run single").
// benchmark=true : all 12 topology/inertia/boundary profiles ("run single benchmark").
inline int runAll(bool benchmark = false) {
    PsoConf conf;                       // loads ./config/psoconf.json
    if (!conf.loaded()) return 1;
    // conf.reload() can be called again (e.g. between runs) to pick up edits.
    const Settings& S = conf.get();
    Global g;
    Engine eng(g, S);
    if (!eng.init()) return 2;

    std::vector<Result> table;
    if (benchmark) {
        int id = 1;
        for (auto topo : { std::string("global"), std::string("ring") })
            for (auto in : { std::string("constant"), std::string("linear_decay"), std::string("jitter") })
                for (auto bd : { std::string("clamp"), std::string("bounce") })
                    table.push_back(eng.run("Profile_" + std::to_string(id++), S.swarmSize, S.iterations,
                                            makeTopology(topo), makeInertia(in), makeBoundary(bd)));
    } else {
        table.push_back(eng.run("single", S.swarmSize, S.iterations,
                                makeTopology(S.topology), makeInertia(S.inertia), makeBoundary(S.boundary)));
    }
    std::sort(table.begin(), table.end(), [](const Result& a, const Result& b) { return a.coverage > b.coverage; });

    std::cout << "\n| Rank | Topology | Inertia | Boundary | Unique coverage | Evals | Time |\n"
              << "|:----:|:---------|:--------|:---------|:---------------:|:-----:|:----:|\n";
    int rank = 1;
    for (const auto& r : table)
        std::cout << "| " << rank++ << " | " << r.topology << " | " << r.inertia << " | " << r.boundary
                  << " | " << std::fixed << std::setprecision(2) << r.coverage * 100.0 << "% | "
                  << r.evals << " | " << std::setprecision(0) << r.ms << " ms |\n";

    eng.describe(table[0].best);
    eng.render(table[0].best);
    return 0;
}

} // namespace pso_radar
