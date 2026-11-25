#pragma once
// cli/commands.h
//
// One struct + one setup_*_cmd() (defines the CLI11 options) + one
// run_*_cmd() (executes, once parsing succeeds) per subcommand.
//
// CHANGED from the original main.cc: every subcommand now owns its own
// option variables. Previously vs_lat/png/single/epsg7755/srcLat/etc were
// declared ONCE at main()'s top level and reused across unrelated
// subcommands (rlos and los literally shared the same srcLat/srcLon/
// tgtLat/tgtLon/antenna/target variables; run's --png flag was the same
// bool as viewshed's --png). That had no runtime effect — CLI11's
// require_subcommand(1) means exactly one subcommand's handler ever runs —
// but it meant renaming or retuning one subcommand's option silently
// touched another's default, with no compiler warning either way. Splitting
// them is a pure reorganization: same flag names, same descriptions, same
// defaults, same behavior.

#include "cli.h"   // CLI::App
#include "format.h"
#include <string>

// ── build ("build") ────────────────────────────────────────────────────
struct BuildArgs {
    bool images = false, timer = false;
    bool all = false, road = false, rail = false;
    bool waterarea = false, waterline = false, ib = false, slope = false;
};
CLI::App* setup_build_cmd(CLI::App& app, BuildArgs& a);
void      run_build_cmd(const BuildArgs& a);

// ── viewshed ────────────────────────────────────────────────────────────
struct ViewshedArgs {
    double lat = 33.9871, lon = 74.7742;
    float  range = 300.0f, antennaHeight = 30.0f;
    bool   curvature = false;
    double loa = 0, roa = 90;
    bool   hasSector = false;
    bool   png = false, svg = false, tiff = false;
    bool   singleBand = false, epsg7755 = false;
    std::string output = "output/viewshed/viewshed.svg";
};
CLI::App* setup_viewshed_cmd(CLI::App& app, ViewshedArgs& a);
// Needs the parsed CLI::App* to check whether --output was explicitly
// given (viewshed_sub->count("--output")), for the sector/360 default
// output filename — every other subcommand here only needs its args struct.
void      run_viewshed_cmd(CLI::App* sub, const ViewshedArgs& a);

// ── radar ───────────────────────────────────────────────────────────────
struct RadarArgs {
    bool singleBand = false, epsg7755 = false;
};
CLI::App* setup_radar_cmd(CLI::App& app, RadarArgs& a);
int       run_radar_cmd(const RadarArgs& a);   // returns exit code

// ── rlos ────────────────────────────────────────────────────────────────
struct RlosArgs {
    double srcLat = 34.017346, srcLon = 74.50587999;
    double tgtLat = 34.16910062, tgtLon = 74.71267172;
    double antennaHeight = 30.0, targetHeight = 1000.0;
    double minElev = 0.5, maxElev = 30.0, maxRangeKm = 300.0;
    bool   detailed = false;
};
CLI::App* setup_rlos_cmd(CLI::App& app, RlosArgs& a);
void      run_rlos_cmd(const RlosArgs& a);

// ── los ─────────────────────────────────────────────────────────────────
struct LosArgs {
    double srcLat = 34.017346, srcLon = 74.50587999;
    double tgtLat = 34.16910062, tgtLon = 74.71267172;
    double antennaHeight = 30.0, targetHeight = 1000.0;
};
CLI::App* setup_los_cmd(CLI::App& app, LosArgs& a);
void      run_los_cmd(const LosArgs& a);

// ── run (GA) ────────────────────────────────────────────────────────────
struct GaArgs {
    int  gen = 0;
    bool runBenchmark = false;
    int  reduction = 50;
    bool png = false;
    bool singleBand = false, epsg7755 = false;
};
CLI::App* setup_ga_cmd(CLI::App& app, GaArgs& a);
int       run_ga_cmd(const GaArgs& a);   // returns exit code

// ── run single (PSO, single radar) ──────────────────────────────────────
// Nested under `run`:
//   <exe> run single             one PSO run (strategy from pso_strategy)
//   <exe> run single benchmark   all 12 strategy profiles, ranked table
// Everything else (radar spec, flags, swarm, output) comes from
// ./config/psoconf.json. Plain `run` stays the GA.
struct PsoSingleArgs { CLI::App* single = nullptr; CLI::App* benchmark = nullptr; };
CLI::App* setup_pso_single_cmd(CLI::App& run_sub, PsoSingleArgs& a);   // returns the `single` sub
// Call BEFORE the plain GA `run` handler: a nested subcommand also makes `run` truthy.
// Returns true (and sets exitCode) if `single` / `single benchmark` was given.
bool      run_pso_single_cmd(const PsoSingleArgs& a, int& exitCode);
