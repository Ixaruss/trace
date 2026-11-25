// cli/commands.cc
//
// Implementations for every subcommand declared in commands.h. Split out
// of main.cc so each subcommand's setup (option definitions) and execution
// (what happens once parsed) sit together and can be read/tested in
// isolation, instead of interleaved with six others inside one function.

#include "commands.h"
#include "../common.h"
#include "../operations.h"
#include "../multi-radar-deployment/config.h"
#include "../multi-radar-deployment/ga.h"
// PSO single-radar siting (header only, everything in namespace pso_radar).
#include "../single-radar-deployment/pso.h"
#include <iostream>
#include <iomanip>
#include <memory>

using std::cout;
using std::endl;
using std::string;

// ═══════════════════════════════════════════════════════════════════════
// build
// ═══════════════════════════════════════════════════════════════════════
CLI::App* setup_build_cmd(CLI::App& app, BuildArgs& a) {
    auto* sub = app.add_subcommand("build", "Build bin files from corresponding maps. config REQUIRED");
    auto* op_group = sub->add_option_group("individual outputs");

    sub->add_flag("-v,--with-images", a.images, "Generates images for bin visualization");
    sub->add_flag("-c,--count-time",  a.timer,  "Counts the time taken for preprocessing");

    op_group->add_flag("-a,--all",          a.all,       "create all bins");
    op_group->add_flag("--road",            a.road,      "create road bin");
    op_group->add_flag("--rail",            a.rail,      "create rail bin");
    op_group->add_flag("--water-areas",     a.waterarea, "create water-areas bin");
    op_group->add_flag("--water-lines",     a.waterline, "create water-lines bin");
    op_group->add_flag("-i, --ib",          a.ib,        "create ib bin");
    op_group->add_flag("--slope",           a.slope,     "create slope bin");
    // FIX: was require_option(0, 1) — zero flags was legal, and zero flags
    // means none of road/rail/waterarea/waterline/ib/slope get set, which
    // leaves `i` at -1 in run_build_cmd -> Global::preprocess treats -1 as
    // "process every layer." So a full build happened whenever NO flag was
    // given, regardless of -a specifically — -a was never actually checked
    // for that decision (see run_build_cmd; -a only reaches the `seperate`
    // parameter). Requiring exactly one flag means omitting all of them is
    // now a CLI11 usage error instead of a silent full build.
    op_group->require_option(1, 1);

    return sub;
}

void run_build_cmd(const BuildArgs& a) {
    int i = -1;
    if      (a.road)      i = 1;
    else if (a.rail)      i = 2;
    else if (a.waterarea) i = 3;
    else if (a.waterline) i = 4;
    else if (a.ib)        i = 5;
    else if (a.slope)     i = 6;
    // a.all -> "seperate" param (per-layer .bin also written standalone),
    // matching Global::preprocess(bool images, bool timer, bool seperate,
    // int individual) exactly as the original called it.
    Global::preprocess(a.images, a.timer, a.all, i);
}

// ═══════════════════════════════════════════════════════════════════════
// viewshed
// ═══════════════════════════════════════════════════════════════════════
CLI::App* setup_viewshed_cmd(CLI::App& app, ViewshedArgs& a) {
    auto* sub = app.add_subcommand("viewshed", "Calculate viewshed for a given location");

    sub->add_option("--lat", a.lat, "Latitude of the origin")->required();
    sub->add_option("--lon", a.lon, "Longitude of the origin")->required();
    sub->add_option("-r,--range", a.range, "Range limit in kilometers");
    sub->add_option("-a,--antenna", a.antennaHeight, "Antenna height in meters");
    sub->add_option("-o,--output", a.output, "Output SVG filename");
    sub->add_flag("-s", a.hasSector, "Use sector");
    sub->add_option("--loa", a.loa, "Left of arc (degrees)");
    sub->add_option("--roa", a.roa, "Right of arc (degrees)");
    sub->add_flag("--single-band", a.singleBand, "create single band tiff instead of default multi-band");
    sub->add_flag("--epsg-7755", a.epsg7755, "reproject tiff to epsg 7755 instead of deafult epsg 4326");
    sub->add_flag("--use-curvature", a.curvature, "use earth curvature corrections");

    auto* vg = sub->add_option_group("output formats available");
    vg->add_flag("--tif", a.tiff, "geotagged tiff");
    vg->add_flag("--svg", a.svg,  "svg vector");
    vg->add_flag("--png", a.png,  "png image (default)");
    vg->require_option(1);

    return sub;
}

void run_viewshed_cmd(CLI::App* sub, const ViewshedArgs& a) {
    Global g;
    g.load({ELEVATION}, OperationType::UNBOUNDED, a.lat, a.lon, a.range);
    auto engine = std::make_unique<ViewshedEngine>(
        g, a.lat, a.lon, a.range, a.antennaHeight, a.curvature
    );

    const OutputFormat fmt = resolve_output_format(a.tiff, a.svg, a.png);
    const string ext = extension_for(fmt);

    string output = a.output;
    if (sub->count("--output") == 0) {
        output = a.hasSector
            ? "output/viewshed/viewshed_sector" + ext
            : "output/viewshed/viewshed_360"    + ext;
    }

    // was hardcoded 50 in both branches of the original; viewshed doesn't
    // expose a --reduce option the way run/sgrid do, so this stays a fixed
    // constant here rather than silently gaining a new flag.
    constexpr int kReduction = 50;
    const OutputCrs   crs  = a.epsg7755   ? OutputCrs::Native7755        : OutputCrs::Wgs84;
    const BandLayout   band = a.singleBand ? BandLayout::SingleBandPalette : BandLayout::Rgba;

    if (a.hasSector) engine->computeViewshed(a.loa, a.roa);
    else             engine->computeViewshed();

    // CHANGED: the original only printed "writing: ..." on the non-sector
    // path (the sector branch's write happened silently) — almost
    // certainly an oversight rather than intentional, since nothing else
    // about the two branches differs after computeViewshed(). Printed in
    // both now.
    cout << "writing: " << output << endl;
    write_viewshed_output(*engine, output, fmt, kReduction, crs, band);
}

// ═══════════════════════════════════════════════════════════════════════
// radar
// ═══════════════════════════════════════════════════════════════════════
CLI::App* setup_radar_cmd(CLI::App& app, RadarArgs& a) {
    auto* sub = app.add_subcommand("radar", "Calculate radar profile coverage.");
    // radar's output path/format/range/etc all come from Rconf (radar
    // config), not CLI flags — only these two rendering options are
    // exposed directly, matching the original.
    (void)a; // singleBand/epsg7755 read in run_radar_cmd via the same struct
    return sub;
}

int run_radar_cmd(const RadarArgs& a) {
    // FIX: this used to run Rconf::load() (and build Global g / RadarEngine
    // off whatever that returns) BEFORE checking Rconf::is_available() —
    // three lines later. When no radar config exists, that's an engine
    // constructed from a default/empty Rconf before the process ever
    // notices anything is wrong. Same "load first, check after" shape as
    // the applyConfig() bug fixed earlier in this codebase. Check first.
    if (!Rconf::is_available()) {
        std::cerr << "[ERROR] radar conf not found\n";
        return -1;
    }
    auto rconf = Rconf::load();

    // NOTE: preserved verbatim from the original, including its handling
    // of a 4-letter ".tiff" extension — substr(size-3,3) on ".tiff" yields
    // "iff" (the leading 't' is dropped), which is exactly why the tiff
    // check below tests for "tif" OR "iff" rather than just "tif". Not
    // touched here; changing it would change which files radar.json can
    // legally point at.
    string output_ext = rconf.output.substr(rconf.output.size() - 3, 3);
    string path        = rconf.output.substr(0, rconf.output.size() - 3);

    Global g;
    g.load({ELEVATION}, OperationType::UNBOUNDED, rconf.srcLat, rconf.srcLon, rconf.range);
    RadarEngine engine(g);

    OutputFormat fmt;
    if      (output_ext == "tif" || output_ext == "iff") fmt = OutputFormat::Tiff;
    else if (output_ext == "svg")                        fmt = OutputFormat::Svg;
    else                                                  fmt = OutputFormat::Png;

    string r_output = path + output_ext;
    auto res = engine.compute();
    cout << "writing: " << r_output << endl;

    write_radar_output(engine, res, r_output, fmt, 50,
                       a.epsg7755   ? OutputCrs::Native7755        : OutputCrs::Wgs84,
                       a.singleBand ? BandLayout::SingleBandPalette : BandLayout::Rgba);
    return 0;
}

// ═══════════════════════════════════════════════════════════════════════
// rlos
// ═══════════════════════════════════════════════════════════════════════
CLI::App* setup_rlos_cmd(CLI::App& app, RlosArgs& a) {
    auto* sub = app.add_subcommand("rlos", "Calculate radar Line of Sight (LOS) between a source and a target.");
    sub->add_option("--src-lat", a.srcLat, "Source (radar) latitude")->required();
    sub->add_option("--src-lon", a.srcLon, "Source (radar) longitude")->required();
    sub->add_option("--tgt-lat", a.tgtLat, "Target latitude")->required();
    sub->add_option("--tgt-lon", a.tgtLon, "Target longitude")->required();
    sub->add_option("-a,--antenna", a.antennaHeight, "Radar antenna height AGL (meters)");
    sub->add_option("-t,--target",  a.targetHeight,  "Target aircraft height AGL (meters)");
    sub->add_option("--min-elev", a.minElev, "Minimum elevation angle (degrees)");
    sub->add_option("--max-elev", a.maxElev, "Maximum elevation angle (degrees)");
    sub->add_option("-r,--range", a.maxRangeKm, "Maximum calculation range (kilometers)");
    sub->add_flag("-d,--detailed", a.detailed, "Execute detailed profile and print cell telemetry");
    return sub;
}

void run_rlos_cmd(const RlosArgs& a) {
    Global g;
    if (!a.detailed) {
        auto res = RadarEngine::radarLOS(
            g, a.srcLat, a.srcLon, a.antennaHeight,
            a.tgtLat, a.tgtLon, a.targetHeight,
            a.minElev, a.maxElev, a.maxRangeKm
        );
        RadarEngine::printLOSResult(res, a.srcLat, a.srcLon, a.tgtLat, a.tgtLon);
        return;
    }

    auto det = RadarEngine::radarLOSDetailed(
        g, a.srcLat, a.srcLon, a.antennaHeight,
        a.tgtLat, a.tgtLon, a.targetHeight,
        a.minElev, a.maxElev, a.maxRangeKm
    );

    cout << "\nDETAILED — " << det.cells.size() << " cells along ray\n";
    cout << "First blocking cell: ";
    if (det.summary.status == LOSStatus::TERRAIN_MASKED) {
        cout << "row=" << det.summary.blockRow
             << " col=" << det.summary.blockCol
             << " dist=" << det.summary.blockDistM / 1000.0 << "km\n";
    } else {
        cout << "none — " << losStatusStr(det.summary.status) << "\n";
    }

    cout << "\nRow    Col    Dist(km)  TerrElev  CorrElev  TerrSlope   TargSlope   Blocks\n";
    int n = std::min((int)det.cells.size(), 10);
    for (int i = 0; i < n; i++) {
        const auto& c = det.cells[i];
        printf("%-6d %-6d %-9.2f %-9.1f %-9.1f %-11.6f %-11.6f %s\n",
               c.row, c.col, c.distM / 1000.0,
               (double)c.terrainElev, (double)c.correctedElev,
               c.terrainSlope, c.targetSlope,
               c.blocks ? "BLOCK" : "ok");
    }
}

// ═══════════════════════════════════════════════════════════════════════
// los
// ═══════════════════════════════════════════════════════════════════════
CLI::App* setup_los_cmd(CLI::App& app, LosArgs& a) {
    auto* sub = app.add_subcommand("los", "Calculate Line of Sight (LOS) between a source and a target.");
    sub->add_option("--src-lat", a.srcLat, "Source (radar) latitude")->required();
    sub->add_option("--src-lon", a.srcLon, "Source (radar) longitude")->required();
    sub->add_option("--tgt-lat", a.tgtLat, "Target latitude")->required();
    sub->add_option("--tgt-lon", a.tgtLon, "Target longitude")->required();
    sub->add_option("-a,--antenna", a.antennaHeight, "Radar antenna height AGL (meters)");
    sub->add_option("-t,--target",  a.targetHeight,  "Target aircraft height AGL (meters)");
    return sub;
}

void run_los_cmd(const LosArgs& a) {
    Global g;
    // Window must be resolved before any Bin access; this subcommand has no
    // --range flag, so size the window from the actual source-target
    // great-circle distance (haversine -- no CRS/poCT exists yet before
    // load()) plus a safety margin.
    {
        const double dLat = Global::deg2rad(a.tgtLat - a.srcLat);
        const double dLon = Global::deg2rad(a.tgtLon - a.srcLon);
        const double hv = std::sin(dLat / 2) * std::sin(dLat / 2) +
                           std::cos(Global::deg2rad(a.srcLat)) * std::cos(Global::deg2rad(a.tgtLat)) *
                           std::sin(dLon / 2) * std::sin(dLon / 2);
        const double srcTgtKm = 2.0 * 6371.0 * std::atan2(std::sqrt(hv), std::sqrt(1.0 - hv));
        g.load({ELEVATION}, OperationType::UNBOUNDED, a.srcLat, a.srcLon, srcTgtKm + 5.0);
    }

    // lineOfVisibilityopt/lineOfVisibility/lineOfVisibilityOptimized were
    // duplicate LOS implementations never updated for windowed/tiled loading
    // and have been removed -- Global::lineOfSight (lineofsight.cc) is the
    // one real two-point LOS check left, already correct for tiled mode.
    // It prints its own src/target elevation + timing breakdown.
    bool visible = g.lineOfSight(a.srcLat, a.srcLon, (int)a.antennaHeight,
                                  a.tgtLat, a.tgtLon, (int)a.targetHeight);
    cout << "Target Status : " << (visible ? "VISIBLE" : "BLOCKED") << "\n";
}

// ═══════════════════════════════════════════════════════════════════════
// run (GA)
// ═══════════════════════════════════════════════════════════════════════
CLI::App* setup_ga_cmd(CLI::App& app, GaArgs& a) {
    auto* sub = app.add_subcommand("run", "run ga cases via config");
    sub->add_option("-g,--gen-simdata", a.gen, "write simulation for top k solutions for ga");
    sub->add_option("--reduce", a.reduction, "reduction percentage for output resolution quality");
    sub->add_flag("--png", a.png, "write output as png as well");
    sub->add_flag("-b, --run-benchmark", a.runBenchmark, "create benchamrk table for all selection, mutation and crossover methods");
    sub->add_flag("--single-band", a.singleBand, "create single band tiff instead of default multi-band");
    sub->add_flag("--epsg-7755", a.epsg7755, "reproject tiff to epsg 7755 instead of deafult epsg 4326");
    return sub;
}

int run_ga_cmd(const GaArgs& a) {
    ThreadPool pool;
    GaEngine ga(pool);

    if (!GaConf::is_available()) return 0;   // matches original: silently does nothing without a GaConf
    GaConf cfg = GaConf::load();

    const OutputCrs  crs  = a.epsg7755   ? OutputCrs::Native7755        : OutputCrs::Wgs84;
    const BandLayout band = a.singleBand ? BandLayout::SingleBandPalette : BandLayout::Rgba;

    if (a.runBenchmark) {
        if      (cfg.homo.is_case)   ga.run_benchmark(BenchmarkCase::HOMOGENOUS);
        else if (cfg.hetero.is_case) ga.run_benchmark(BenchmarkCase::HETEROGENOUS);
        else if (cfg.autoG.is_case)  ga.run_benchmark(BenchmarkCase::AUTO);
        return 0;
    }
    if (a.gen != 0) {
        ga.run_native_random_simulation(a.gen);
        return 0;
    }
    if (cfg.homo.is_case) {
        ga.run_homogenous(a.reduction, a.png, crs, band);
        return 0;
    }
    if (cfg.hetero.is_case) {
        ga.run_heterogenous(a.reduction, a.png, crs, band);
        return 0;
    }
    if (cfg.autoG.is_case) {
        ga.setAutoOperators(SelectionMethod::TOURNAMENT,
                            AutoCrossoverMethod::RANDOM_SUBSET,
                            /*mRate=*/0.10, /*cRate=*/0.85);
        ga.run_auto(a.reduction, a.png, crs, band);
        return 0;
    }

    std::cerr << "[ga] Config found but no is_case flag is set — exiting\n";
    return -1;
}

// ═══════════════════════════════════════════════════════════════════════
// run single (PSO)
// ═══════════════════════════════════════════════════════════════════════
CLI::App* setup_pso_single_cmd(CLI::App& run_sub, PsoSingleArgs& a) {
    a.single = run_sub.add_subcommand("single",
        "find the best location for a single radar with PSO (config: ./config/psoconf.json)");
    a.benchmark = a.single->add_subcommand("benchmark",
        "run every PSO strategy profile (topology x inertia x boundary) and print a ranking");
    return a.single;
}

bool run_pso_single_cmd(const PsoSingleArgs& a, int& exitCode) {
    if (a.benchmark && *a.benchmark) { exitCode = pso_radar::runAll(true);  return true; }
    if (a.single    && *a.single)    { exitCode = pso_radar::runAll(false); return true; }
    return false;
}
