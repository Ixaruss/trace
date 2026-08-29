#include "cli.h"
#include "../common.h"
#include "../operations.h"
#include "../ga/config.h"
#include "../ga/ga.h"
#include "../iad/coverage.h"
int main(int argc, char** argv) {

    CLI::App app{"Geospatial Analysis and Radar CLI Tool"};
    app.require_subcommand(1);

    auto pre_sub = app.add_subcommand("build", "Build bin files from corresponding maps. config REQUIRED");
    bool images = false, timer = false, s = false, road=false, rail = false, waterarea = false, waterline = false, ib = false, slope=false;
    auto* op_group = pre_sub->add_option_group("individual outputs");
    pre_sub->add_flag("-v,--with-images",images,"Generates images for bin visualization");
    pre_sub->add_flag("-c,--count-time",timer,"Counts the time taken for preprocessing");
    op_group->add_flag("-a,--all",s,"create all bins");
    op_group->add_flag("--road",road,"create road bin");
    op_group->add_flag("--rail",rail,"create rail bin");
    op_group->add_flag("--water-areas",waterarea,"create water-areas bin");
    op_group->add_flag("--water-lines",waterline,"create water-lines bin");
    op_group->add_flag("-i, --ib",ib,"create ib bin");
    op_group->add_flag("--slope",slope,"create slope bin");
    op_group->require_option(0,1);

    auto viewshed_sub = app.add_subcommand("viewshed", "Calculate viewshed for a given location");
    double vs_lat = 33.9871;
    double vs_lon = 74.7742;
    float vs_range = 300.0f;
    float vs_ant_height = 30.0f;
    bool vs_curvature = false;
    double loa = 0;
    double roa = 90;
    bool has_sector = false;
    bool png = false, svg = false, tiff = false;
    bool single = false;
    bool epsg7755 = false;
    bool use_amsl = false;
    std::string vs_output = "output/viewshed/viewshed.svg";

    viewshed_sub->add_option("--lat", vs_lat, "Latitude of the origin")->required();
    viewshed_sub->add_option("--lon", vs_lon, "Longitude of the origin")->required();
    viewshed_sub->add_option("-r,--range", vs_range, "Range limit in kilometers");
    viewshed_sub->add_option("-a,--antenna", vs_ant_height, "Antenna height in meters");
    viewshed_sub->add_option("-o,--output", vs_output, "Output SVG filename");
    viewshed_sub->add_flag("-s", has_sector, "Use sector");
    viewshed_sub->add_option("--loa", loa, "Left of arc (degrees)");
    viewshed_sub->add_option("--roa", roa, "Right of arc (degrees)");
    viewshed_sub->add_flag("--single-band", single, "create single band tiff instead of default multi-band");
    viewshed_sub->add_flag("--epsg-7755", epsg7755, "reproject tiff to epsg 7755 instead of deafult epsg 4326");

    viewshed_sub->add_flag("--use-curvature", vs_curvature, "use earth curvature corrections");

    auto vg = viewshed_sub->add_option_group("output formats available");
    vg->add_flag("--tif", tiff,"geotagged tiff");
    vg->add_flag("--svg",svg, "svg vector");
    vg->add_flag("--png", png, "png image (default)");
    vg->require_option(1);
    auto radar_sub = app.add_subcommand("radar", "Calculate radar profile coverage.");

    double r_lat = 33.9871;
    double r_lon = 74.7742;
    double r_ant_agl = 50.0;
    double r_tgt_agl = 1000.0;
    double r_min_elev = 10.0;
    double r_max_elev = 30.0;
    double r_max_range = 300.0;

    std::string r_output = "radar_output.svg";
    bool r_fast = false;


    radar_sub->add_option("--lat", r_lat, "Latitude of radar base")->required();
    radar_sub->add_option("--lon", r_lon, "Longitude of radar base")->required();
    radar_sub->add_option("-a,--antenna", r_ant_agl, "Antenna height AGL (meters)");
    radar_sub->add_option("-t,--target", r_tgt_agl, "Target/threat aircraft height AGL (meters)");
    radar_sub->add_option("--min-elev", r_min_elev, "Minimum elevation angle (degrees)");
    radar_sub->add_option("--max-elev", r_max_elev, "Maximum elevation angle (degrees)");
    radar_sub->add_option("-r,--range", r_max_range, "Maximum radar range (kilometers)");
    radar_sub->add_option("-o,--output", r_output, "Output SVG filename");
    radar_sub->add_flag("-s", has_sector, "Use sector");
    radar_sub->add_option("--loa", loa, "Left of arc (degrees)");
    radar_sub->add_option("--roa", roa, "Right of arc (degrees)");
    radar_sub->add_flag("-f,--fast", r_fast, "Execute calculation in fast mode (without image)");
    radar_sub->add_flag("--use-curvature", vs_curvature, "use earth curvature corrections");
    radar_sub->add_flag("--single-band", single, "create single band tiff instead of default multi-band");
    radar_sub->add_flag("--epsg-7755", epsg7755, "reproject tiff to epsg 7755 instead of deafult epsg 4326");
    radar_sub->add_flag("--use-amsl", use_amsl, "use AMSL height instead of AGL for target/threat height");

    auto rg = radar_sub->add_option_group("output formats available");
    rg->add_flag("--tif", tiff,"geotagged tiff");
    rg->add_flag("--svg",svg, "svg vector");
    rg->add_flag("--png", png, "png image (default)");
    rg->require_option(1);

    double srcLat = 34.017346;
    double srcLon = 74.50587999;
    double tgtLat = 34.16910062;
    double tgtLon = 74.71267172;
    double los_ant_height = 30.0;
    double los_tgt_height = 1000.0;
    double los_min_elev = 0.5;
    double los_max_elev = 30.0;
    double los_max_range = 300.0;
    bool detailed = false;
    auto rlos_sub = app.add_subcommand("rlos", "Calculate radar Line of Sight (LOS) between a source and a target.");
    rlos_sub->add_option("--src-lat", srcLat, "Source (radar) latitude")->required();
    rlos_sub->add_option("--src-lon", srcLon, "Source (radar) longitude")->required();
    rlos_sub->add_option("--tgt-lat", tgtLat, "Target latitude")->required();
    rlos_sub->add_option("--tgt-lon", tgtLon, "Target longitude")->required();
    rlos_sub->add_option("-a,--antenna", los_ant_height, "Radar antenna height AGL (meters)");
    rlos_sub->add_option("-t,--target", los_tgt_height, "Target aircraft height AGL (meters)");
    rlos_sub->add_option("--min-elev", los_min_elev, "Minimum elevation angle (degrees)");
    rlos_sub->add_option("--max-elev", los_max_elev, "Maximum elevation angle (degrees)");
    rlos_sub->add_option("-r,--range", los_max_range, "Maximum calculation range (kilometers)");

    rlos_sub->add_flag("-d,--detailed", detailed, "Execute detailed profile and print cell telemetry");

    auto los_sub = app.add_subcommand("los", "Calculate Line of Sight (LOS) between a source and a target.");
    los_sub->add_option("--src-lat", srcLat, "Source (radar) latitude")->required();
    los_sub->add_option("--src-lon", srcLon, "Source (radar) longitude")->required();
    los_sub->add_option("--tgt-lat", tgtLat, "Target latitude")->required();
    los_sub->add_option("--tgt-lon", tgtLon, "Target longitude")->required();
    los_sub->add_option("-a,--antenna", los_ant_height, "Radar antenna height AGL (meters)");
    los_sub->add_option("-t,--target", los_tgt_height, "Target aircraft height AGL (meters)");
    los_sub->add_flag("-d,--detailed", detailed, "Execute detailed profile and print cell telemetry");

    auto run  = app.add_subcommand("run", "run ga cases via config");
    int gen = 0;
    bool run_benchmark = false;
    int reduction_factor = 50;
     run->add_option("-g,--gen-simdata",gen, "write simulation for top k solutions for ga");
     run->add_option("--reduce",reduction_factor, "reduction percentage for output resolution quality");
     run-> add_flag("--png", png,"write output as png as well");
     run->add_flag("-b, --run-benchmark", run_benchmark, "create benchamrk table for all selection, mutation and crossover methods");
     run->add_flag("--single-band", single, "create single band tiff instead of default multi-band");
     run->add_flag("--epsg-7755", epsg7755, "reproject tiff to epsg 7755 instead of deafult epsg 4326");

     auto sgrid  = app.add_subcommand("sgrid", "create a survaillence grid cooresponding to iadconf");


     sgrid->add_flag("--single-band", single, "create single band tiff instead of default multi-band");
     sgrid->add_flag("--epsg-7755", epsg7755, "reproject tiff to epsg 7755 instead of deafult epsg 4326");

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        return app.exit(e);
    }


    if (*viewshed_sub) {
        Global g;
        auto engine = std::make_unique<ViewshedEngine>(
            g, vs_lat, vs_lon, vs_range, vs_ant_height, vs_curvature
        );
        string output_ext;
        if(tiff){
            output_ext = ".tif";
        }else if(svg) {
            output_ext = ".svg";
        } else {
            output_ext = ".png";
        }
        if (viewshed_sub->count("--output") == 0 && has_sector) {
            vs_output = "output/viewshed/viewshed_sector" + output_ext;
        } else if (viewshed_sub->count("--output") == 0) {
            vs_output = "output/viewshed/viewshed_360" + output_ext;
        }
         if (has_sector) {
             engine->computeViewshed(loa,roa);
             if(tiff){
                  engine->writeViewshedToGeoTIFF(vs_output,50,epsg7755 ? OutputCrs::Native7755 : OutputCrs::Wgs84, single ? BandLayout::SingleBandPalette : BandLayout::Rgba);
             }else if(svg) {
                  engine->writeViewshedToSVG(vs_output);
             } else {
                 engine->writeViewshedToPNG(vs_output);
             }
         } else {
             engine->computeViewshed();
             cout << "writing: " << vs_output <<endl;
             if(tiff){
                  engine->writeViewshedToGeoTIFF(vs_output,50,epsg7755 ? OutputCrs::Native7755 : OutputCrs::Wgs84, single ? BandLayout::SingleBandPalette : BandLayout::Rgba);
             }else if(svg) {
                  engine->writeViewshedToSVG(vs_output);
             } else {
                 engine->writeViewshedToPNG(vs_output);
             }
         }
    }

    // Route 2: Execution for Radar
    else if (*radar_sub) {

        Global g;
        g.bin.init({ ELEVATION });
        RadarEngine engine(
            g, r_lat, r_lon, r_ant_agl, r_tgt_agl, r_min_elev, r_max_elev, r_max_range, vs_curvature,use_amsl
        );
        string output_ext;
        if(tiff){
            output_ext = ".tif";
        }else if(svg) {
            output_ext = ".svg";
        } else {
            output_ext = ".png";
        }
        if (r_fast) {
            if (has_sector) {

                engine.computeFast(loa,roa);
            } else {
                engine.computeFast();
            }
        } else {
            if (radar_sub->count("--output") == 0 && has_sector) {
                r_output = "output/radar/radar_sector" +  output_ext;
            } else if (radar_sub->count("--output") == 0) {
                r_output = "output/radar/radar_360" + output_ext;
            }

            if (has_sector) {
                auto res = engine.compute(loa,roa);
                if(tiff){
                     engine.writeCoverageToGeoTIFF(res,r_output,50,epsg7755 ? OutputCrs::Native7755 : OutputCrs::Wgs84, single ? BandLayout::SingleBandPalette : BandLayout::Rgba);
                }else if(svg) {
                     engine.writeToSVG(res,r_output);
                } else {
                    engine.writeToPNG(res,r_output);
                }
            } else {
                auto res = engine.compute();
                cout << "writing: " << r_output <<endl;
                if(tiff){
                     engine.writeCoverageToGeoTIFF(res,r_output,50,epsg7755 ? OutputCrs::Native7755 : OutputCrs::Wgs84, single ? BandLayout::SingleBandPalette : BandLayout::Rgba);
                }else if(svg) {
                     engine.writeToSVG(res,r_output);
                } else {
                    engine.writeToPNG(res,r_output);
                }
            }
        }
    }
    else if (*pre_sub) {
        int i = -1;
        if (road) i = 1;
        else if (rail) i = 2;
        else if(waterarea) i =3;
        else if (waterline) i =4;
        else if (ib) i =5;
        else if (slope) i =6;
        Global::preprocess(images,timer,s,i);
    }
    else if (*rlos_sub) {
        Global g;
        if (!detailed) {
            // 1. Fast LOS Workflow
            auto res = RadarEngine::radarLOS(
                g,
                srcLat, srcLon, los_ant_height,
                tgtLat, tgtLon, los_tgt_height,
                los_min_elev, los_max_elev, los_max_range
            );

            RadarEngine::printLOSResult(res, srcLat, srcLon, tgtLat, tgtLon);
        }
        else {
            // 2. Detailed LOS Workflow (for visualization / debugging)
            auto det = RadarEngine::radarLOSDetailed(
                g,
                srcLat, srcLon, los_ant_height,
                tgtLat, tgtLon, los_tgt_height,
                los_min_elev, los_max_elev, los_max_range
            );

            std::cout << "\nDETAILED — " << det.cells.size() << " cells along ray\n";
            std::cout << "First blocking cell: ";
            if (det.summary.status == LOSStatus::TERRAIN_MASKED) {
                std::cout << "row=" << det.summary.blockRow
                          << " col=" << det.summary.blockCol
                          << " dist=" << det.summary.blockDistM / 1000.0 << "km\n";
            } else {
                std::cout << "none — " << losStatusStr(det.summary.status) << "\n";
            }

            // Print first 10 cells telemetry table
            std::cout << "\nRow    Col    Dist(km)  TerrElev  CorrElev  TerrSlope   TargSlope   Blocks\n";
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
    }
    else if (*los_sub) {
        Global g;
        auto res = g.lineOfVisibilityopt(srcLat, srcLon, los_ant_height, tgtLat, tgtLon, los_tgt_height);
        auto target = res.back();
        float srcTerrain = g.bin.getElevation(srcLat, srcLon);
        double srcElev = srcTerrain + los_ant_height;
           if (detailed) {
               std::cout << "=========================================\n";
                  std::cout << "Target Status : " << (target.visible ? "VISIBLE" : "BLOCKED") << "\n";
                  std::cout << "Total Distance: " << std::fixed << std::setprecision(2) << target.dist << " meters\n";
                  std::cout << "Final Angle   : " << target.angle << "°\n";
                  std::cout << "Max Obstacle  : " << target.maxAngle << "°\n";
                  std::cout << "Src elev      : " << srcElev << "m\n";
                  std::cout << "Target Elev   : " << target.elev << "m\n";
                  std::cout << "=========================================\n\n";

               size_t totalCells = res.size();
                 size_t printCount = (totalCells < 10) ? totalCells : 10;
                 size_t startIndex = totalCells - printCount;

                 std::cout << "--- Printing Last " << printCount << " Path Cells ---\n";

                 // Table Header
                 std::cout << std::setw(6)  << "Idx"
                           << std::setw(12) << "Coord (X,Y)"
                           << std::setw(10) << "Dist(m)"
                           << std::setw(9)  << "Elev(m)"
                           << std::setw(9)  << "Angle"
                           << std::setw(10) << "MaxAng"
                           << std::setw(9)  << "Visible" << "\n";
                 std::cout << std::string(65, '-') << "\n";

                 // Print rows
                 for (size_t i = startIndex; i < totalCells; ++i) {
                     const auto& cell = res[i];
                     std::string coordStr = "(" + std::to_string(cell.x) + "," + std::to_string(cell.y) + ")";

                     std::cout << std::setw(6)  << i
                               << std::setw(12) << coordStr
                               << std::setw(10) << std::fixed << std::setprecision(1) << cell.dist
                               << std::setw(9)  << cell.elev
                               << std::setw(9)  << std::setprecision(2) << cell.angle
                               << std::setw(10) << cell.maxAngle
                               << std::setw(9)  << (cell.visible ? "Yes" : "No") << "\n";
                 }
                 std::cout << std::string(65, '-') << "\n";
           }else {
               std::cout << "Target Status : " << (target.visible ? "VISIBLE" : "BLOCKED") << "\n";
               std::cout << "Total Distance: " << std::fixed << std::setprecision(2) << target.dist << " meters\n";
               std::cout << "Total cells : " << res.size() << "\n";
           }
    }else if (*run) {
        ThreadPool pool;
        GaEngine ga(pool);
        if (GaConf::is_available()) {
            GaConf cfg = GaConf::load();
            if(run_benchmark) {
                ga.run_benchmark();
            }else if (gen != 0) {
                 ga.run_native_random_simulation(gen);
            } else if (cfg.homo.is_case) {
                // Parameters (pop, rates, weights, output) already applied in
                // GaEngine constructor via applyConfig(); just call the runner.
                ga.run_homogenous(reduction_factor, png, epsg7755 ? OutputCrs::Native7755 : OutputCrs::Wgs84, single ? BandLayout::SingleBandPalette : BandLayout::Rgba);
            }
            else if (cfg.hetero.is_case) {
                // Fleet built inside run_heterogenous() from cfg.hetero.radars
                ga.run_heterogenous(reduction_factor, png, epsg7755 ? OutputCrs::Native7755 : OutputCrs::Wgs84, single ? BandLayout::SingleBandPalette : BandLayout::Rgba);
            }
            else if (cfg.autoG.is_case) {
                ga.setAutoOperators(SelectionMethod::TOURNAMENT,
                                    AutoCrossoverMethod::RANDOM_SUBSET,
                                    /*mRate=*/0.10, /*cRate=*/0.85);
                // Fleet built inside run_auto() from cfg.autoG.radars
                ga.run_auto(reduction_factor, png, epsg7755 ? OutputCrs::Native7755 : OutputCrs::Wgs84, single ? BandLayout::SingleBandPalette : BandLayout::Rgba);
            }
            else {
                std::cerr << "[ga] Config found but no is_case flag is set — "
                             "exiting\n";
                return -1;
            }
        }
    }else if (sgrid) {
        Global g;
        g.bin.init();
        Iad::get_iad_coverage_tif(g, epsg7755 ? OutputCrs::Native7755 : OutputCrs::Wgs84, single ? BandLayout::SingleBandPalette : BandLayout::Rgba);
    }
    return 0;
}
