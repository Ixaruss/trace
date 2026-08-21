#include "coverage.h"
#include "config.h"

using namespace std;
namespace {
bool gridToLatLon(Global& g, int r, int c, double& lat, double& lon) {
    static OGRCoordinateTransformation* inv = nullptr;
    if (!inv) {
        OGRSpatialReference src, dst;
        src.importFromEPSG(7755);
        dst.importFromEPSG(4326);
        src.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
        dst.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
        inv = OGRCreateCoordinateTransformation(&src, &dst);
        if (!inv) return false;
    }
    double x = g.conf.grid.origin_x + (double)c * RESOLUTION;
    double y = g.conf.grid.origin_y - (double)r * RESOLUTION;
    if (!inv->Transform(1, &x, &y)) return false;
    lon = x; lat = y;
    return true;
}
}

void Iad::get_iad_coverage_tif(Global &g, OutputCrs opt, BandLayout band, string output_location) {
    IadConf conf;
    if (IadConf::is_available()) conf = IadConf::load();
    else { cerr << "iadconf.json not found -- exiting\n"; return; }

    if (conf.radars.empty()) { cerr << "[iad] no radars in config\n"; return; }

    const int rows = g.conf.grid.rows;
    const int cols = g.conf.grid.cols;

    vector<GridPoint>    placements;
    vector<IadRadarConf> used;
    placements.reserve(conf.radars.size());
    used.reserve(conf.radars.size());

    for (const auto& x : conf.radars) {
        int r = -1, c = -1;
        const bool ok = g.bin.getGridCoords(x.lat, x.lon, r, c);

        cout << "[iad] " << (x.label.empty() ? "<unlabelled>" : x.label)
             << " lat=" << x.lat << " lon=" << x.lon
             << " -> (" << r << "," << c << ")\n";

        if (!ok || r < 0 || c < 0 || r >= rows || c >= cols) {
            cerr << "[iad] SKIP: maps outside the " << rows << "x" << cols
                 << " grid. Check that the DEM actually covers this lat/lon "
                    "and that both are in EPSG:4326 degrees.\n";
            continue;
        }

        placements.push_back({r, c, 0});
        used.push_back(x);
    }

    if (placements.empty()) {
       double latNW, lonNW, latSE, lonSE;
if (gridToLatLon(g, 0, 0, latNW, lonNW) &&
    gridToLatLon(g, rows - 1, cols - 1, latSE, lonSE)) {
    cerr << "[iad] grid corners: NW (" << latNW << "," << lonNW
         << ")  SE (" << latSE << "," << lonSE << ")\n";
}
return;
    }

    calculate_radar_coverage_multi(
        g,
        placements,
        used,
        opt,
        band,
        conf.target,
        conf.output_path,
        /*reduction=*/0,false,
        conf.use_curv);
}
