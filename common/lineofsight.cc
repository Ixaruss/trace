#include "../common.h"
#include <cmath>
#include <cstddef>
#include <chrono>
#include <gdal.h>
#include <ogr_geometry.h>
#include <ogr_spatialref.h>
#include <ogr_srs_api.h>
using namespace std;

// int get_dist(double lat0, double long0, double lat1,double long1){
//     GDALAllRegister();
//     OGRSpatialReference src;
//     src.importFromEPSG(7755);
//     src.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
//     OGRSpatialReference target;
//     src.importFromEPSG(7755);
//     src.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
//     OGRCoordinateTransformation *poct = OGRCreateCoordinateTransformation(&src,&target);
//     if(poct ==nullptr) {
//         std::cerr << "transformation failed";
//         return 1;
//     }
//     OGRPoint p1(lat0,long0);
//     OGRPoint p2(lat1,long1);
//     if(p1.transform((poct)))
// }
bool Global::lineOfSight(double lat0, double long0, int h0,
                          double lat1, double long1, int h1) {
    // no local Global — use this->bin which is already loaded
    int x0, y0, x1, y1;
    auto start = chrono::steady_clock::now();

    if (!bin.getGridCoords(lat0, long0, x0, y0) ||
        !bin.getGridCoords(lat1, long1, x1, y1)) {
        cout << "no\n";
        return false;
    }

    double dx       = x1 - x0;
    double dy       = y1 - y0;
    // FIX: was `sqrt(dx*dx+dy*dy) / RESOLUTION`. dx/dy are already in CELL
    // units (x0,y0,x1,y1 come from getGridCoords), so this needs to
    // MULTIPLY by the cell size to get metres, not divide by it — dividing
    // produced a dimensionless, meaningless number. Also cell_m (what the
    // grid was actually built at), not RESOLUTION.
    double distance = sqrt(dx*dx + dy*dy) * bin.cell_m;

    int src_elev  = bin.getElevation(lat0, long0) + h0;
    int targ_elev = bin.getElevation(lat1, long1) + h1;

    cout << "src: "    << src_elev  << "\n";
    cout << "target: " << targ_elev << "\n";

    auto cells = Global::bresenham(x0, y0, x1, y1);
    double delta = (double)(targ_elev - src_elev) / cells.size();

    cout << "delta: " << delta
         << "  cells: " << cells.size()
         << "  dist: "  << distance << "\n";

    for (size_t i = 0; i < cells.size(); i++) {
        int lineh = src_elev + (int)(delta * i);
        // FIX: was a direct elev_matrix[row*conf.grid.cols+col] index —
        // that's the FULL-grid layout, only correct by coincidence in
        // load_full mode (data_window == full grid there). In tiled mode
        // elev_matrix holds only data_window, so this read went out of
        // bounds. getElevationRC() indexes via data_window and is correct
        // in both layouts.
        int cellh = bin.getElevationRC(cells[i].first, cells[i].second);
        if (cellh <= -9999) {
            cerr << "[los] WARN: cell (" << cells[i].first << ", "
                 << cells[i].second << ") is not resident or has no data — "
                    "treating as blocking\n";
            return false;
        }
        if (cellh > lineh) return false;
    }

    auto end   = chrono::steady_clock::now();
    int  total = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    cout << "total time: " << total << "ms\n";
    return true;
}
