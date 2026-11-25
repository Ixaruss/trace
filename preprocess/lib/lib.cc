#include "../lib.h"
#include <iostream>
#include <cstdlib>

using namespace std;

// The grid resolution actually in force, in metres.
//
// Cached in a function-local static: Config::load() parses JSON, so this
// must not be called per-cell in an inner loop without hoisting. The
// static initialiser is thread-safe from C++11 onward, and the resolution
// cannot change during a run, so caching for the process lifetime is safe.
//
// No fallback to the compile-time RESOLUTION: that constant is what
// produced a 100m grid silently computed over 30m data. If config isn't
// available or grid.resolution isn't set, every caller needs to know that
// loudly, not receive a plausible-looking wrong number.
int lib::grid_resolution() {
    static const int res = []() -> int {
        if (!Config::is_available()) {
            std::cerr << "[lib] ERROR: grid_resolution() requires a config — "
                         "none available. Aborting.\n";
            std::exit(EXIT_FAILURE);
        }
        Config cfg = Config::load();
        if (cfg.grid.resolution <= 0) {
            std::cerr << "[lib] ERROR: grid_resolution() — config grid.resolution "
                         "is not set. Aborting.\n";
            std::exit(EXIT_FAILURE);
        }
        return cfg.grid.resolution;
    }();
    return res;
}

// Currently unused within this codebase (grid dims come from cfg.grid.rows/
// cols directly wherever a config exists — see preprocess.cc). Kept for
// any external caller that still derives dims from extent/resolution; same
// no-fallback policy applies.
int lib::calc_rows(int resolution) {
    if (!Config::is_available()) {
        std::cerr << "[lib] ERROR: calc_rows() requires a config (grid.min_y/"
                     "max_y) — none available. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
    Config cfg = Config::load();
    return (int)ceil((cfg.grid.max_y - cfg.grid.min_y) / resolution);
}

int lib::calc_cols(int resolution) {
    if (!Config::is_available()) {
        std::cerr << "[lib] ERROR: calc_cols() requires a config (grid.min_x/"
                     "max_x) — none available. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
    Config cfg = Config::load();
    return (int)ceil((cfg.grid.max_x - cfg.grid.min_x) / resolution);
}

vector<double> lib::get_min_max_coordiantes(string path, int type) {
    GDALAllRegister();

    GDALDataset *poDS = (GDALDataset*) GDALOpenEx(path.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL);
    if (poDS == NULL) {
        std::cerr << "Open failed.\n";
    }

    OGRLayer *poLayer = poDS->GetLayer(0);
    if (poLayer == NULL) {
        std::cerr << "Layer retrieval failed.\n";
        GDALClose(poDS);
    }

    const OGRSpatialReference *poSourceSRS = poLayer->GetSpatialRef();
    if (poSourceSRS == NULL) {
        std::cerr << "Layer has no defined coordinate system.\n";
        GDALClose(poDS);
    }

    OGRSpatialReference oTargetSRS;
    oTargetSRS.importFromEPSG(4326);

    #if GDAL_VERSION_NUM >= GDAL_COMPUTE_VERSION(3,0,0)
    oTargetSRS.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    #endif

    OGRCoordinateTransformation *poCT = OGRCreateCoordinateTransformation(poSourceSRS, &oTargetSRS);
    if (poCT == NULL) {
        std::cerr << "Transformation creation failed.\n";
        GDALClose(poDS);
    }
    vector<double> coordinates(4);
    OGREnvelope oEnvelope;
    if (poLayer->GetExtent(&oEnvelope, 1) == OGRERR_NONE) {


        if(type == 0 ){
            double left_lon = oEnvelope.MinX;
            double bottom_lat = oEnvelope.MinY;
            double right_lon = oEnvelope.MaxX;
            double top_lat = oEnvelope.MaxY;
            if (poCT->Transform(1, &left_lon, &bottom_lat) && poCT->Transform(1, &right_lon, &top_lat)) {
                coordinates[0] = left_lon;
                coordinates[1] = right_lon;
                coordinates[2] = bottom_lat;
                coordinates[3] = top_lat;
            } else {
                std::cerr << "Coordinate transform failed.\n";
            }
        }else {
            coordinates[0] = oEnvelope.MinX;
            coordinates[1] = oEnvelope.MaxX;
            coordinates[2]= oEnvelope.MinY;
            coordinates[3] = oEnvelope.MaxY;
        }
    }

    OGRCoordinateTransformation::DestroyCT(poCT);
    GDALClose(poDS);
    return coordinates;
}
