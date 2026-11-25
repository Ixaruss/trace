#include "../bin.h"
#include <vector>
#include "../../multi-radar-deployment/h5.h"
using namespace std;
using namespace lib;

bool Bin::getGridCoords(double lat, double lng, int& out_r, int& out_c) {
    // CHANGED: was init({ELEVATION}) — that's window-RESOLVING and warned
    // on every call once a real window was already set via load()/init()
    // elsewhere. ensureLoaded() only loads under the existing window.
    ensureLoaded({ ELEVATION });
    if (!elevation_loaded || !poCT) return false;

    double x = lng, y = lat;
    if (!poCT->Transform(1, &x, &y)) return false;

    // cell_m, not RESOLUTION: with 30m data and RESOLUTION=100 this used to
    // land every lat/lng in the wrong cell (a 100m grid laid over 30m data),
    // which is why the radar source was reported outside a window that had
    // been built — correctly, at 30m — centred on it.
    out_c = (int)std::floor((x - conf.grid.origin_x) / cell_m);
    out_r = (int)std::floor((conf.grid.origin_y - y)  / cell_m);

    // NOTE: this validates against the FULL grid, not the resident window —
    // a coordinate can be a legitimate grid cell while not being loaded.
    // Use data_window.contains(r, c) to test residency.
    return (out_r >= 0 && out_r < conf.grid.rows &&
            out_c >= 0 && out_c < conf.grid.cols);
}

bool Bin::getLatLngFromGrid(int r, int c, double& out_lat, double& out_lng) {
    ensureLoaded({ ELEVATION });
    if (!elevation_loaded || !poCT) return false;

    if (r < 0 || r >= conf.grid.rows || c < 0 || c >= conf.grid.cols)
        return false;

    // Use cell center, since the forward transform floors to the containing cell
    double x = conf.grid.origin_x + (c + 0.5) * cell_m;
    double y = conf.grid.origin_y - (r + 0.5) * cell_m;

    OGRCoordinateTransformation* poCTInv = poCT->GetInverse();
    if (!poCTInv) return false;

    bool ok = poCTInv->Transform(1, &x, &y);
    OCTDestroyCoordinateTransformation(poCTInv);
    if (!ok) return false;

    out_lng = x;
    out_lat = y;
    return true;
}

// ── window-aware raw accessors ─────────────────────────────────────────────
// elev_matrix / slope_matrix hold only op_window, so global (r,c) must go
// through data_window.index(). Indexing them with (r * grid.cols + c) is only
// correct in the FULL case and silently reads the wrong cell otherwise.
int16_t Bin::getElevationRC(int r, int c) {
    if (!elevation_loaded) return -9999;
    if (!data_window.contains(r, c)) return -9999;
    return elev_matrix[data_window.index(r, c)];
}

// Window-aware constraints read. vec_matrix is now windowed like
// elev_matrix, so (r * VEC_COLS + c) is no longer a valid index — that
// only worked while the whole grid was resident. Returns 0 (no bits set)
// for a cell outside the loaded window, which reads as "unconstrained";
// callers that care about the difference should check data_window.contains.
uint8_t Bin::getConstraintRC(int r, int c) {
    if (!vector_loaded) return 0;
    if (!data_window.contains(r, c)) return 0;
    return vec_matrix[data_window.index(r, c)];
}

float Bin::getSlopeValueRC(int r, int c) {
    if (!slope_loaded) return -9999.0f;
    if (!data_window.contains(r, c)) return -9999.0f;
    return slope_matrix[data_window.index(r, c)];
}

float Bin::getSlopeValue(double lat, double lng) {
    ensureLoaded({ SLOPE_ });
    int r, c;
    if (!getGridCoords(lat, lng, r, c)) return -9999.0f;
    return getSlopeValueRC(r, c);
}

int16_t Bin::getElevation(double lat, double lng) {
    // Was a private re-implementation of lat/lng -> cell (with its own
    // RESOLUTION divide and its own op_window indexing). Now one conversion
    // (getGridCoords, cell_m) and one lookup (getElevationRC, data_window),
    // shared with every other query — correct in tiled and flat mode alike.
    int r, c;
    if (!getGridCoords(lat, lng, r, c)) return -9999;
    if (!data_window.contains(r, c)) {
        cerr << "[bin] WARN: (" << lat << ", " << lng << ") -> cell ("
             << r << ", " << c << ") is not resident - rows ["
             << data_window.r0 << " -> " << (data_window.r0 + data_window.rows)
             << ") cols [" << data_window.c0 << " -> " << (data_window.c0 + data_window.cols)
             << ")\n";
        return -9999;
    }
    return getElevationRC(r, c);
}

bool Bin::getSlope(double lat, double lng) {
    int r, c;
    if (!getGridCoords(lat, lng, r, c)) return false;
    return (getConstraintRC(r, c) & (1 << BIT_SLOPE)) != 0;
}

bool Bin::getStatusRoad(double lat, double lng) {
    int r, c;
    if (!getGridCoords(lat, lng, r, c)) return false;
    return (getConstraintRC(r, c) & (1 << BIT_ROAD)) != 0;
}

bool Bin::getStatusRail(double lat, double lng) {
    int r, c;
    if (!getGridCoords(lat, lng, r, c)) return false;
    return (getConstraintRC(r, c) & (1 << BIT_RAIL)) != 0;
}

bool Bin::getStatusWaterArea(double lat, double lng) {
    int r, c;
    if (!getGridCoords(lat, lng, r, c)) return false;
    return (getConstraintRC(r, c) & (1 << BIT_WATER_AREA)) != 0;
}

bool Bin::getStatusWaterLine(double lat, double lng) {
    int r, c;
    if (!getGridCoords(lat, lng, r, c)) return false;
    return (getConstraintRC(r, c) & (1 << BIT_WATER_LINE)) != 0;
}

bool Bin::getStatusIB(double lat, double lng) {
    int r, c;
    if (!getGridCoords(lat, lng, r, c)) return false;
    return (getConstraintRC(r, c) & (1 << BIT_IB)) != 0;
}

vector<int> Bin::getStatusAll(double lat, double lng) {
    ensureLoaded({});   // vec + elev
    vector<int> res;
    int r, c;
    if (!getGridCoords(lat, lng, r, c)) return res;

    uint8_t cell = getConstraintRC(r, c);
    res.push_back((cell >> BIT_ROAD)       & 1);
    res.push_back((cell >> BIT_RAIL)       & 1);
    res.push_back((cell >> BIT_WATER_AREA) & 1);
    res.push_back((cell >> BIT_WATER_LINE) & 1);
    res.push_back((cell >> BIT_IB)         & 1);
    res.push_back((cell >> BIT_SLOPE)      & 1);
    return res;
}

bool Bin::latlon_to_projected(double lat, double lon, double& px, double& py) {
    px = lon; py = lat;
    if (!poCT->Transform(1, &px, &py)) {
        cerr << "[proj] ERROR: could not transform ("
             << lat << ", " << lon << ")\n";
        return false;
    }
    return true;
}

bool Bin::point_in_polygon(double px, double py,
                            const vector<pair<double,double>>& proj_poly) {
    int  n      = (int)proj_poly.size();
    bool inside = false;
    for (int i = 0, j = n - 1; i < n; j = i++) {
        double xi = proj_poly[i].first,  yi = proj_poly[i].second;
        double xj = proj_poly[j].first,  yj = proj_poly[j].second;
        if ((yi > py) != (yj > py)) {
            double x_intersect = (xj - xi) * (py - yi) / (yj - yi) + xi;
            if (px < x_intersect) inside = !inside;
        }
    }
    return inside;
}

BoundingBox Bin::getMinimumBoundingBox() {
    if (conf.area_op.size() < 3) {
        cerr << "[mbr] ERROR: area_op must have at least 3 points\n";
        return { 0, 0, 0, 0 };
    }
    double min_lat =  1e9, max_lat = -1e9;
    double min_lon =  1e9, max_lon = -1e9;
    for (auto& [lat, lon] : conf.area_op) {
        min_lat = min(min_lat, lat);
        max_lat = max(max_lat, lat);
        min_lon = min(min_lon, lon);
        max_lon = max(max_lon, lon);
    }
    cout << "[mbr] BoundingBox → lat: [" << min_lat << " → " << max_lat
         << "] lon: [" << min_lon << " → " << max_lon << "]\n";
    return { min_lat, max_lat, min_lon, max_lon };
}

#ifndef _FILE_OFFSET_BITS
  #define _FILE_OFFSET_BITS 64
#endif

// ── feasible grid ──────────────────────────────────────────────────────────
// Scan bounds now come from the Bin's resident window rather than being
// recomputed here, so this can never scan cells that aren't loaded:
//
//   FULL       → whole grid (window is the whole grid)
//   UNBOUNDED  → the square window around the center point
//   AOI        → the window (= polygon AABB), plus a per-cell polygon test
//
// `opt == 1` still forces a full-grid scan for backward compatibility, but
// it will refuse if the Bin only holds a partial window.
FeasibleGrid Bin::getFeasibleGrid(Bin& b, short op, short opt, short allCells) {
    FeasibleGrid result;

    if (!b.vector_loaded) {
        cerr << "[feasible] ERROR: vec_matrix not loaded\n";
        return result;
    }

    // ── grid bounds ────────────────────────────────────────────────────────
    int r_min, c_min, r_max, c_max;
    bool use_polygon = false;
    vector<pair<double,double>> proj_poly;

    if (!b.window_valid) {
        cerr << "[feasible] ERROR: no window resolved — call Global::load() "
                "(or Bin::init()) first\n";
        return result;
    }

    // vec_matrix must describe the same grid as the elevation/config grid.
    // If map.bin/map_tiles were built at a different resolution, fail loudly
    // here rather than scanning garbage.
    if (b.GRID_ROWS > 0 && b.GRID_COLS > 0 &&
        (b.VEC_ROWS != b.GRID_ROWS || b.VEC_COLS != b.GRID_COLS)) {
        cerr << "[feasible] ERROR: constraints grid does not match the elevation/config grid\n"
             << "  map (vec)   : " << b.VEC_COLS  << " x " << b.VEC_ROWS  << "\n"
             << "  config/elev : " << b.GRID_COLS << " x " << b.GRID_ROWS << "\n"
             << "  Rebuild the constraints (preprocess) at the current resolution.\n";
        return result;
    }

    // ── scan bounds ──────────────────────────────────────────────────────
    // ONE code path for tiled and load_full. Bounds come from op_window (the
    // region of interest); every cell READ goes through getConstraintRC,
    // which indexes via data_window (the resident layout). In load_full mode
    // data_window is the whole grid, so those reads are exactly the legacy
    // r*VEC_COLS+c — no separate branch needed, and none to drift.
    const bool force_full = (opt == 1);

    if (force_full) {
        // A full-grid scan needs the full grid RESIDENT, which is a property
        // of data_window, not op_window: in load_full mode the whole grid is
        // in memory even for an AOI run, so opt=1 is legitimate there.
        if (!b.data_window.full(b.GRID_ROWS, b.GRID_COLS)) {
            cerr << "[feasible] ERROR: opt=1 asks for a full-grid scan but only a "
                 << b.data_window.cols << " x " << b.data_window.rows
                 << " window is resident — load with FULL, or use block_size=0\n";
            return result;
        }
        r_min = 0;           c_min = 0;
        r_max = b.GRID_ROWS; c_max = b.GRID_COLS;
        cout << "[feasible] opt=1 — scanning entire grid: "
             << b.GRID_COLS << " x " << b.GRID_ROWS << "\n";
    } else {
        r_min = b.op_window.r0;
        c_min = b.op_window.c0;
        r_max = b.op_window.r0 + b.op_window.rows;
        c_max = b.op_window.c0 + b.op_window.cols;

        switch (b.op_type) {
            case FULL:
                cout << "[feasible] FULL — scanning entire grid: "
                     << (c_max - c_min) << " x " << (r_max - r_min) << "\n";
                break;
            case UNBOUNDED:
                cout << "[feasible] UNBOUNDED — scanning rows ["
                     << r_min << " → " << r_max << ") cols ["
                     << c_min << " → " << c_max << ")\n";
                break;
            case AOI: {
                proj_poly.reserve(b.conf.area_op.size());
                double pminx = 1e30, pmaxx = -1e30, pminy = 1e30, pmaxy = -1e30;
                for (auto& [lat, lon] : b.conf.area_op) {
                    double px, py;
                    if (!b.latlon_to_projected(lat, lon, px, py)) {
                        cerr << "[feasible] ERROR: failed to project polygon vertex\n";
                        return result;
                    }
                    proj_poly.push_back({ px, py });
                    pminx = min(pminx, px); pmaxx = max(pmaxx, px);
                    pminy = min(pminy, py); pmaxy = max(pmaxy, py);
                }
                use_polygon = true;

                // Tighten the scan to the polygon's own bounding box. op_window
                // carries a margin (the asset range, e.g. 75km) so radars placed
                // near the AOI edge can see outward — but no cell in that margin
                // can pass the polygon test, so scanning it only burns time.
                // Your 75km-margin AOI window is 62.6M cells; its polygon box
                // is ~8.4M. Result is identical — those cells were all
                // "outside polygon" anyway.
                const double cm = (double)b.cell_m;
                const int pc0 = (int)std::floor((pminx - b.conf.grid.origin_x) / cm);
                const int pc1 = (int)std::ceil ((pmaxx - b.conf.grid.origin_x) / cm) + 1;
                const int pr0 = (int)std::floor((b.conf.grid.origin_y - pmaxy) / cm);
                const int pr1 = (int)std::ceil ((b.conf.grid.origin_y - pminy) / cm) + 1;
                r_min = max(r_min, pr0); r_max = min(r_max, pr1);
                c_min = max(c_min, pc0); c_max = min(c_max, pc1);

                cout << "[feasible] AOI — polygon box rows [" << r_min << " → " << r_max
                     << ") cols [" << c_min << " → " << c_max << ")"
                     << "  (region of interest rows [" << b.op_window.r0 << " → "
                     << (b.op_window.r0 + b.op_window.rows) << ") incl. margin)\n";
                break;
            }
        }
    }

    if (r_max <= r_min || c_max <= c_min) {
        cerr << "[feasible] ERROR: empty scan region (rows " << r_min << "→" << r_max
             << ", cols " << c_min << "→" << c_max << ")\n";
        return result;
    }

    bool need_infra = b.flag_road || b.flag_rail;

    // Resolution the grid was BUILT at, hoisted out of the per-cell loop.
    // Converts (r, c) back to projected metres for the polygon test. This
    // MUST match the resolution used to build the window: if it's the
    // compile-time RESOLUTION (100) against 30m data, every cell centre is
    // placed 3.33x too far from the origin and the entire scan reports
    // "outside polygon" — which is exactly the zero-feasible-cells failure
    // that then crashes the GA.
    const double CELL = (double)b.cell_m;

    result.origin_r = r_min;
    result.origin_c = c_min;
    result.height   = r_max - r_min;
    result.width    = c_max - c_min;

    long long outside_count  = 0;
    long long water_count    = 0;
    long long ib_count       = 0;
    long long slope_count    = 0;
    long long no_infra_count = 0;
    const long long scanned  = (long long)(r_max - r_min) * (c_max - c_min);

    cout << "[feasible] Scanning " << scanned << " cells @ " << b.cell_m << " m...\n";

    for (int r = r_min; r < r_max; r++) {
        for (int c = c_min; c < c_max; c++) {

            if (use_polygon) {
                double cx = b.conf.grid.origin_x + (c + 0.5) * CELL;
                double cy = b.conf.grid.origin_y - (r + 0.5) * CELL;
                if (!Bin::point_in_polygon(cx, cy, proj_poly)) {
                    outside_count++; continue;
                }
            }

            // data_window indexing — correct in both tiled and load_full mode
            const uint8_t cell = b.getConstraintRC(r, c);

            if (!allCells) {
                if (b.flag_water_line && (cell & (1 << BIT_WATER_LINE))) { water_count++; continue; }
                if (b.flag_water_area && (cell & (1 << BIT_WATER_AREA))) { water_count++; continue; }
                if (b.flag_ib         && (cell & (1 << BIT_IB)))         { ib_count++;    continue; }
                if (b.flag_slope      && (cell & (1 << BIT_SLOPE)))      { slope_count++; continue; }

                if (need_infra) {
                    bool road_ok = b.flag_road && (cell & (1 << BIT_ROAD));
                    bool rail_ok = b.flag_rail && (cell & (1 << BIT_RAIL));
                    if (!road_ok && !rail_ok) { no_infra_count++; continue; }
                }
            }
            result.cells.push_back({ r, c, cell });
        }
    }
    if (op == WRITE)
        writeBin("./bin/areaofoperation.bin", result.cells);

    cout << "[feasible] Outside polygon : " << outside_count  << "\n";
    cout << "[feasible] Water cells     : " << water_count    << "\n";
    cout << "[feasible] IB cells        : " << ib_count       << "\n";
    cout << "[feasible] Slope rejected  : " << slope_count    << "\n";
    cout << "[feasible] No infra        : " << no_infra_count << "\n";
    cout << "[feasible] Feasible        : " << result.cells.size()
         << " / " << scanned << "\n";

    // Empty is fatal downstream in EITHER mode (the GA indexes cells[0] and
    // draws from an empty distribution), so this is no longer gated on
    // load_full. Name the likeliest cause from the counters.
    if (result.cells.empty()) {
        cerr << "[feasible] WARN: no cells survived — callers must not proceed.\n";
        if (use_polygon && outside_count == scanned)
            cerr << "  EVERY cell fell outside the AOI polygon. The polygon and the\n"
                 << "  scan disagree about geometry — almost always a resolution\n"
                 << "  mismatch (cell_m=" << b.cell_m << ", config grid.resolution="
                 << b.conf.grid.resolution << "). Rebuild after changing either.\n";
        else if (no_infra_count > 0 && no_infra_count + outside_count == scanned)
            cerr << "  Every in-polygon cell lacked road/rail — the road/rail bits look\n"
                 << "  empty; check preprocess's 'final hits' for those layers.\n";
        else
            cerr << "  flags: road=" << b.flag_road << " rail=" << b.flag_rail
                 << " water_area=" << b.flag_water_area << " water_line=" << b.flag_water_line
                 << " ib=" << b.flag_ib << " slope=" << b.flag_slope << "\n";
    }

    return result;
}

int Bin::getTerrianResolution() {
    if (cell_m > 0) return cell_m;
    // No silent fallback to the compile-time RESOLUTION: a caller asking
    // for the grid's resolution and getting a plausible-looking wrong
    // number back (100 instead of 30) is exactly the bug class this pass
    // removes. Fail loudly instead of returning a guess.
    std::cerr << "[bin] ERROR: getTerrianResolution() has no resolution — "
                 "config unavailable or grid.resolution unset. Aborting.\n";
    std::exit(EXIT_FAILURE);
}

void Bin::getonlyFeasibleGrid(Bin& b){
    auto x = getFeasibleGrid(b).cells;
    cout <<"bin size: " << x.size() <<endl;
    for (int i =0; i<x.size();i++) {
        cout << "cell " << i << ": " << "coordinates: " << x[i].c << x[i].r << endl;
    }
     writeBin("./pso.bin", x);
}

void Bin::write_heatmap_bin_h5 (){
    Bin b;
    // AOI: only the operational polygon's bounding box is brought in
    b.init({CONSTRAINTS}, AOI);
    auto bin = Bin::getFeasibleGrid(b,0,0,WITHALL);
    auto h5_path = "output/heatmap/aoi_cells.h5";
    H5E_auto2_t old_func;
    void*       old_data;
    H5Eget_auto2(H5E_DEFAULT, &old_func, &old_data);
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);          // silence
    hid_t raw_fid = H5Fopen(h5_path, H5F_ACC_RDWR, H5P_DEFAULT);
    H5Eset_auto2(H5E_DEFAULT, old_func, old_data);        // restore

    if (raw_fid < 0)
        raw_fid = H5Fcreate(h5_path,
                             H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (raw_fid < 0)
        throw std::runtime_error(
            "write_h5_indiviual: cannot open or create file ");

    h5_lib::H5FileGuard fid{ raw_fid };

    H5G_info_t root_info{};
    H5Gget_info(fid.id, &root_info);
    const int run_n = static_cast<int>(root_info.nlinks) + 1;

    h5_lib::H5GroupGuard data{
        H5Gcreate2(fid.id, "aoi",
                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) };
    if (data.id < 0)
        throw std::runtime_error("write_h5_indiviual: failed to create 'data' group");
    vector<int> r,c;
    for (auto x: bin.cells) {
        r.push_back(x.r);
        c.push_back(x.c);
    }
    h5_lib::write_i32(data.id, "r",        r);
    h5_lib::write_i32(data.id, "c",        c);
}
pair<double,double> Bin::getMapOrigins() {
    if (!Config::is_available()) {
        std::cerr << "[bin] ERROR: getMapOrigins() requires a config — none "
                     "available. Aborting rather than returning the "
                     "compiled-in VEC_ORIGIN_X/Y placeholder.\n";
        std::exit(EXIT_FAILURE);
    }
    return { conf.grid.origin_x, conf.grid.origin_y };
}
