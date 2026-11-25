#include "../bin.h"
#include <cstdint>
#include <filesystem>
#include <string>
using namespace std;
using namespace lib;

// Cells of slack added around a requested window. Engines clamp their
// working radius to the distance from the source cell to the nearest edge,
// so an exactly-sized window costs a cell of range.
static const int WINDOW_MARGIN_CELLS = 2;

// ── shared tiled/flat loader for elev, slope and vec ──────────────────────
// All three matrices (elev.bin/elev_tiles, slope.bin/slope_tiles,
// map.bin/map_tiles) are written by the SAME tiling decision — cfg.grid.
// block_size — and must be loaded by the same decision, with the same
// fallback-and-explain behaviour if the preferred layout is missing or
// doesn't match the grid. This one function is that decision; load_elev,
// load_slope and load_vec each just plug in their element type and path.
template <typename T>
static bool read_tiled(Bin* b, const string& path, const char* label,
                       const GridWindow& win, vector<T>& out, string& why) {
    const string dir = tilegrid::tile_dir_for(path);
    tilegrid::Grid g;
    if (!g.discover(dir)) { why = "no complete tile folder at " + dir; return false; }
    g.print(label);
    if (g.rows != b->GRID_ROWS || g.cols != b->GRID_COLS) {
        why = "tile folder is " + to_string(g.cols) + " x " + to_string(g.rows) +
              " but the grid is " + to_string(b->GRID_COLS) + " x " + to_string(b->GRID_ROWS) +
              " (built at a different resolution/extent?)";
        return false;
    }
    if (g.item_size != (int)sizeof(T)) {
        why = "tiles are " + to_string(g.item_size) + " bytes/cell, expected " + to_string(sizeof(T));
        return false;
    }
    if (!tilegrid::read_window<T>(g, win.r0, win.c0, win.rows, win.cols, out)) {
        out.clear();
        why = "tile read failed (see the [tiles] ERROR line above)";
        return false;
    }
    return true;
}

template <typename T>
static bool read_flat(Bin* b, const string& path, const GridWindow& win,
                      vector<T>& out, string& why) {
    namespace fs = std::filesystem;
    error_code ec;
    if (!fs::is_regular_file(path, ec)) { why = "no flat file at " + path; return false; }
    const uintmax_t have = fs::file_size(path, ec);
    const uintmax_t need = (uintmax_t)b->GRID_ROWS * (uintmax_t)b->GRID_COLS * sizeof(T);
    if (ec || have != need) {
        why = path + " is " + to_string((double)have / 1e6) + " MB but the " +
              to_string(b->GRID_COLS) + " x " + to_string(b->GRID_ROWS) + " grid needs " +
              to_string((double)need / 1e6) + " MB (built at a different resolution/extent?)";
        return false;
    }
    if (!tilegrid::read_window_flat<T>(path, b->GRID_ROWS, b->GRID_COLS,
                                        win.r0, win.c0, win.rows, win.cols, out)) {
        out.clear();
        why = "flat read failed (see the [tiles] ERROR line above)";
        return false;
    }
    return true;
}

// Tries the layout cfg.grid.block_size says to prefer, falls back to the
// other one on failure, and prints exactly which layout was actually used
// (or, on total failure, why both attempts failed). Returns true iff `out`
// now holds the requested window.
template <typename T>
static bool load_matrix(Bin* b, const string& path, const char* label,
                        int block_size, const GridWindow& win, vector<T>& out) {
    const bool prefer_tiles = (block_size > 0);
    string why1, why2; bool ok = false; const char* used = "";

    cout << "[bin] " << label << ": block_size=" << block_size << " -> "
         << (prefer_tiles ? "tiled layout" : "single flat bin") << "\n";

    if (prefer_tiles) {
        if ((ok = read_tiled<T>(b, path, label, win, out, why1))) used = "tiles";
        else {
            cerr << "[bin] WARN: " << label << " tiles unusable - " << why1
                 << "\n[bin]       trying the flat file instead\n";
            if ((ok = read_flat<T>(b, path, win, out, why2))) used = "flat (fallback)";
        }
    } else {
        if ((ok = read_flat<T>(b, path, win, out, why1))) used = "flat";
        else {
            cerr << "[bin] WARN: " << label << " flat file unusable - " << why1
                 << "\n[bin]       trying tiles instead\n";
            if ((ok = read_tiled<T>(b, path, label, win, out, why2))) used = "tiles (fallback)";
        }
    }

    if (!ok) {
        out.clear();
        cerr << "[bin] ERROR: " << label << " matrix could not be loaded.\n"
             << "  " << (prefer_tiles ? "tiles" : "flat ") << " : " << why1 << "\n"
             << "  " << (prefer_tiles ? "flat " : "tiles") << " : " << why2 << "\n"
             << "  Re-run preprocess/build_elevation_matrix with the same grid.block_size you"
                " are running with (0 = single flat bin, >0 = tiled folder).\n";
        return false;
    }
    cout << "[bin] " << label << " loaded [" << used << "] - window " << win.cols << " x "
         << win.rows << " (" << (out.size() * sizeof(T)) / 1e6 << " MB)\n";
    return true;
}

// ── vector matrix (map.bin / map_tiles), windowed ─────────────────────────
// The constraints matrix can sit on disk in TWO layouts, and both must keep
// working end to end (preprocess writes → this loads → GA consumes):
//
//     cfg.grid.block_size  > 0   →  map_tiles/  folder of fixed-size tiles
//     cfg.grid.block_size == 0   →  one flat map.bin, exactly as before tiling
//
// preprocess picks the layout from block_size; this loader picks it from the
// SAME value, so a leftover folder/file from an earlier run of the OTHER
// layout can never shadow the one that was just built. (Previously this
// looked for map_tiles/ first no matter what block_size said — a stale or
// mismatched tile folder made the load fail even though a perfectly good
// flat map.bin was sitting next to it, and vice versa.)
//
// If the preferred layout is missing or doesn't match the grid, the other
// one is tried before giving up, and every failure says why.
//
// vec_matrix is WINDOW-relative like elev_matrix. Every read must go through
// data_window.index(r, c) — see Bin::getConstraintRC.

// Tiled read. `why` is filled on failure.
static void load_vec(Bin* b, const string& param_path, const GridWindow& win) {
    if (b->vector_loaded) return;

    if (!Config::is_available()) {
        cerr << "[bin] ERROR: load_vec requires a config (fcase flags, "
                "grid.block_size) — none available.\n";
        return;
    }
    auto conf          = Config::load();
    b->flag_road       = conf.fcase.roads;
    b->flag_rail       = conf.fcase.rails;
    b->flag_water_area = conf.fcase.waterbodies;
    b->flag_water_line = conf.fcase.waterlines;
    b->flag_slope      = conf.fcase.slope;
    b->flag_ib         = conf.fcase.ib;
    const int block_size = conf.grid.block_size;

    if (b->GRID_ROWS <= 0 || b->GRID_COLS <= 0) {
        cerr << "[bin] ERROR: load_vec called before the grid dimensions were "
                "resolved - go through Bin::init()/Global::load()\n";
        return;
    }

    // VEC_ROWS/COLS describe the FULL grid the window is cut from, and that
    // is by definition the config grid: load_matrix's size checks reject any
    // file/tile set that isn't. Assigned unconditionally (not "only if 0")
    // so a value left over from a different resolution can't survive into
    // getFeasibleGrid's grid-match guard.
    b->VEC_ROWS = b->GRID_ROWS;
    b->VEC_COLS = b->GRID_COLS;

    b->vector_loaded = load_matrix<uint8_t>(b, param_path, "map", block_size, win, b->vec_matrix);
}

// load_vec_full / load_elev_full removed: load_full mode now loads through
// the same load_vec/load_elev/load_slope as tiled mode, with data_window set
// to the whole grid (see Bin::init). They were duplicates that had already
// drifted — neither loaded slope, and load_elev_full exit()ed the process on
// a missing file where the windowed path reports and continues.

// ── window resolution ──────────────────────────────────────────────────────
// Turns an OperationType into a concrete sub-rectangle of the global grid.
// Nothing is read here — this only decides what WILL be read.
bool Bin::computeWindow(OperationType type,
                        double center_lat, double center_lng, double range_km,
                        GridWindow& out) {
    // FIX: conf was previously read as whatever Bin() captured at
    // construction time. That's fine for a Global built fresh per request,
    // but with Global::load()/reload() a single Bin can now be reused
    // across many requests with DIFFERENT AOI polygons — reading a stale
    // conf.area_op here would silently compute the window for the WRONG
    // polygon. Refreshing here keeps AOI correct regardless of how long
    // this Bin has been alive.
    // Config is required — every grid dimension, origin, resolution and
    // block_size decision downstream comes from it. The old ELEV_ROWS/
    // ELEV_COLS/ELEV_ORIGIN_X/Y/RESOLUTION fallback let the whole pipeline
    // run silently against compiled-in placeholders (a 100m grid) while
    // the actual data was 30m — which is what produced a window computed
    // against the wrong grid. Fail loudly instead.
    if (!Config::is_available()) {
        cerr << "[bin] ERROR: computeWindow requires a config (grid dims, "
                "origin, resolution, block_size) — none is available. "
                "Refusing to fall back to compiled-in defaults.\n";
        return false;
    }
    conf = Config::load();

    const Config cfg = conf;
    if (cfg.grid.rows <= 0 || cfg.grid.cols <= 0) {
        cerr << "[bin] ERROR: config grid.rows/grid.cols are not set (or <= 0) — "
                "run build_elevation_matrix()/preprocess() first.\n";
        return false;
    }
    if (cfg.grid.resolution <= 0) {
        cerr << "[bin] ERROR: config grid.resolution is not set (or <= 0).\n";
        return false;
    }

    const int    grid_rows = cfg.grid.rows;
    const int    grid_cols = cfg.grid.cols;
    const double ox        = cfg.grid.origin_x;
    const double oy        = cfg.grid.origin_y;
    const int    resol     = cfg.grid.resolution;

    GRID_ROWS = grid_rows;
    GRID_COLS = grid_cols;

    // Authoritative geometry + layout for everything downstream.
    cell_m    = resol;
    load_full = (cfg.grid.block_size == 0);

    // RESOLUTION (lib.h) is a compile-time constant. Every geometry
    // computation in this codebase has been converted to use cell_m
    // (this value) instead — this warning exists only in case a stray
    // reference to RESOLUTION survives somewhere it shouldn't.
    if (resol != (int)RESOLUTION) {
        static bool warned = false;
        if (!warned) {
            warned = true;
            cerr << "[bin] NOTE: compile-time RESOLUTION=" << (int)RESOLUTION
                 << " differs from config grid.resolution=" << resol
                 << ". This is harmless AS LONG AS nothing still reads "
                    "RESOLUTION directly for geometry — everything in this "
                    "codebase should be using Bin::cell_m.\n";
        }
    }

    if (type == FULL) {
        out = { 0, 0, grid_rows, grid_cols };
        cout << "[bin] window FULL — " << grid_cols << " x " << grid_rows << "\n";
        return true;
    }

    // projected-space bounds of the requested region
    double min_x, max_x, min_y, max_y;

    if (type == AOI) {
        if (conf.area_op.size() < 3) {
            cerr << "[bin] ERROR: AOI requested but conf.area_op has < 3 points\n";
            return false;
        }
        min_x =  1e30; max_x = -1e30;
        min_y =  1e30; max_y = -1e30;
        for (auto& [lat, lon] : conf.area_op) {
            double px, py;
            if (!latlon_to_projected(lat, lon, px, py)) {
                cerr << "[bin] ERROR: could not project AOI vertex\n";
                return false;
            }
            min_x = min(min_x, px); max_x = max(max_x, px);
            min_y = min(min_y, py); max_y = max(max_y, py);
        }

        // Optional margin, in km, passed via range_km. Assets are PLACED
        // inside the AOI but SEE outward from it — a 75km radar sitting on
        // the AOI boundary needs terrain 75km beyond that boundary to
        // compute its coverage. Loading only the polygon's bounding box
        // leaves no data out there, and the coverage code then silently
        // truncates the radius to the window edge, which looks like radars
        // of different sizes even when they all share one range.
        // Pass the largest asset range as range_km to avoid that.
        if (range_km > 0.0) {
            const double margin = range_km * 1000.0
                                + WINDOW_MARGIN_CELLS * (double)resol;
            min_x -= margin; max_x += margin;
            min_y -= margin; max_y += margin;
            cout << "[bin] window AOI — " << conf.area_op.size()
                 << " polygon vertices, +" << range_km << " km margin\n";
        } else {
            cout << "[bin] window AOI — " << conf.area_op.size()
                 << " polygon vertices (no margin — asset coverage will be "
                    "truncated at the AOI edge)\n";
        }

    } else { // UNBOUNDED
        if (range_km <= 0.0) {
            cerr << "[bin] ERROR: UNBOUNDED requested but range_km <= 0\n";
            return false;
        }
        double cx, cy;
        if (!latlon_to_projected(center_lat, center_lng, cx, cy)) {
            cerr << "[bin] ERROR: could not project center ("
                 << center_lat << ", " << center_lng << ")\n";
            return false;
        }
        // Pad by WINDOW_MARGIN_CELLS. Engines clamp their radius to the
        // distance from the source to the nearest window edge
        // (maxPossible = min(centerX, cols - centerX - 1, ...)), so a
        // window sized to exactly +/- range_km trips a spurious one-cell
        // clamp. The margin keeps the full requested radius usable.
        const double half = range_km * 1000.0 + WINDOW_MARGIN_CELLS * (double)resol;
        min_x = cx - half; max_x = cx + half;
        min_y = cy - half; max_y = cy + half;
        cout << "[bin] window UNBOUNDED — center (" << center_lat << ", " << center_lng
             << ") +/- " << range_km << " km\n";
    }

    // projected → grid. Note y is flipped: max_y is the TOP (row 0 side).
    long long c_lo = (long long)floor((min_x - ox) / resol);
    long long c_hi = (long long)ceil ((max_x - ox) / resol);
    long long r_lo = (long long)floor((oy - max_y) / resol);
    long long r_hi = (long long)ceil ((oy - min_y) / resol);

    // clamp to the grid; a window that misses entirely is an error
    c_lo = max(0LL, min(c_lo, (long long)grid_cols));
    c_hi = max(0LL, min(c_hi, (long long)grid_cols));
    r_lo = max(0LL, min(r_lo, (long long)grid_rows));
    r_hi = max(0LL, min(r_hi, (long long)grid_rows));

    if (c_hi <= c_lo || r_hi <= r_lo) {
        cerr << "[bin] ERROR: requested region does not overlap the grid\n";
        return false;
    }

    out.r0   = (int)r_lo;
    out.c0   = (int)c_lo;
    out.rows = (int)(r_hi - r_lo);
    out.cols = (int)(c_hi - c_lo);

    cout << "[bin] window rows [" << out.r0 << " → " << (out.r0 + out.rows)
         << "] cols [" << out.c0 << " → " << (out.c0 + out.cols) << "]"
         << "  = " << out.cols << " x " << out.rows
         << " (" << (out.count() * 100.0 / ((double)grid_rows * grid_cols))
         << "% of grid)\n";
    return true;
}

// ── elevation / slope, windowed ────────────────────────────────────────────
// Prefers the tile folder; falls back to a seek-per-row read of the flat
// .bin so windowing works even before a tiled rebuild.
static void load_elev(Bin* b, const string& elev_path, const GridWindow& win) {
    if (b->elevation_loaded) return;
    if (b->GRID_ROWS <= 0 || b->GRID_COLS <= 0) {
        cerr << "[bin] ERROR: load_elev called before the grid dimensions were "
                "resolved - go through Bin::init()/Global::load()\n";
        return;
    }
    if (!Config::is_available()) {
        cerr << "[bin] ERROR: load_elev requires a config (grid.block_size) — none available.\n";
        return;
    }
    const int block_size = Config::load().grid.block_size;
    b->elevation_loaded = load_matrix<int16_t>(b, elev_path, "elev", block_size, win, b->elev_matrix);
}

static void load_slope(Bin* b, const string& slope_path, const GridWindow& win) {
    if (b->slope_loaded) return;
    if (b->GRID_ROWS <= 0 || b->GRID_COLS <= 0) {
        cerr << "[bin] ERROR: load_slope called before the grid dimensions were "
                "resolved - go through Bin::init()/Global::load()\n";
        return;
    }
    if (!Config::is_available()) {
        cerr << "[bin] ERROR: load_slope requires a config (grid.block_size) — none available.\n";
        return;
    }
    const int block_size = Config::load().grid.block_size;
    // Slope was previously optional-on-failure ("slope disabled") rather than
    // a hard error, since a lot of GA scoring runs fine without it. Keep that:
    // log through the same load_matrix() path, but don't propagate failure
    // as an error the way elev/vec do.
    b->slope_loaded = load_matrix<float>(b, slope_path, "slope", block_size, win, b->slope_matrix);
    if (!b->slope_loaded)
        cerr << "[bin] WARN: slope unavailable for this window - slope constraint disabled\n";
}

// ── unload / reload ────────────────────────────────────────────────────────
// Frees every windowed matrix and clears the window so a subsequent init()
// resolves a fresh one.
//
// CHANGED: vec_matrix used to be kept here on the grounds that it was
// full-grid and window-independent. It is WINDOWED now (see load_vec), so a
// vec_matrix kept across a window change holds the OLD region while
// op_window describes the NEW one — load_vec's early-out on vector_loaded
// would then skip the reload and getConstraintRC would index the old buffer
// with new-window offsets (wrong cells, or past the end). It is therefore
// always dropped. `drop_vec` is accepted only so existing callers compile.
void Bin::unload(bool /*drop_vec*/) {
    elev_matrix.clear();  elev_matrix.shrink_to_fit();
    slope_matrix.clear(); slope_matrix.shrink_to_fit();
    vec_matrix.clear();   vec_matrix.shrink_to_fit();
    elevation_loaded = false;
    slope_loaded     = false;
    vector_loaded    = false;

    window_valid = false;
    op_window    = GridWindow{};
    data_window  = GridWindow{};
    cout << "[bin] unloaded (elev, slope, vec)\n";
}

// Move to a different region without tearing down the Bin. Any raw index
// or (r,c) captured against the previous window is invalid after this —
// engines holding gridCols/gridRows or a translated centre must be rebuilt.
void Bin::reload(set<LoadOptions> options,
                 OperationType type,
                 double center_lat, double center_lng, double range_km) {
    cout << "[bin] reload → type=" << (int)type;
    if (type == UNBOUNDED)
        cout << " center=(" << center_lat << ", " << center_lng
             << ") range=" << range_km << "km";
    cout << "\n";

    // vec is windowed, so if it was resident it has to follow the window to
    // the new region — otherwise reload({ELEVATION}, ...) would leave the
    // caller with elevation for the new area and no constraints at all.
    const bool had_vec = vector_loaded;
    unload(/*drop_vec=*/true);
    if (had_vec && !options.empty() && !options.count(ALL) && !options.count(CONSTRAINTS))
        options.insert(CONSTRAINTS);

    init(options, type, center_lat, center_lng, range_km);
}

// Shared by init() (which may resolve a NEW window first) and
// ensureLoaded() (which never touches the window) — both end up loading
// exactly the same way once a window is known.
static void dispatchLoad(Bin* b, set<LoadOptions> options,
                         const string& param_path, const string& elev_path,
                         const string& slope_path) {
    // Loads into data_window — the RESIDENT layout — not op_window. In tiled
    // mode they're identical. In load_full mode data_window is the whole
    // grid, so the very same loaders read the flat bins in full: one code
    // path for both layouts, and slope now loads in flat mode too (the old
    // load_full branch only ever loaded vec + elev).
    const GridWindow& w = b->data_window;
    if (options.empty()) {
        load_vec (b, param_path, w);
        load_elev(b, elev_path, w);

    } else if (options.count(ALL)) {
        load_vec  (b, param_path, w);
        load_elev (b, elev_path,  w);
        load_slope(b, slope_path, w);

    } else {
        if (options.count(CONSTRAINTS)) load_vec  (b, param_path, w);
        if (options.count(ELEVATION))   load_elev (b, elev_path,  w);
        if (options.count(SLOPE_))      load_slope(b, slope_path, w);
    }
}

// Both callers (ensureLoaded, init) only ever reach this after window_valid
// is true, which itself now requires Config::is_available() (computeWindow
// hard-errors otherwise) — so this can't actually observe a missing config
// in practice. Errors instead of falling back to DEFAULT_PARAM_MATRIX/etc
// anyway: those constants describe a 100m grid, and silently pointing a
// 30m run at them is exactly the class of bug this whole pass is removing.
static bool resolvePaths(string& param_path, string& elev_path, string& slope_path) {
    if (!Config::is_available()) {
        cerr << "[bin] ERROR: resolvePaths requires a config (output.matrix_file/"
                "elevation_file/slope_file) — none available.\n";
        return false;
    }
    Config cfg  = Config::load();
    param_path  = cfg.output.matrix_file;
    elev_path   = cfg.output.elevation_file;
    slope_path  = cfg.output.slope_file;
    return true;
}

// Load under whatever window is ALREADY resolved — never resolves or
// changes the window itself. Use this for "make sure X is loaded" checks
// (as the per-call getters below do); use init()/Global::load() to pick
// or change the window. Errors out if no window has been resolved yet,
// rather than silently defaulting to FULL — a caller that hasn't called
// load() yet has a bug worth surfacing, not papering over.
void Bin::ensureLoaded(set<LoadOptions> options) {
    if (!window_valid) {
        cerr << "[bin] ERROR: ensureLoaded() called with no window resolved — "
                "call Global::load() (or Bin::init()) with an OperationType "
                "first\n";
        return;
    }
    string param_path, elev_path, slope_path;
    if (!resolvePaths(param_path, elev_path, slope_path)) return;
    dispatchLoad(this, options, param_path, elev_path, slope_path);
}

void Bin::init(set<LoadOptions> options,
               OperationType type,
               double center_lat, double center_lng, double range_km) {
    if (!Config::is_available()) {
        cerr << "[bin] config not found, cant load anything." << endl;
        return;
    }

    // resolve paths once
    string param_path, elev_path, slope_path;
    if (!resolvePaths(param_path, elev_path, slope_path)) return;

    // Decide the window before touching any file. If a window was already
    // resolved by a previous init() on this Bin, keep it — re-resolving
    // under a different type would silently invalidate indices into the
    // matrices already loaded.
    if (!window_valid) {
        if (!computeWindow(type, center_lat, center_lng, range_km, op_window)) {
            cerr << "[bin] ERROR: could not resolve window — nothing loaded\n";
            return;
        }
        op_type      = type;
        window_valid = true;

        // What actually becomes resident. op_window is the region of
        // interest; data_window is the buffer layout every index uses.
        //   tiled     : load just the region  -> data_window = op_window
        //   load_full : flat bins load whole  -> data_window = full grid
        // With data_window = full grid, index(r,c) == r*GRID_COLS + c, i.e.
        // the legacy layout — so no caller needs its own load_full branch.
        if (load_full) data_window = GridWindow{ 0, 0, GRID_ROWS, GRID_COLS };
        else           data_window = op_window;

        cout << "[bin] resident layout : "
             << (load_full ? "FULL grid (block_size=0, flat bins)" : "window (tiled)")
             << "  " << data_window.cols << " x " << data_window.rows
             << "  @ " << cell_m << " m/cell\n";
    } else if (type != op_type) {
        cerr << "[bin] WARN: init() called with type " << (int)type
             << " but this Bin already holds a " << (int)op_type
             << " window — keeping the existing window\n";
    }

    dispatchLoad(this, options, param_path, elev_path, slope_path);
}
