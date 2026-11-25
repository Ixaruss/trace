#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <vector>
#include <algorithm>
#include <string>
#include <png.h>
#include "bin.h"
#include "lib.h"
#include "tile.h"
#include "../common.h"
#include <chrono>
#include <limits>
#include <stdio.h>
using namespace std;

// ══════════════════════════════════════════════════════════════════════════
// STREAMING PREPROCESS
//
// Nothing grid-sized is ever held in RAM. Previously this file allocated:
//   matrix          calloc(N,1)          N bytes, held for the whole run
//   rasterize_layer GDAL MEM band + buf  2N bytes, per job
//   expand_buffer   distance array       N..4N bytes, per job
// At 30m over a subcontinental extent N is ~2.06e10, so `matrix` alone is
// ~20.6GB and the jobs loop peaked near 60GB.
//
// Now:
//   - rasterization targets a disk-backed GTiff; GDAL streams it through
//     its own bounded block cache instead of materialising the grid
//   - the distance transform runs as two streaming row passes IN PLACE on
//     that raster (forward top->bottom, backward bottom->top), needing
//     three rows of RAM rather than a full distance array
//   - the bit matrix accumulates into a second disk raster, OR'd a row at
//     a time, then is written out as tiles
//
// Peak RAM is now GDAL's cache (bounded below) plus a handful of rows.
// Peak DISK is the matrix raster plus one job raster at a time.
// ══════════════════════════════════════════════════════════════════════════

// Bound GDAL's block cache explicitly — it defaults to a share of physical
// RAM, which on a big machine would quietly undo the point of streaming.
static const GIntBig GDAL_CACHE_BYTES = 512LL * 1024 * 1024;

static bool make_dir(const string& path) {
#ifdef _WIN32
    int ret = system(("mkdir \"" + path + "\" 2>nul").c_str());
#else
    int ret = system(("mkdir -p \"" + path + "\"").c_str());
#endif
    (void)ret;
    return true;
}

// ── disk-backed uint8 grid raster ─────────────────────────────────────────
// Striped with one row per block: access here is strictly row-wise (both
// distance-transform passes and every combine pass), so one-row strips keep
// GDAL's cache working set minimal and avoid reading/rewriting whole tiles
// for a single row. Uncompressed because the DT rewrites each row after
// reading it, and rewriting compressed strips in place is a good way to
// corrupt or balloon the file. SPARSE_OK means untouched blocks cost no
// disk and read back as zero.
static GDALDataset* create_u8_raster(const string& path,
                                     int ROWS, int COLS,
                                     double ox, double oy, int cell,
                                     GDALDataType dt = GDT_Byte) {
    GDALDriver* drv = GetGDALDriverManager()->GetDriverByName("GTiff");
    if (!drv) { cerr << "[preprocess] ERROR: GTiff driver unavailable\n"; return nullptr; }

    char** opts = nullptr;
    opts = CSLSetNameValue(opts, "BIGTIFF",   "YES");
    opts = CSLSetNameValue(opts, "TILED",     "NO");
    opts = CSLSetNameValue(opts, "BLOCKYSIZE","1");
    opts = CSLSetNameValue(opts, "COMPRESS",  "NONE");
    opts = CSLSetNameValue(opts, "SPARSE_OK", "TRUE");

    GDALDataset* ds = drv->Create(path.c_str(), COLS, ROWS, 1, dt, opts);
    CSLDestroy(opts);
    if (!ds) {
        cerr << "[preprocess] ERROR: could not create raster " << path << "\n";
        return nullptr;
    }
    double gt[6] = { ox, (double)cell, 0, oy, 0, -(double)cell };
    ds->SetGeoTransform(gt);
    return ds;
}

static inline bool row_read(GDALDataset* ds, int r, uint8_t* buf, int COLS) {
    return ds->GetRasterBand(1)->RasterIO(GF_Read, 0, r, COLS, 1,
                                          buf, COLS, 1, GDT_Byte, 0, 0) == CE_None;
}
static inline bool row_write(GDALDataset* ds, int r, const uint8_t* buf, int COLS) {
    return ds->GetRasterBand(1)->RasterIO(GF_Write, 0, r, COLS, 1,
                                          const_cast<uint8_t*>(buf), COLS, 1,
                                          GDT_Byte, 0, 0) == CE_None;
}

// ── streaming PNG writer ─────────────────────────────────────────────────
struct PngRowWriter {
    FILE*       f    = nullptr;
    png_structp png  = nullptr;
    png_infop   info = nullptr;
    bool        ok   = false;

    bool open(const string& path, int rows, int cols) {
        f = fopen(path.c_str(), "wb");
        if (!f) { cerr << "[preprocess] ERROR: cannot open " << path << " for write\n"; return false; }

        png = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
        if (!png) { cerr << "[preprocess] ERROR: png_create_write_struct failed\n"; fclose(f); f = nullptr; return false; }

        info = png_create_info_struct(png);
        if (!info) {
            cerr << "[preprocess] ERROR: png_create_info_struct failed\n";
            png_destroy_write_struct(&png, nullptr); fclose(f); f = nullptr; return false;
        }
        if (setjmp(png_jmpbuf(png))) {
            cerr << "[preprocess] ERROR: libpng error opening " << path << "\n";
            png_destroy_write_struct(&png, &info); fclose(f); f = nullptr; return false;
        }

        png_init_io(png, f);
        png_set_IHDR(png, info,
                     (png_uint_32)cols, (png_uint_32)rows,
                     8, PNG_COLOR_TYPE_RGB,
                     PNG_INTERLACE_NONE,
                     PNG_COMPRESSION_TYPE_DEFAULT,
                     PNG_FILTER_TYPE_DEFAULT);
        png_write_info(png, info);
        ok = true;
        return true;
    }

    bool write_row(const uint8_t* rgb_row) {
        if (!ok) return false;
        if (setjmp(png_jmpbuf(png))) {
            cerr << "[preprocess] ERROR: libpng error writing row\n";
            ok = false;
            return false;
        }
        png_write_row(png, const_cast<png_bytep>(rgb_row));
        return true;
    }

    void close() {
        if (!ok) { if (f) fclose(f); f = nullptr; return; }
        png_write_end(png, nullptr);
        png_destroy_write_struct(&png, &info);
        fclose(f);
        f  = nullptr;
        ok = false;
    }
};

// Stream a 0/1 mask raster straight to a PNG, one row at a time.
static bool mask_raster_to_png(GDALDataset* ds, const string& path,
                               int ROWS, int COLS) {
    PngRowWriter pw;
    if (!pw.open(path, ROWS, COLS)) return false;

    vector<uint8_t> row((size_t)COLS);
    vector<uint8_t> rgb((size_t)COLS * 3);
    for (int r = 0; r < ROWS; r++) {
        if (!row_read(ds, r, row.data(), COLS)) break;
        for (int c = 0; c < COLS; c++) {
            uint8_t v = row[c] ? 255 : 0;
            rgb[c*3] = v; rgb[c*3+1] = v; rgb[c*3+2] = v;
        }
        if (!pw.write_row(rgb.data())) break;
        if (r % 5000 == 0 && r > 0)
            cout << "[preprocess] mask png row " << r << "/" << ROWS << "\n";
    }
    pw.close();
    cout << "[preprocess] png → " << path << "\n";
    return true;
}

// Stream a mask raster to a flat .bin (unchanged on-disk format).
static bool mask_raster_to_bin(GDALDataset* ds, const string& path,
                               int ROWS, int COLS) {
    FILE* f = fopen(path.c_str(), "wb");
    if (!f) { cerr << "[preprocess] ERROR: cannot open " << path << "\n"; return false; }
    vector<uint8_t> row((size_t)COLS);
    for (int r = 0; r < ROWS; r++) {
        if (!row_read(ds, r, row.data(), COLS)) { fclose(f); return false; }
        if (fwrite(row.data(), 1, (size_t)COLS, f) != (size_t)COLS) {
            cerr << "[preprocess] ERROR: short write on " << path << "\n";
            fclose(f); return false;
        }
    }
    fclose(f);
    cout << "[preprocess] bin → " << path
         << " (" << ((double)ROWS * COLS) / 1e6 << " MB)\n";
    return true;
}

// ── rasterize one layer into a disk-backed raster ─────────────────────────
static bool rasterize_layer_to_raster(const string& shp_path,
                                       int ROWS, int COLS,
                                       double ox, double oy, int cell_size,
                                       bool boundary_only,
                                       const string& out_path,
                                       size_t& burned_out,
                                       GDALDataType dt = GDT_Byte) {
    burned_out = 0;

    GDALDataset* shp_ds = (GDALDataset*)GDALOpenEx(
        shp_path.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr);
    if (!shp_ds) { cerr << "[preprocess] ERROR: cannot open " << shp_path << "\n"; return false; }

    OGRLayer* layer = shp_ds->GetLayer(0);
    if (!layer) {
        cerr << "[preprocess] ERROR: no layer in " << shp_path << "\n";
        GDALClose(shp_ds); return false;
    }

    OGREnvelope env;
    if (layer->GetExtent(&env) == OGRERR_NONE) {
        cout << "[preprocess] shp extent: ("
             << env.MinX << ", " << env.MinY << ") - ("
             << env.MaxX << ", " << env.MaxY << ")\n";
    }

    vector<OGRGeometryH> geoms;
    OGRFeature* feat;
    layer->ResetReading();
    while ((feat = layer->GetNextFeature()) != nullptr) {
        OGRGeometry* geom = feat->GetGeometryRef();
        if (geom) {
            if (boundary_only) {
                OGRGeometry* bnd = geom->Boundary();
                if (bnd) geoms.push_back((OGRGeometryH)bnd);
            } else {
                geoms.push_back((OGRGeometryH)geom->clone());
            }
        }
        OGRFeature::DestroyFeature(feat);
    }
    GDALClose(shp_ds);

    if (geoms.empty()) {
        cerr << "[preprocess] WARN: no geometries in " << shp_path << "\n";
        return false;
    }
    cout << "[preprocess] loaded " << geoms.size() << " geoms\n";
    cout << "[preprocess] raster: " << COLS << " x " << ROWS
         << " origin=(" << ox << ", " << oy << ")  → " << out_path << "\n";

    // CHANGED: was a GDAL "MEM" dataset (a full COLS*ROWS byte buffer in
    // RAM) plus a second same-size caller buffer that RasterIO copied into.
    // Rasterizing into a disk-backed GTiff lets GDAL page blocks through
    // its bounded cache instead, so nothing grid-sized is resident.
    GDALDataset* ds = create_u8_raster(out_path, ROWS, COLS, ox, oy, cell_size, dt);
    if (!ds) {
        for (auto g : geoms) OGR_G_DestroyGeometry(g);
        return false;
    }

    int band_list[1] = { 1 };
    vector<double> burn_vals(geoms.size(), 1.0);
    CPLErr err = (CPLErr)GDALRasterizeGeometries(
        ds, 1, band_list,
        (int)geoms.size(), geoms.data(),
        nullptr, nullptr, burn_vals.data(),
        nullptr, nullptr, nullptr);

    for (auto g : geoms) OGR_G_DestroyGeometry(g);

    if (err != CE_None) {
        cerr << "[preprocess] ERROR: GDALRasterizeGeometries failed (CPLErr=" << err << ")\n";
        GDALClose(ds);
        return false;
    }

    // count burned cells by streaming, not by holding the grid
    {
        vector<uint8_t> row((size_t)COLS);
        for (int r = 0; r < ROWS; r++) {
            if (!row_read(ds, r, row.data(), COLS)) break;
            for (int c = 0; c < COLS; c++) if (row[c]) burned_out++;
        }
    }
    cout << "[preprocess] burned cells after rasterize: " << burned_out << "\n";

    // MUST close before any reopen — GTiff buffers writes and only
    // guarantees a flush on GDALClose.
    GDALClose(ds);
    return true;
}

// ── streaming distance transform, in place on the raster ──────────────────
// Two-pass city-block chamfer, exactly the algorithm the in-RAM version
// used, with distances saturated at radius+1 so they fit in the same uint8
// raster the mask lives in. Saturation cannot change the thresholded
// result: any cell whose true distance d <= radius has every value along
// its optimal path <= d <= radius < INF, so no clamping occurs on that
// path; cells with d > radius end up >= radius+1 either way.
//
// Forward pass reads row r (raw 0/1 on the first pass), computes the
// running distance using the previous row, and writes it back to row r —
// the raw value is fully consumed by then, so one raster serves as both
// input and output. Backward pass does the same in reverse and thresholds
// to 0/1 on the way out. RAM: three rows.
template <typename T> static inline bool row_read_t(GDALDataset* ds, int r, T* buf, int COLS, GDALDataType dt) {
    return ds->GetRasterBand(1)->RasterIO(GF_Read, 0, r, COLS, 1, buf, COLS, 1, dt, 0, 0) == CE_None;
}
template <typename T> static inline bool row_write_t(GDALDataset* ds, int r, const T* buf, int COLS, GDALDataType dt) {
    return ds->GetRasterBand(1)->RasterIO(GF_Write, 0, r, COLS, 1, const_cast<T*>(buf), COLS, 1, dt, 0, 0) == CE_None;
}

// The element type T is the raster's own type: uint8_t (max radius 254) or
// uint16_t (max radius 65534). Wider radii need the wider raster because the
// saturated distance INF = radius + 1 has to fit in it. The 16-bit raster is
// only created for the layers that need it (see the jobs loop) — at 30 m that
// is the 10 km IB buffer (333 cells), which used to be SKIPPED with an error,
// silently leaving IB out of the constraints map.
template <typename T>
static bool expand_raster_t(const string& path, int ROWS, int COLS,
                            int radius, size_t& expanded_out, GDALDataType dt) {
    expanded_out = 0;

    GDALDataset* ds = (GDALDataset*)GDALOpenEx(
        path.c_str(), GDAL_OF_RASTER | GDAL_OF_UPDATE, nullptr, nullptr, nullptr);
    if (!ds) { cerr << "[preprocess] ERROR: cannot reopen " << path << " for update\n"; return false; }

    const int MAXR = (int)std::numeric_limits<T>::max() - 1;
    if (radius > MAXR) {
        cerr << "[preprocess] ERROR: buffer radius " << radius << " cells exceeds the " << MAXR
             << " representable in this raster type.\n";
        GDALClose(ds);
        return false;
    }

    const int INF = radius + 1;
    vector<T> prev((size_t)COLS), cur((size_t)COLS), nxt((size_t)COLS);

    if (radius <= 0) {
        cout << "[preprocess] radius=0 — skipping expansion\n";
        for (int r = 0; r < ROWS; r++) {
            if (!row_read_t<T>(ds, r, cur.data(), COLS, dt)) break;
            for (int c = 0; c < COLS; c++) if (cur[c]) expanded_out++;
        }
        GDALClose(ds);
        return true;
    }

    // forward: top → bottom, left → right
    for (int r = 0; r < ROWS; r++) {
        if (!row_read_t<T>(ds, r, cur.data(), COLS, dt)) { GDALClose(ds); return false; }
        for (int c = 0; c < COLS; c++) {
            int v = cur[c] ? 0 : INF;
            if (c > 0) { int t = (int)cur[c-1] + 1; if (t < v) v = t; }
            if (r > 0) { int t = (int)prev[c]  + 1; if (t < v) v = t; }
            if (v > INF) v = INF;
            cur[c] = (T)v;
        }
        if (!row_write_t<T>(ds, r, cur.data(), COLS, dt)) { GDALClose(ds); return false; }
        prev.swap(cur);
        if (r % 20000 == 0 && r > 0)
            cout << "[preprocess] dt forward row " << r << "/" << ROWS << "\n";
    }

    // backward: bottom → top, right → left, thresholding on the way out
    bool have_next = false;
    for (int r = ROWS - 1; r >= 0; r--) {
        if (!row_read_t<T>(ds, r, cur.data(), COLS, dt)) { GDALClose(ds); return false; }
        for (int c = COLS - 1; c >= 0; c--) {
            int v = cur[c];
            if (c < COLS-1) { int t = (int)cur[c+1] + 1; if (t < v) v = t; }
            if (have_next)  { int t = (int)nxt[c]   + 1; if (t < v) v = t; }
            if (v > INF) v = INF;
            cur[c] = (T)v;
        }
        nxt.assign(cur.begin(), cur.end());   // distances seen by the row above
        have_next = true;

        for (int c = 0; c < COLS; c++) {
            cur[c] = (cur[c] <= radius) ? 1 : 0;
            if (cur[c]) expanded_out++;
        }
        if (!row_write_t<T>(ds, r, cur.data(), COLS, dt)) { GDALClose(ds); return false; }
        if (r % 20000 == 0 && r < ROWS-1)
            cout << "[preprocess] dt backward row " << r << "/" << ROWS << "\n";
    }

    cout << "[preprocess] expanded → " << expanded_out << " cells\n";
    GDALClose(ds);
    return true;
}

// Which raster element type a given buffer radius (in cells) needs.
static inline GDALDataType raster_type_for_radius(int radius) {
    return radius > 254 ? GDT_UInt16 : GDT_Byte;
}

static bool expand_raster(const string& path, int ROWS, int COLS,
                          int radius, size_t& expanded_out) {
    return radius > 254
        ? expand_raster_t<uint16_t>(path, ROWS, COLS, radius, expanded_out, GDT_UInt16)
        : expand_raster_t<uint8_t >(path, ROWS, COLS, radius, expanded_out, GDT_Byte);
}

// OR a job's 0/1 mask into the accumulating bit matrix, a row at a time.
static bool or_mask_into_matrix(GDALDataset* mat, GDALDataset* job,
                                int ROWS, int COLS, int bit,
                                long long& total, long long& hits) {
    vector<uint8_t> mrow((size_t)COLS), jrow((size_t)COLS);
    for (int r = 0; r < ROWS; r++) {
        if (!row_read(mat, r, mrow.data(), COLS)) return false;
        if (!row_read(job, r, jrow.data(), COLS)) return false;
        for (int c = 0; c < COLS; c++) {
            total++;
            if (jrow[c]) { mrow[c] |= (uint8_t)(1 << bit); hits++; }
        }
        if (!row_write(mat, r, mrow.data(), COLS)) return false;
    }
    return true;
}

void Global::preprocess(bool images, bool timer, bool seperate, int individual) {
    GDALAllRegister();
    GDALSetCacheMax64(GDAL_CACHE_BYTES);

    chrono::time_point<chrono::steady_clock, chrono::duration<long, ratio<1, 1000000000>>> start,end;
    if (timer){ cout<< ">>> Starting timer count" << endl;
         start = chrono::steady_clock::now();
    }
    if (!Config::is_available()) {
        cerr << "[preprocess] ERROR: no config available — road/rail/water/ib "
                "shapefile paths, distances, grid extent and block_size all "
                "come from it. Refusing to fall back to compiled-in DEFAULT_* "
                "paths, which describe a different (100m) grid.\n";
        return;
    }
    string mpath = Config::load().output.elevation_file;

    if (images) make_dir("./output/preprocess");
    if(!filesystem::exists(mpath) && !filesystem::exists(tilegrid::tile_dir_for(mpath)))
        if ( individual == -1 || individual == 6)  Bin::build_elevation_matrix();

    string road_shp, rail_shp, warea_shp, wline_shp, ib_shp, elev_bin_path, opath, slope_path;
    int road_bit, rail_bit, waterarea_bit, waterlines_bit, ib_bit, slope_bit;
    int dist_road, dist_rail, dist_water, dist_ib;
    int ELEV_R, ELEV_C;
    double ox, oy;
    double min_x = 0, max_x = 0, min_y = 0, max_y = 0;
    int block_size = 0;
    long long total = 0;
    long long feasible = 0;

    {
        Config cfg    = Config::load();
        road_shp      = cfg.input.road;
        rail_shp      = cfg.input.rail;
        warea_shp     = cfg.input.water_area;
        wline_shp     = cfg.input.water_lines;
        ib_shp        = cfg.input.ib;
        elev_bin_path = cfg.output.elevation_file;
        opath         = cfg.output.matrix_file;

        dist_road  = cfg.distance.from_road  > 0 ? cfg.distance.from_road  : 1000;
        dist_rail  = cfg.distance.from_rail  > 0 ? cfg.distance.from_rail  : 1000;
        dist_water = cfg.distance.from_water > 0 ? cfg.distance.from_water : 1000;
        dist_ib    = cfg.distance.from_ib    > 0 ? cfg.distance.from_ib    : 10000;

        slope_bit  = 5;
        slope_path = cfg.output.slope_file;

        ox = cfg.grid.min_x;
        oy = cfg.grid.max_y;
        min_x = cfg.grid.min_x; max_x = cfg.grid.max_x;
        min_y = cfg.grid.min_y; max_y = cfg.grid.max_y;
        ELEV_R = cfg.grid.rows;
        ELEV_C = cfg.grid.cols;
        block_size = cfg.grid.block_size;

        if (ELEV_R <= 0 || ELEV_C <= 0 || cfg.grid.resolution <= 0) {
            cerr << "[preprocess] ERROR: config grid.rows/grid.cols/grid.resolution "
                    "are not set — run build_elevation_matrix() first.\n";
            return;
        }
    }

    road_bit       = BIT_ROAD;
    rail_bit       = BIT_RAIL;
    waterarea_bit  = BIT_WATER_AREA;
    waterlines_bit = BIT_WATER_LINE;
    ib_bit         = BIT_IB;
    slope_bit      = 5;

    auto idx = opath.find_last_of("\\/");
    string path;
    if(std::string::npos != idx) path = opath.substr(0,idx+1);
    make_dir(path);

    // Grid dims come straight from config — cfg.grid.rows/cols ARE the
    // authoritative dims (checked non-zero above), not something to
    // recompute from extents/RESOLUTION as a fallback.
    const int RES  = lib::grid_resolution();
    const int ROWS = ELEV_R;
    const int COLS = ELEV_C;
    const size_t N = (size_t)ROWS * COLS;

    const bool tiled_out = (block_size > 0);

    cout << "[preprocess] resolution : " << RES << " m/cell"
         << (RES != (int)RESOLUTION ? "  (NOTE: compile-time RESOLUTION differs)" : "")
         << "\n";
    cout << "[preprocess] Grid       : " << ROWS << " x " << COLS << " (" << N / 1e6 << "M cells)\n";
    cout << "[preprocess] origin     : (" << ox << ", " << oy << ")\n";
    cout << "[preprocess] dist_road  : " << dist_road  << "m (" << dist_road  / RES << " cells)\n";
    cout << "[preprocess] dist_rail  : " << dist_rail  << "m (" << dist_rail  / RES << " cells)\n";
    cout << "[preprocess] dist_water : " << dist_water << "m (" << dist_water / RES << " cells)\n";
    cout << "[preprocess] dist_ib    : " << dist_ib    << "m (" << dist_ib    / RES << " cells)\n";
    cout << "[preprocess] elev grid  : " << ELEV_C << " x " << ELEV_R << "\n";
    cout << "[preprocess] output     : "
         << (tiled_out ? ("tiled, block=" + to_string(block_size)) : string("flat bin")) << "\n\n";

    // Disk pre-flight. Temps are the accumulating matrix raster plus one job
    // raster at a time; the final map is another N bytes in EITHER layout
    // (tiles hold the same cells as the flat file). Per-layer flat bins add N
    // each — and note `seperate || individual` below is true for the default
    // individual == -1, so a plain full run writes all five.
    {
        namespace fs = std::filesystem;
        error_code sec;
        const int per_layer_bins = (seperate || individual != 0) ? 5 : 0;
        const bool wide_job = (dist_road / RES > 254) || (dist_rail / RES > 254) ||
                              (dist_water / RES > 254) || (dist_ib / RES > 254);
        // matrix raster + final map + one job raster (2N if any layer needs 16-bit)
        const double need_gb = ((double)N * ((2 + (wide_job ? 2 : 1)) + per_layer_bins)) / 1e9;
        const auto sp = fs::space(path.empty() ? fs::path(".") : fs::path(path), sec);
        if (!sec) {
            const double have_gb = (double)sp.available / 1e9;
            cout << "[preprocess] disk       : need ~" << need_gb << " GB, "
                 << have_gb << " GB free\n";
            if (have_gb < need_gb)
                cerr << "[preprocess] WARN: free space is below the estimate — a full disk "
                        "mid-write leaves map.bin/map_tiles incomplete and GA will not load it.\n";
        }
    }

    // ── the accumulating bit matrix lives on disk ────────────────────────
    const string mat_tmp = path + "matrix_tmp.tif";
    GDALDataset* mat = create_u8_raster(mat_tmp, ROWS, COLS, ox, oy, RES);
    if (!mat) return;

    struct LayerJob {
        string name;
        string shp;
        int    distance;
        bool   boundary_only;
        int    index;
    };
    vector<LayerJob> jobs;
    if(individual != -1) {
        switch (individual) {
            case 1: jobs.push_back({ "road",      road_shp,  dist_road,  false, road_bit }); break;
            case 2: jobs.push_back({ "rail",      rail_shp,  dist_rail,  false, rail_bit }); break;
            case 3: jobs.push_back({ "waterarea", warea_shp, dist_water, false, waterarea_bit }); break;
            case 4: jobs.push_back({ "waterline", wline_shp, dist_water, false, waterlines_bit }); break;
            case 5: jobs.push_back({ "ib",        ib_shp,    dist_ib,    true,  ib_bit  }); break;
        }
    } else {
        jobs.push_back({ "road",      road_shp,  dist_road,  false, road_bit });
        jobs.push_back({ "rail",      rail_shp,  dist_rail,  false, rail_bit });
        jobs.push_back({ "waterarea", warea_shp, dist_water, false, waterarea_bit });
        jobs.push_back({ "waterline", wline_shp, dist_water, false, waterlines_bit });
        jobs.push_back({ "ib",        ib_shp,    dist_ib,    true,  ib_bit  });
    }

    int elevTime = 0;

    if(individual == 6 || individual == -1) {

    // ── elevation + slope ────────────────────────────────────────────────
    cout << "\n=== elevation + slope ===\n";

    size_t elev_n = (size_t)ELEV_R * ELEV_C;
    bool want_png = (images || individual == 6);
    static const size_t TARGET_CHUNK_BYTES = 64ull * 1024 * 1024;

    tilegrid::Grid elev_grid, slope_grid;
    const bool elev_tiled  = elev_grid.discover(tilegrid::tile_dir_for(elev_bin_path));
    const bool slope_tiled = slope_grid.discover(tilegrid::tile_dir_for(slope_path));

    // elevation: PNG only, contributes no bits
    if (want_png) {
        if (elev_tiled) {
            elev_grid.print("elev");
            PngRowWriter pw;
            if (pw.open("./output/preprocess/elev.png", ELEV_R, ELEV_C)) {
                vector<int16_t> row(ELEV_C);
                vector<uint8_t> rgb((size_t)ELEV_C * 3);
                tilegrid::RowReader rd;
                for (int tr = 0; tr < elev_grid.tile_rows; tr++) {
                    if (!rd.open(elev_grid, tr)) break;
                    for (int rr = 0; rr < rd.height; rr++) {
                        if (!rd.read_row(row.data(), sizeof(int16_t))) break;
                        for (int c = 0; c < ELEV_C; c++)
                            Bin::get_elevation_color(row[c], rgb[c*3], rgb[c*3+1], rgb[c*3+2]);
                        pw.write_row(rgb.data());
                    }
                    rd.close();
                }
                pw.close();
                cout << "[preprocess] png → ./output/preprocess/elev.png\n";
            }
        } else {
            FILE* ef = fopen(elev_bin_path.c_str(), "rb");
            if (!ef) cerr << "[preprocess] WARN: elevation binary missing — skipping\n";
            else {
                PngRowWriter pw;
                if (pw.open("./output/preprocess/elev.png", ELEV_R, ELEV_C)) {
                    int chunk_rows = (int)max<size_t>(1, TARGET_CHUNK_BYTES / ((size_t)ELEV_C * sizeof(int16_t)));
                    if (chunk_rows > ELEV_R) chunk_rows = ELEV_R;
                    vector<int16_t> ebuf((size_t)chunk_rows * ELEV_C);
                    vector<uint8_t> rgb((size_t)ELEV_C * 3);
                    for (int r0 = 0; r0 < ELEV_R; r0 += chunk_rows) {
                        int rows_this = min(chunk_rows, ELEV_R - r0);
                        size_t n = (size_t)rows_this * ELEV_C;
                        if (fread(ebuf.data(), sizeof(int16_t), n, ef) != n) break;
                        for (int rr = 0; rr < rows_this; rr++) {
                            const int16_t* row = ebuf.data() + (size_t)rr * ELEV_C;
                            for (int c = 0; c < ELEV_C; c++)
                                Bin::get_elevation_color(row[c], rgb[c*3], rgb[c*3+1], rgb[c*3+2]);
                            pw.write_row(rgb.data());
                        }
                    }
                    pw.close();
                    cout << "[preprocess] png → ./output/preprocess/elev.png\n";
                }
                fclose(ef);
            }
        }
    }

    // slope: PNG + burn the slope bit straight into the matrix raster
    {
        PngRowWriter pw;
        const bool png_ok = want_png && pw.open("./output/preprocess/slope.png", ELEV_R, ELEV_C);
        vector<uint8_t> rgb((size_t)ELEV_C * 3);
        vector<uint8_t> mrow((size_t)COLS);
        vector<float>   srow(ELEV_C);
        bool any = false;

        auto burn_slope_row = [&](int r) {
            if (r >= ROWS) return;
            if (!row_read(mat, r, mrow.data(), COLS)) return;
            const int lim = min(COLS, ELEV_C);
            for (int c = 0; c < lim; c++)
                if (srow[c] > 45.0f) mrow[c] |= (uint8_t)(1 << slope_bit);
            row_write(mat, r, mrow.data(), COLS);
        };

        if (slope_tiled) {
            slope_grid.print("slope");
            any = true;
            tilegrid::RowReader rd;
            int gr = 0;
            for (int tr = 0; tr < slope_grid.tile_rows; tr++) {
                if (!rd.open(slope_grid, tr)) break;
                for (int rr = 0; rr < rd.height; rr++, gr++) {
                    if (!rd.read_row(srow.data(), sizeof(float))) break;
                    if (png_ok) {
                        for (int c = 0; c < ELEV_C; c++) {
                            float s = srow[c];
                            if (s < 0.0f) s = 0.0f;
                            if (s > 90.0f) s = 90.0f;
                            uint8_t v = (uint8_t)(s / 90.0f * 255.0f);
                            rgb[c*3] = v; rgb[c*3+1] = v; rgb[c*3+2] = v;
                        }
                        pw.write_row(rgb.data());
                    }
                    burn_slope_row(gr);
                }
                rd.close();
                cout << "[preprocess] slope tile row " << (tr+1) << "/" << slope_grid.tile_rows << "\n";
            }
        } else {
            FILE* sf = fopen(slope_path.c_str(), "rb");
            if (!sf) cerr << "[preprocess] WARN: slope binary missing — skipping\n";
            else {
                size_t expected = elev_n * sizeof(float);
                fseek64(sf, 0, SEEK_END);
                size_t fb = ftell64(sf);
                fseek64(sf, 0, SEEK_SET);
                if (fb != expected) {
                    cerr << "[preprocess] ERROR: slope size mismatch — expected "
                         << expected/1e6 << " MB, got " << fb/1e6 << " MB\n";
                } else {
                    any = true;
                    for (int r = 0; r < ELEV_R; r++) {
                        if (fread(srow.data(), sizeof(float), ELEV_C, sf) != (size_t)ELEV_C) break;
                        if (png_ok) {
                            for (int c = 0; c < ELEV_C; c++) {
                                float s = srow[c];
                                if (s < 0.0f) s = 0.0f;
                                if (s > 90.0f) s = 90.0f;
                                uint8_t v = (uint8_t)(s / 90.0f * 255.0f);
                                rgb[c*3] = v; rgb[c*3+1] = v; rgb[c*3+2] = v;
                            }
                            pw.write_row(rgb.data());
                        }
                        burn_slope_row(r);
                    }
                }
                fclose(sf);
            }
        }
        if (png_ok) { pw.close(); cout << "[preprocess] png → ./output/preprocess/slope.png\n"; }
        if (any) cout << "[preprocess] written slope bit\n";
    }
    // FIX: this closing brace used to be at the very end of the function
    // (right before Global::preprocess's own closing brace) instead of
    // here. That put the constraint-layer loop (`for (auto& job : jobs)`),
    // the map.bin/map_tiles write, and the GDALClose(mat)/remove(mat_tmp)
    // cleanup ALL inside `if (individual == 6 || individual == -1)` — the
    // gate that's only supposed to cover elevation+slope. Every individual
    // layer (--road/--rail/--water-areas/--water-lines/-i, i.e.
    // individual == 1..5) fell outside that condition, so none of them
    // ever reached the loop that actually rasterizes their shapefile —
    // and `mat`'s GDAL handle plus matrix_tmp.tif were never closed or
    // removed on those runs either, since that cleanup was unreachable
    // too. Closing the gate here (right after elevation+slope) lets
    // everything below run for every value of `individual`; the pieces
    // that genuinely are -1-only or 6-only already have their own
    // `if (individual == -1)` / `if (individual == 6)` checks further
    // down (e.g. line ~784's map.bin write, line ~868's elevTime print).
    }   // closes: if (individual == 6 || individual == -1)

    if(timer){
       end = chrono::steady_clock::now();
       elevTime = chrono::duration_cast<chrono::seconds>(end - start).count();
       start = chrono::steady_clock::now();
    }

    // ── constraint layers ────────────────────────────────────────────────
    for (auto& job : jobs) {
        cout << "\n=== " << job.name << " (buffer=" << job.distance << "m) ===\n";
        if (!lib::is_7755(job.shp)) {
            cout << "[preprocess] reprojecting " << job.shp << "...\n";
            if (lib::reproject(job.shp) != 0) {
                cerr << "[preprocess] ERROR: reprojection failed\n";
                GDALClose(mat); remove(mat_tmp.c_str());
                return;
            }
        }

        const string job_tmp = path + job.name + "_tmp.tif";
        size_t burned = 0;
        const int job_radius = job.distance / RES;
        const GDALDataType job_dt = raster_type_for_radius(job_radius);
        if (job_dt == GDT_UInt16)
            cout << "[preprocess] radius " << job_radius << " cells > 254 — using a 16-bit "
                    "job raster for this layer (~" << (2.0 * N / 1e9) << " GB temp)\n";
        if (!rasterize_layer_to_raster(job.shp, ROWS, COLS, ox, oy, RES,
                                        job.boundary_only, job_tmp, burned, job_dt)) {
            cerr << "[preprocess] SKIP: " << job.name << "\n";
            remove(job_tmp.c_str());
            continue;
        }

        cout << "[preprocess] raw hits   : " << burned
             << " (" << (burned * 100.0 / N) << "%)\n";
        if (burned == 0)
            cerr << "[preprocess] WARN: no cells burned — geotransform or SHP issue\n";

        size_t expanded = 0;
        if (!expand_raster(job_tmp, ROWS, COLS, job_radius, expanded)) {
            remove(job_tmp.c_str());
            continue;
        }

        GDALDataset* jds = (GDALDataset*)GDALOpenEx(
            job_tmp.c_str(), GDAL_OF_RASTER, nullptr, nullptr, nullptr);
        if (!jds) { remove(job_tmp.c_str()); continue; }

        long long final_hits = 0;
        if (!or_mask_into_matrix(mat, jds, ROWS, COLS, job.index, total, final_hits))
            cerr << "[preprocess] WARN: matrix accumulation failed for " << job.name << "\n";

        feasible += final_hits;
        cout << "[preprocess] final hits : " << final_hits
             << " (" << (final_hits * 100.0 / N) << "%)\n";

        if (seperate || individual)
            mask_raster_to_bin(jds, path + job.name + ".bin", ROWS, COLS);
        if (images)
            mask_raster_to_png(jds, "./output/preprocess/" + job.name + ".png", ROWS, COLS);

        GDALClose(jds);
        remove(job_tmp.c_str());   // free the disk before the next job
    }

    // ── write the matrix out ─────────────────────────────────────────────
    // block_size > 0  → map_tiles/ (same layout as elev_tiles/)
    // block_size == 0 → ONE flat map.bin, the format that predates tiling
    // Either way the result is re-checked with the SAME code the loader uses,
    // so "preprocess said OK" and "GA can load it" can't disagree.
    bool map_ok = true;
    if(individual == -1) {
        GDALFlushCache(mat);
        namespace fs = std::filesystem;
        error_code fec;
        const string tdir = tilegrid::tile_dir_for(opath);

        if (tiled_out) {
            make_dir(tdir);
            cout << "\n[preprocess] writing tiles → " << tdir << "\n";

            tilegrid::TileRowWriter<uint8_t> tw;
            vector<uint8_t> row((size_t)COLS);
            int cur_tile_row = -1;
            bool bad = false;

            for (int r = 0; r < ROWS && !bad; r++) {
                const int tr = r / block_size;
                if (tr != cur_tile_row) {
                    const int row0 = tr * block_size;
                    const int h    = min(block_size, ROWS - row0);
                    if (!tw.open(tdir, "map", tr, row0, h, COLS, block_size)) { bad = true; break; }
                    cur_tile_row = tr;
                }
                if (!row_read(mat, r, row.data(), COLS)) { bad = true; break; }
                if (!tw.write_row(row.data()))           { bad = true; break; }
            }
            tw.close();

            if (!bad) {
                bad = !tilegrid::write_grid_meta(tdir, ROWS, COLS, block_size,
                                                 (int)sizeof(uint8_t), RES,
                                                 min_x, max_x, min_y, max_y);
            }
            if (!bad) {
                tilegrid::Grid vg;   // what the loader will see
                if (!vg.discover(tdir) || vg.rows != ROWS || vg.cols != COLS ||
                    vg.item_size != (int)sizeof(uint8_t)) {
                    cerr << "[preprocess] ERROR: map tiles written but not loadable as a "
                         << COLS << " x " << ROWS << " uint8 grid\n";
                    bad = true;
                } else {
                    cout << "[preprocess] map tiles → " << tdir << " ("
                         << vg.tiles.size() << " tiles, " << N / 1e6 << " MB total) — verified\n";
                }
            }
            map_ok = !bad;

            if (fs::exists(opath, fec))
                cout << "[preprocess] NOTE: a flat " << opath << " from an earlier block_size=0 run "
                        "is still on disk; block_size>0 makes the loader use the tiles, "
                        "so it is only a fallback.\n";
        } else {
            map_ok = mask_raster_to_bin(mat, opath, ROWS, COLS);
            if (map_ok) {
                const uintmax_t got = fs::file_size(opath, fec);
                if (fec || got != (uintmax_t)N) {
                    cerr << "[preprocess] ERROR: " << opath << " is " << got
                         << " bytes, expected " << (uintmax_t)N << "\n";
                    map_ok = false;
                } else {
                    cout << "[preprocess] map.bin verified — " << COLS << " x " << ROWS
                         << " single flat bin\n";
                }
            }
            if (fs::exists(tdir, fec))
                cout << "[preprocess] NOTE: " << tdir << " from an earlier tiled run is still on "
                        "disk; block_size=0 makes the loader use the flat bin, so it is only a "
                        "fallback.\n";
        }

        if (images)
            mask_raster_to_png(mat, "./output/preprocess/combined.png", ROWS, COLS);
    }

    GDALClose(mat);
    remove(mat_tmp.c_str());

    int constraintTime = 0;
    if (timer){
        end = chrono::steady_clock::now();
        constraintTime = chrono::duration_cast<chrono::seconds>(end - start).count();
    }

    if (timer) {
        if (individual == -1) cout << "\n>>> Elevation Build took: " << elevTime << " seconds"<< (images ? " (with png writing)" : "") << endl;
        cout << "\n>>> Combined bin took: " << constraintTime << " seconds" << (images ? " (with png writing)" : "") << endl;
    }
    cout<< "total cells: " << total <<endl;
    cout<< "feasible cells: " << feasible << "i.e. " << ((double)feasible / total) * 100 << "%" <<endl;
    if (!map_ok)
        cerr << "\n>>> FAILED: the constraints map (map.bin / map_tiles) was NOT written "
                "correctly — GA cannot load it. See the errors above.\n";
    else
        cout << "\n>>> Completed All Tasks Succesfully..\n";
}

void init_prequisites (){
    filesystem::create_directories("./output/ga");
    filesystem::create_directories("./output/heatmap");
    filesystem::create_directories("./output/radar");
    filesystem::create_directories("./output/viewshed");
    filesystem::create_directories("./output/path");
    filesystem::create_directories("./output/pso");
    filesystem::create_directories("./output/radarfdr");
    filesystem::create_directories("./output/simdata");
    filesystem::create_directories("./output/gbadws");
    filesystem::create_directories("./output/preprocess");
}
