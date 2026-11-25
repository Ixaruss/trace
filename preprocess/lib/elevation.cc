#include "../bin.h"
#include <cstdint>
#include <gdal_utils.h>
#include <gdalwarper.h>
#include <gdal_alg.h>
#include <vector>
#include <algorithm>
#include <string>

using namespace std;

static bool make_dir(const string& path) {
#ifdef _WIN32
    int ret = system(("mkdir \"" + path + "\" 2>nul").c_str());
#else
    int ret = system(("mkdir -p \"" + path + "\"").c_str());
#endif
    (void)ret;
    return true;
}

// ── tiled output layout ────────────────────────────────────────────────────
//
// When a block size is given, elevation/slope are written as a FOLDER of
// fixed-size tiles instead of one monolithic .bin, so a consumer can read
// only the tiles overlapping its area of interest.
//
// All metadata lives in the filenames — no manifest to parse, no sidecar to
// keep in sync. A loader builds its whole index from a directory listing.
//
//   <prefix>_tr<TILEROW>_tc<TILECOL>_r<ROW0>_c<COL0>_h<ROWS>_w<COLS>.bin
//
//     tr / tc   tile index in the tile grid (AOI → tile is integer division)
//     r  / c    global row/col of this tile's top-left cell
//     h  / w    ACTUAL cell dims — edge tiles are smaller than the block
//
// Each tile is plain row-major with stride = w (NOT the global width), so
// reading one tile is a single contiguous fread of h*w*sizeof(T) bytes.
//
// Alongside the tiles, a marker file carries grid-level facts in its name:
//
//   __grid_R<ROWS>_C<COLS>_B<BLOCK>_TR<NTILEROWS>_TC<NTILECOLS>_I<ITEMSIZE>_res<M>.meta
//
// so a loader can validate a folder (and learn the element size) without
// opening a single tile. Its contents repeat the same facts as text plus
// the projected extent.

static string tile_dir_for(const string& bin_path, const char* suffix) {
    // "./bin/elev.bin" → "./bin/elev_tiles/"
    string base = bin_path;
    auto dot = base.find_last_of('.');
    auto sep = base.find_last_of("\\/");
    if (dot != string::npos && (sep == string::npos || dot > sep))
        base = base.substr(0, dot);
    return base + suffix + "/";
}

// Streaming writer for one tile-row.
//
// Rows within a tile are contiguous in the file, and the producer emits
// rows strictly top-to-bottom, so we can hold this tile-row's files open
// and simply append each row's slice to the matching file. That keeps peak
// memory at ONE row regardless of block size — an earlier version buffered
// a full `block`-row band, which at block=10240 meant ~2GB (and ~8.6GB on a
// 30m grid), defeating the purpose of tiling in the first place.
//
// Open handles = number of tile columns (e.g. 4 at 100m, ~14 at 30m).
template <typename T>
struct TileRowWriter {
    vector<FILE*> files;
    vector<int>   col0;
    vector<int>   width;
    bool          ok = false;

    bool open(const string& dir, const char* prefix,
              int tile_row, int row0, int rows_in_band,
              int OUT_COLS, int block) {
        close();
        const int n_tile_cols = (OUT_COLS + block - 1) / block;
        files.assign(n_tile_cols, nullptr);
        col0.assign(n_tile_cols, 0);
        width.assign(n_tile_cols, 0);

        for (int tc = 0; tc < n_tile_cols; tc++) {
            const int c0 = tc * block;
            const int w  = min(block, OUT_COLS - c0);
            col0[tc]  = c0;
            width[tc] = w;

            char name[320];
            snprintf(name, sizeof(name),
                     "%s_tr%05d_tc%05d_r%09d_c%09d_h%05d_w%05d.bin",
                     prefix, tile_row, tc, row0, c0, rows_in_band, w);

            const string p = dir + name;
            files[tc] = fopen(p.c_str(), "wb");
            if (!files[tc]) {
                cerr << "[elevation] ERROR: could not open tile " << p << " for writing\n";
                close();
                return false;
            }
        }
        ok = true;
        return true;
    }

    bool write_row(const T* full_row) {
        if (!ok) return false;
        for (size_t tc = 0; tc < files.size(); tc++) {
            const size_t w   = (size_t)width[tc];
            const size_t got = fwrite(full_row + col0[tc], sizeof(T), w, files[tc]);
            if (got != w) {
                cerr << "[elevation] ERROR: short write on tile column " << tc << "\n";
                ok = false;
                return false;
            }
        }
        return true;
    }

    void close() {
        for (auto* f : files) if (f) fclose(f);
        files.clear(); col0.clear(); width.clear();
        ok = false;
    }

    ~TileRowWriter() { close(); }
};

static bool write_grid_meta(const string& dir,
                             int OUT_ROWS, int OUT_COLS,
                             int block, int item_size, int resol,
                             double min_x, double max_x,
                             double min_y, double max_y) {
    const int n_tile_rows = (OUT_ROWS + block - 1) / block;
    const int n_tile_cols = (OUT_COLS + block - 1) / block;

    char name[320];
    snprintf(name, sizeof(name),
             "__grid_R%d_C%d_B%d_TR%d_TC%d_I%d_res%d.meta",
             OUT_ROWS, OUT_COLS, block, n_tile_rows, n_tile_cols, item_size, resol);

    const string p = dir + name;
    FILE* f = fopen(p.c_str(), "wb");
    if (!f) {
        cerr << "[elevation] ERROR: could not write grid meta " << p << "\n";
        return false;
    }
    fprintf(f, "rows=%d\ncols=%d\nblock=%d\ntile_rows=%d\ntile_cols=%d\n"
               "item_size=%d\nresolution=%d\n"
               "min_x=%.6f\nmax_x=%.6f\nmin_y=%.6f\nmax_y=%.6f\nepsg=7755\n",
            OUT_ROWS, OUT_COLS, block, n_tile_rows, n_tile_cols,
            item_size, resol, min_x, max_x, min_y, max_y);
    fclose(f);
    cout << "[elevation] meta → " << p << "\n";
    return true;
}

// Horn's formula, nodata-aware, operating on a 3-row sliding window
// (north/center/south row pointers) rather than the whole grid, so the
// caller never needs to hold more than a few rows in memory at once.
// Nodata neighbors are skipped — falls back to one-sided diff.
// If both sides of an axis are nodata → slope = -9999.0f
// N_row/S_row may be nullptr (off the top/bottom edge of the grid).
static float compute_slope_row(const int16_t* N_row,
                                 const int16_t* C_row,
                                 const int16_t* S_row,
                                 int c, int COLS, int resolution) {
    int16_t center = C_row[c];
    if (center == -9999) return -9999.0f;

    auto val = [&](const int16_t* row, int cc) -> int16_t {
        if (!row || cc < 0 || cc >= COLS) return -9999;
        return row[cc];
    };

    int16_t Wv = val(C_row, c - 1);
    int16_t Ev = val(C_row, c + 1);
    int16_t Nv = val(N_row, c);
    int16_t Sv = val(S_row, c);

    // ── dz/dx (west→east) ─────────────────────────────────────────────────
    float dzdx;
    bool e_ok = (Ev != -9999);
    bool w_ok = (Wv != -9999);

    if      (e_ok && w_ok)  dzdx = (Ev - Wv) / (2.0f * resolution);
    else if (e_ok)          dzdx = (Ev - center) / (float)resolution;
    else if (w_ok)          dzdx = (center - Wv) / (float)resolution;
    else                    return -9999.0f;   // both sides nodata

    // ── dz/dy (north→south) ───────────────────────────────────────────────
    float dzdy;
    bool n_ok = (Nv != -9999);
    bool s_ok = (Sv != -9999);

    if      (n_ok && s_ok)  dzdy = (Nv - Sv) / (2.0f * resolution);
    else if (n_ok)          dzdy = (center - Nv) / (float)resolution;
    else if (s_ok)          dzdy = (Sv - center) / (float)resolution;
    else                    return -9999.0f;   // both sides nodata

    // ── Horn's slope in degrees ────────────────────────────────────────────
    return atan(sqrt(dzdx * dzdx + dzdy * dzdy)) * (180.0f / (float)M_PI);
}

// block_size: edge length in cells of each output tile. Pass <= 0 to keep
// the legacy single-file output (one flat elev.bin / slope.bin).
//
// Sizing guide (tile bytes = block^2 * itemsize; slope's 4-byte float is
// the binding constraint since both folders share one block size):
//   block  8192 → elev 134MB, slope 268MB
//   block 10240 → elev 210MB, slope 419MB   <- good default for 200-500MB
//   block 11264 → elev 254MB, slope 507MB
// Peak RAM does NOT depend on this value — see TileRowWriter above.
int Bin::build_elevation_matrix(bool split, int block_size) {
    if (!Config::is_available()) {
        cerr << "[elevation] ERROR: no config available — output paths, input "
                "elevation path and grid resolution all come from it. Refusing "
                "to fall back to compiled-in DEFAULT_*/RESOLUTION placeholders.\n";
        return -1;
    }
    Config cfg0 = Config::load();
    string opath = cfg0.output.elevation_file;
    string ipath = cfg0.input.elevation;
    string spath = cfg0.output.slope_file;
    int    resol = cfg0.grid.resolution;
    if (resol <= 0) {
        cerr << "[elevation] ERROR: config grid.resolution is not set (or <= 0).\n";
        return -1;
    }

    // block_size < 0 means "caller didn't say" -> follow cfg.grid.block_size so
    // elevation/slope use the same layout preprocess will use for map.bin.
    if (block_size < 0) block_size = cfg0.grid.block_size;
    if (block_size < 0) block_size = 0;
    auto idx = opath.find_last_of("\\/");
    string path;

    if(std::string::npos != idx){
        path = opath.substr(0,idx+1);
    }
    make_dir(path);
    GDALAllRegister();

    const bool chunked = (block_size > 0);
    if (chunked)
        cout << "[elevation] Tiled output enabled, block = " << block_size << " cells\n";
    else
        cout << "[elevation] Tiled output disabled — writing single flat bin(s)\n";

    // ── 1. open source TIF ─────────────────────────────────────────────────
    cout << "[elevation] Opening: " << ipath << "\n";
    GDALDataset* src_ds = (GDALDataset*)GDALOpenEx(
        ipath.c_str(), GDAL_OF_RASTER, nullptr, nullptr, nullptr);
    if (!src_ds) {
        cerr << "[elevation] ERROR: could not open " << ipath << "\n";
        return -1;
    }
    cout << "[elevation] Source size: "
         << src_ds->GetRasterXSize() << " x " << src_ds->GetRasterYSize()
         << "  bands: " << src_ds->GetRasterCount() << "\n";

    // ── 1b. read source nodata value from the band itself ──────────────────
    int    has_src_nd = 0;
    double src_nd      = src_ds->GetRasterBand(1)->GetNoDataValue(&has_src_nd);
    if (has_src_nd)
        cout << "[elevation] Source NoData: " << src_nd << "\n";
    else
        cout << "[elevation] Source has no NoData value set — warp will not mask any source pixels\n";

    // ── 2. build source/dest SRS ────────────────────────────────────────────
    const char* src_wkt = src_ds->GetProjectionRef();
    OGRSpatialReference src_srs, dst_srs;
    if (src_srs.importFromWkt(src_wkt) != OGRERR_NONE) {
        cerr << "[elevation] ERROR: could not parse source WKT\n";
        GDALClose(src_ds);
        return -1;
    }
    dst_srs.importFromEPSG(7755);
    src_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    dst_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

    char* dst_wkt_str = nullptr;
    dst_srs.exportToWkt(&dst_wkt_str);

    // ── 3. derive EPSG:7755 extent + output size via GDALSuggestedWarpOutput
    char** transform_opts = nullptr;
    transform_opts = CSLSetNameValue(transform_opts, "DST_SRS", dst_wkt_str);
    void* sugg_transformer_arg = GDALCreateGenImgProjTransformer2(
        (GDALDatasetH)src_ds, nullptr, transform_opts);
    CSLDestroy(transform_opts);
    if (!sugg_transformer_arg) {
        cerr << "[elevation] ERROR: could not create transformer for extent estimation\n";
        CPLFree(dst_wkt_str);
        GDALClose(src_ds);
        return -1;
    }

    double sugg_gt[6];
    int    sugg_cols = 0, sugg_rows = 0;
    if (GDALSuggestedWarpOutput(src_ds, GDALGenImgProjTransform, sugg_transformer_arg,
                                 sugg_gt, &sugg_cols, &sugg_rows) != CE_None) {
        cerr << "[elevation] ERROR: GDALSuggestedWarpOutput failed\n";
        GDALDestroyGenImgProjTransformer(sugg_transformer_arg);
        CPLFree(dst_wkt_str);
        GDALClose(src_ds);
        return -1;
    }
    GDALDestroyGenImgProjTransformer(sugg_transformer_arg);

    double proj_min_x = sugg_gt[0];
    double proj_max_y = sugg_gt[3];
    double proj_max_x = sugg_gt[0] + sugg_cols * sugg_gt[1];
    double proj_min_y = sugg_gt[3] + sugg_rows * sugg_gt[5];

    cout << "[elevation] Projected extent (EPSG:7755):\n"
         << "  X: " << proj_min_x << " → " << proj_max_x << "\n"
         << "  Y: " << proj_min_y << " → " << proj_max_y << "\n";

    // ── 4. compute output grid at same resol as vector matrix ─────────
    int    OUT_COLS = (int)ceil((proj_max_x - proj_min_x) / resol);
    int    OUT_ROWS = (int)ceil((proj_max_y - proj_min_y) / resol);
    double ORIGIN_X = proj_min_x;
    double ORIGIN_Y = proj_max_y;

    cout << "[elevation] Output grid: " << OUT_COLS << " x " << OUT_ROWS
         << " @ " << resol << "m/cell\n";

    if (chunked) {
        const int ntr = (OUT_ROWS + block_size - 1) / block_size;
        const int ntc = (OUT_COLS + block_size - 1) / block_size;
        cout << "[elevation] Tile grid  : " << ntc << " x " << ntr
             << " = " << (ntr * ntc) << " tiles of "
             << block_size << "x" << block_size << "\n";
        cout << "[elevation] Open handles per tile row: " << ntc << "\n";
    }

    // ── 5. create warp destination as a tiled GTiff on disk, not MEM ───────
    string warped_tmp_path = path + "elev_warp_tmp.tif";

    GDALDriver* dst_driver = GetGDALDriverManager()->GetDriverByName("GTiff");
    char** create_opts = nullptr;
    create_opts = CSLSetNameValue(create_opts, "BIGTIFF", "YES");
    create_opts = CSLSetNameValue(create_opts, "TILED", "YES");
    create_opts = CSLSetNameValue(create_opts, "BLOCKXSIZE", "512");
    create_opts = CSLSetNameValue(create_opts, "BLOCKYSIZE", "512");
    create_opts = CSLSetNameValue(create_opts, "COMPRESS", "DEFLATE");
    create_opts = CSLSetNameValue(create_opts, "NUM_THREADS", "ALL_CPUS");

    GDALDataset* dst_ds = dst_driver->Create(
        warped_tmp_path.c_str(), OUT_COLS, OUT_ROWS, 1, GDT_Int16, create_opts);
    CSLDestroy(create_opts);
    if (!dst_ds) {
        cerr << "[elevation] ERROR: failed to create warp destination " << warped_tmp_path << "\n";
        CPLFree(dst_wkt_str);
        GDALClose(src_ds);
        return -1;
    }

    double gt[6] = { ORIGIN_X, (double)resol, 0, ORIGIN_Y, 0, -(double)resol };
    dst_ds->SetGeoTransform(gt);
    dst_ds->SetProjection(dst_wkt_str);
    dst_ds->GetRasterBand(1)->SetNoDataValue(-9999);
    CPLFree(dst_wkt_str);

    // ── 6. warp with multithreading ────────────────────────────────────────
    cout << "[elevation] Warping to EPSG:7755 @ " << resol << "m...\n";

    GDALWarpOptions* warp_opts      = GDALCreateWarpOptions();
    warp_opts->hSrcDS               = src_ds;
    warp_opts->hDstDS               = dst_ds;
    warp_opts->nBandCount           = 1;
    warp_opts->panSrcBands          = (int*)CPLMalloc(sizeof(int));
    warp_opts->panDstBands          = (int*)CPLMalloc(sizeof(int));
    warp_opts->panSrcBands[0]       = 1;
    warp_opts->panDstBands[0]       = 1;
    warp_opts->eResampleAlg         = GRA_Bilinear;
    warp_opts->padfDstNoDataReal    = (double*)CPLMalloc(sizeof(double));
    warp_opts->padfDstNoDataReal[0] = -9999.0;

    if (has_src_nd) {
        warp_opts->padfSrcNoDataReal    = (double*)CPLMalloc(sizeof(double));
        warp_opts->padfSrcNoDataReal[0] = src_nd;
    } else {
        warp_opts->padfSrcNoDataReal = nullptr;
    }

    warp_opts->papszWarpOptions = CSLSetNameValue(warp_opts->papszWarpOptions, "INIT_DEST", "NO_DATA");
    warp_opts->papszWarpOptions = CSLSetNameValue(warp_opts->papszWarpOptions, "NUM_THREADS", "ALL_CPUS");
    warp_opts->pTransformerArg      = GDALCreateGenImgProjTransformer(
        src_ds, nullptr, dst_ds, nullptr, FALSE, 0, 1);
    warp_opts->pfnTransformer       = GDALGenImgProjTransform;

    GDALWarpOperation warp_op;
    if (warp_op.Initialize(warp_opts) != CE_None) {
        cerr << "[elevation] ERROR: warp init failed\n";
        GDALDestroyGenImgProjTransformer(warp_opts->pTransformerArg);
        GDALDestroyWarpOptions(warp_opts);
        GDALClose(dst_ds); GDALClose(src_ds);
        return -1;
    }
    if ((CPLErr)warp_op.ChunkAndWarpImage(0, 0, OUT_COLS, OUT_ROWS) != CE_None) {
        cerr << "[elevation] ERROR: warp failed\n";
        GDALDestroyGenImgProjTransformer(warp_opts->pTransformerArg);
        GDALDestroyWarpOptions(warp_opts);
        GDALClose(dst_ds); GDALClose(src_ds);
        return -1;
    }
    GDALDestroyGenImgProjTransformer(warp_opts->pTransformerArg);
    GDALDestroyWarpOptions(warp_opts);
    GDALClose(src_ds);

    // dst_ds MUST be closed before reopening for read — GDAL's GTiff writer
    // buffers tiles and only guarantees a flush to disk on GDALClose().
    GDALClose(dst_ds);
    cout << "[elevation] Warp done\n";

    dst_ds = (GDALDataset*)GDALOpenEx(
        warped_tmp_path.c_str(), GDAL_OF_RASTER, nullptr, nullptr, nullptr);
    if (!dst_ds) {
        cerr << "[elevation] ERROR: could not reopen warped raster " << warped_tmp_path << "\n";
        return -1;
    }
    GDALRasterBand* band = dst_ds->GetRasterBand(1);

    // ── 7. set up output sinks (tiled folders, or legacy flat files) ───────
    string elev_dir, slope_dir;
    FILE*  f_elev  = nullptr;
    FILE*  f_slope = nullptr;

    if (chunked) {
        elev_dir = tile_dir_for(opath, "_tiles");
        make_dir(elev_dir);
        cout << "[elevation] Tile dir   : " << elev_dir << "\n";
        if (split) {
            slope_dir = tile_dir_for(spath, "_tiles");
            make_dir(slope_dir);
            cout << "[elevation] Tile dir   : " << slope_dir << "\n";
        }
    } else {
        f_elev = fopen(opath.c_str(), "wb");
        if (!f_elev) {
            cerr << "[elevation] ERROR: could not open " << opath << " for writing\n";
            GDALClose(dst_ds);
            return -1;
        }
        if (split) {
            f_slope = fopen(spath.c_str(), "wb");
            if (!f_slope) {
                cerr << "[elevation] ERROR: could not open " << spath << " for writing\n";
                fclose(f_elev);
                GDALClose(dst_ds);
                return -1;
            }
        }
    }

    // ── 8. stream elevation + slope with a 3-row sliding window ────────────
    cout << "[elevation] Computing slope (Horn's formula, nodata-aware, streamed)...\n";

    size_t total_cells  = (size_t)OUT_ROWS * OUT_COLS;
    size_t raw_nodata    = 0, raw_zero = 0, nodata_slope = 0;
    int16_t raw_min = 32767, raw_max = -32768;

    vector<int16_t> row_buf[3];
    for (int i = 0; i < 3; i++) row_buf[i].resize(OUT_COLS);
    vector<float>    slope_row(OUT_COLS);
    vector<ElevCell> combined_row;
    if (!split) combined_row.resize(OUT_COLS);

    // tile-row writers (chunked mode only) — hold open the files of the
    // tile row currently being produced, and roll over every block_size rows
    TileRowWriter<ElevCell> tw_combined;
    TileRowWriter<int16_t>  tw_elev;
    TileRowWriter<float>    tw_slope;
    int  cur_tile_row  = -1;
    bool write_failed  = false;

    auto read_row = [&](int r, int16_t* buf) {
        if (r < 0 || r >= OUT_ROWS) {
            std::fill(buf, buf + OUT_COLS, (int16_t)-9999);
            return;
        }
        if (band->RasterIO(GF_Read, 0, r, OUT_COLS, 1,
                            buf, OUT_COLS, 1, GDT_Int16, 0, 0) != CE_None) {
            cerr << "[elevation] WARN: row read failed at r=" << r << "\n";
            std::fill(buf, buf + OUT_COLS, (int16_t)-9999);
        }
    };

    read_row(-1, row_buf[0].data());
    read_row(0,  row_buf[1].data());
    read_row(1,  row_buf[2].data());

    for (int r = 0; r < OUT_ROWS && !write_failed; r++) {
        if (r % 5000 == 0 && r > 0)
            cout << "[elevation] slope row " << r << "/" << OUT_ROWS << "\n";

        // roll over to the next tile row when we cross a block boundary
        if (chunked) {
            const int tr = r / block_size;
            if (tr != cur_tile_row) {
                const int row0         = tr * block_size;
                const int rows_in_band = min(block_size, OUT_ROWS - row0);
                bool okopen = true;
                if (!split) {
                    okopen = tw_combined.open(elev_dir, "elev", tr, row0,
                                               rows_in_band, OUT_COLS, block_size);
                } else {
                    okopen = tw_elev.open(elev_dir, "elev", tr, row0,
                                           rows_in_band, OUT_COLS, block_size);
                    if (okopen)
                        okopen = tw_slope.open(slope_dir, "slope", tr, row0,
                                                rows_in_band, OUT_COLS, block_size);
                }
                if (!okopen) { write_failed = true; break; }
                cur_tile_row = tr;
            }
        }

        const int16_t* Nrow = row_buf[0].data();
        const int16_t* Crow = row_buf[1].data();
        const int16_t* Srow = row_buf[2].data();

        for (int c = 0; c < OUT_COLS; c++) {
            int16_t center = Crow[c];
            if (center == -9999) raw_nodata++;
            if (center == 0)     raw_zero++;
            raw_min = min(raw_min, center);
            raw_max = max(raw_max, center);

            float s = compute_slope_row(Nrow, Crow, Srow, c, OUT_COLS, resol);
            if (s == -9999.0f) nodata_slope++;

            if (!split) {
                combined_row[c].elevation = center;
                combined_row[c].slope     = s;
            } else {
                slope_row[c] = s;
            }
        }

        if (!chunked) {
            if (!split) {
                fwrite(combined_row.data(), sizeof(ElevCell), OUT_COLS, f_elev);
            } else {
                fwrite(Crow, sizeof(int16_t), OUT_COLS, f_elev);
                fwrite(slope_row.data(), sizeof(float), OUT_COLS, f_slope);
            }
        } else {
            if (!split) {
                if (!tw_combined.write_row(combined_row.data())) write_failed = true;
            } else {
                if (!tw_elev.write_row(Crow))                     write_failed = true;
                if (!write_failed && !tw_slope.write_row(slope_row.data()))
                    write_failed = true;
            }
        }

        // slide the window: [N,C,S] → [C,S,<next>]
        std::swap(row_buf[0], row_buf[1]);
        std::swap(row_buf[1], row_buf[2]);
        read_row(r + 2, row_buf[2].data());
    }

    tw_combined.close();
    tw_elev.close();
    tw_slope.close();
    if (f_elev)  fclose(f_elev);
    if (f_slope) fclose(f_slope);
    GDALClose(dst_ds);
    remove(warped_tmp_path.c_str());   // done with the intermediate warp file

    if (write_failed) {
        cerr << "[elevation] ERROR: aborted — tile write failed\n";
        return -1;
    }

    cout << "[elevation] RAW elevation: nodata=" << raw_nodata
         << " zero=" << raw_zero
         << " min=" << raw_min << " max=" << raw_max << "\n";
    cout << "[elevation] Slope done. Nodata slope cells: " << nodata_slope
         << " (" << (nodata_slope * 100.0 / total_cells) << "%)\n";

    // ── 9. grid-level metadata marker(s) + summary ─────────────────────────
    if (chunked) {
        const int n_tiles = ((OUT_ROWS + block_size - 1) / block_size) *
                            ((OUT_COLS + block_size - 1) / block_size);
        if (!split) {
            write_grid_meta(elev_dir, OUT_ROWS, OUT_COLS, block_size,
                            (int)sizeof(ElevCell), resol,
                            proj_min_x, proj_max_x, proj_min_y, proj_max_y);
            cout << "[elevation] Written → " << elev_dir << " ("
                 << n_tiles << " tiles, "
                 << (total_cells * sizeof(ElevCell)) / 1e6 << " MB total, "
                 << sizeof(ElevCell) << " bytes/cell)\n";
        } else {
            write_grid_meta(elev_dir, OUT_ROWS, OUT_COLS, block_size,
                            (int)sizeof(int16_t), resol,
                            proj_min_x, proj_max_x, proj_min_y, proj_max_y);
            write_grid_meta(slope_dir, OUT_ROWS, OUT_COLS, block_size,
                            (int)sizeof(float), resol,
                            proj_min_x, proj_max_x, proj_min_y, proj_max_y);
            cout << "[elevation] Written → " << elev_dir << " ("
                 << n_tiles << " tiles, "
                 << (total_cells * sizeof(int16_t)) / 1e6 << " MB total)\n"
                 << "[elevation] Written → " << slope_dir << " ("
                 << n_tiles << " tiles, "
                 << (total_cells * sizeof(float)) / 1e6 << " MB total)\n";
        }
    } else {
        if (!split) {
            cout << "[elevation] Written → " << opath << " ("
                 << (total_cells * sizeof(ElevCell)) / 1e6 << " MB, "
                 << sizeof(ElevCell) << " bytes/cell)\n";
        } else {
            cout << "[elevation] Written → " << opath << " ("
                 << (total_cells * sizeof(int16_t)) / 1e6 << " MB)\n"
                 << "[elevation] Written → " << spath << " ("
                 << (total_cells * sizeof(float)) / 1e6 << " MB)\n";
        }
    }

    // ── 10. save derived grid params to config ────────────────────────────
    // Config was required at entry, so this always runs — not an
    // if(available)/else(warn) fallback branch.
    {
        Config cfg                = Config::load();
        cfg.grid.rows             = OUT_ROWS;
        cfg.grid.cols             = OUT_COLS;
        cfg.grid.origin_x         = (float)ORIGIN_X;
        cfg.grid.origin_y         = (float)ORIGIN_Y;
        cfg.grid.min_x            = proj_min_x;
        cfg.grid.max_x            = proj_max_x;
        cfg.grid.min_y            = proj_min_y;
        cfg.grid.max_y            = proj_max_y;
        cfg.save();
        cout << "[elevation] Config updated:\n"
             << "  elev grid : " << OUT_COLS << " x " << OUT_ROWS << "\n"
             << "  min_x     : " << proj_min_x << "\n"
             << "  max_x     : " << proj_max_x << "\n"
             << "  min_y     : " << proj_min_y << "\n"
             << "  max_y     : " << proj_max_y << "\n";
    }

    return 0;
}
