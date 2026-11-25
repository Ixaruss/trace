#ifndef TILEGRID_H_
#define TILEGRID_H_

// ══════════════════════════════════════════════════════════════════════════
// Tiled grid reader — shared by preprocess (sequential full-grid streaming)
// and the Bin loader (windowed AOI reads).
//
// Layout written by Bin::build_elevation_matrix(split, block_size):
//
//   __grid_R<ROWS>_C<COLS>_B<BLOCK>_TR<NTR>_TC<NTC>_I<ITEMSIZE>_res<M>.meta
//   <prefix>_tr<TR>_tc<TC>_r<ROW0>_c<COL0>_h<H>_w<W>.bin
//
// All metadata lives in the filenames, so discovery is a directory listing
// and nothing has to be parsed out of file contents. Each tile is row-major
// with stride = w (the TILE's width, not the grid's).
// ══════════════════════════════════════════════════════════════════════════

#include <cstdio>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <filesystem>
#include <iostream>

#ifdef _WIN32
  #define TG_FSEEK _fseeki64
  #define TG_FTELL _ftelli64
#else
  #define TG_FSEEK fseeko
  #define TG_FTELL ftello
#endif

namespace tilegrid {

struct TileInfo {
    int         tr = 0, tc = 0;
    int         r0 = 0, c0 = 0;
    int         h  = 0, w  = 0;
    std::string path;
};

// "./bin/elev.bin" → "./bin/elev_tiles/"
inline std::string tile_dir_for(const std::string& bin_path,
                                 const char* suffix = "_tiles") {
    std::string base = bin_path;
    auto dot = base.find_last_of('.');
    auto sep = base.find_last_of("\\/");
    if (dot != std::string::npos && (sep == std::string::npos || dot > sep))
        base = base.substr(0, dot);
    return base + suffix + "/";
}

struct Grid {
    bool        valid      = false;
    int         rows       = 0;
    int         cols       = 0;
    int         block      = 0;
    int         tile_rows  = 0;
    int         tile_cols  = 0;
    int         item_size  = 0;
    int         resolution = 0;
    std::string dir;
    std::vector<TileInfo> tiles;   // [tr * tile_cols + tc]

    const TileInfo* at(int tr, int tc) const {
        if (tr < 0 || tr >= tile_rows || tc < 0 || tc >= tile_cols) return nullptr;
        const TileInfo& t = tiles[(size_t)tr * tile_cols + tc];
        return t.w > 0 ? &t : nullptr;
    }

    // Quiet false when the folder simply isn't there — caller falls back
    // to the flat .bin.
    bool discover(const std::string& d) {
        namespace fs = std::filesystem;
        dir = d;
        std::error_code ec;
        if (!fs::exists(d, ec) || !fs::is_directory(d, ec)) return false;

        bool got_meta = false;
        std::vector<TileInfo> found;

        for (const auto& ent : fs::directory_iterator(d, ec)) {
            if (ec) break;
            if (!ent.is_regular_file()) continue;
            const std::string fn = ent.path().filename().string();

            if (fn.rfind("__grid_", 0) == 0) {
                int R=0,C=0,B=0,TR=0,TC=0,I=0,RES=0;
                if (sscanf(fn.c_str(),
                           "__grid_R%d_C%d_B%d_TR%d_TC%d_I%d_res%d.meta",
                           &R,&C,&B,&TR,&TC,&I,&RES) == 7) {
                    rows = R; cols = C; block = B;
                    tile_rows = TR; tile_cols = TC;
                    item_size = I; resolution = RES;
                    got_meta = true;
                }
                continue;
            }

            const auto p = fn.find("_tr");
            if (p == std::string::npos) continue;
            TileInfo t;
            if (sscanf(fn.c_str() + p, "_tr%d_tc%d_r%d_c%d_h%d_w%d.bin",
                       &t.tr, &t.tc, &t.r0, &t.c0, &t.h, &t.w) == 6) {
                t.path = d + fn;
                found.push_back(std::move(t));
            }
        }

        if (!got_meta) {
            std::cerr << "[tiles] WARN: no __grid_*.meta in " << d << "\n";
            return false;
        }
        if (found.empty()) {
            std::cerr << "[tiles] WARN: no tiles found in " << d << "\n";
            return false;
        }

        tiles.assign((size_t)tile_rows * tile_cols, TileInfo{});
        size_t placed = 0;
        for (auto& t : found) {
            if (t.tr < 0 || t.tr >= tile_rows || t.tc < 0 || t.tc >= tile_cols) {
                std::cerr << "[tiles] WARN: tile outside declared grid: " << t.path << "\n";
                continue;
            }
            tiles[(size_t)t.tr * tile_cols + t.tc] = t;
            placed++;
        }

        const size_t expect = (size_t)tile_rows * tile_cols;
        if (placed != expect) {
            std::cerr << "[tiles] ERROR: incomplete tile folder — " << placed
                      << "/" << expect << " tiles in " << d << "\n";
            return false;
        }

        valid = true;
        return true;
    }

    void print(const char* tag) const {
        std::cout << "[tiles] " << tag << " : " << dir << "\n"
                  << "[tiles]   grid  : " << cols << " x " << rows
                  << "  block=" << block
                  << "  item=" << item_size << "B"
                  << "  tiles=" << tile_cols << "x" << tile_rows
                  << "  res=" << resolution << "m\n";
    }
};

// ── windowed read ─────────────────────────────────────────────────────────
// Reads the global sub-rectangle [r0, r0+wrows) x [c0, c0+wcols) into `out`
// (row-major, stride = wcols). Only tiles overlapping the window are opened,
// and only the overlapping rows/cols are read out of each. Memory is the
// window, nothing more.
template <typename T>
bool read_window(const Grid& g,
                 int r0, int c0, int wrows, int wcols,
                 std::vector<T>& out) {
    if (!g.valid) return false;
    if ((int)sizeof(T) != g.item_size) {
        std::cerr << "[tiles] ERROR: item size mismatch — folder says "
                  << g.item_size << "B, caller wants " << sizeof(T) << "B\n";
        return false;
    }
    if (r0 < 0 || c0 < 0 || wrows <= 0 || wcols <= 0 ||
        r0 + wrows > g.rows || c0 + wcols > g.cols) {
        std::cerr << "[tiles] ERROR: window (" << r0 << "," << c0 << " "
                  << wrows << "x" << wcols << ") outside grid "
                  << g.rows << "x" << g.cols << "\n";
        return false;
    }

    out.assign((size_t)wrows * wcols, T{});

    const int tr_lo = r0 / g.block;
    const int tr_hi = (r0 + wrows - 1) / g.block;
    const int tc_lo = c0 / g.block;
    const int tc_hi = (c0 + wcols - 1) / g.block;

    size_t tiles_touched = 0;

    for (int tr = tr_lo; tr <= tr_hi; tr++) {
        for (int tc = tc_lo; tc <= tc_hi; tc++) {
            const TileInfo* t = g.at(tr, tc);
            if (!t) {
                std::cerr << "[tiles] ERROR: missing tile tr=" << tr << " tc=" << tc << "\n";
                return false;
            }

            // overlap of window and tile, in global coords
            const int gr_lo = std::max(r0, t->r0);
            const int gr_hi = std::min(r0 + wrows, t->r0 + t->h);
            const int gc_lo = std::max(c0, t->c0);
            const int gc_hi = std::min(c0 + wcols, t->c0 + t->w);
            if (gr_lo >= gr_hi || gc_lo >= gc_hi) continue;

            const int run = gc_hi - gc_lo;          // cells per row to read

            FILE* f = fopen(t->path.c_str(), "rb");
            if (!f) {
                std::cerr << "[tiles] ERROR: cannot open " << t->path << "\n";
                return false;
            }

            for (int gr = gr_lo; gr < gr_hi; gr++) {
                const int64_t local_r = gr - t->r0;
                const int64_t local_c = gc_lo - t->c0;
                const int64_t off = (local_r * (int64_t)t->w + local_c) * (int64_t)sizeof(T);

                if (TG_FSEEK(f, off, SEEK_SET) != 0) {
                    std::cerr << "[tiles] ERROR: seek failed in " << t->path << "\n";
                    fclose(f);
                    return false;
                }
                T* dst = out.data() + (size_t)(gr - r0) * wcols + (gc_lo - c0);
                if (fread(dst, sizeof(T), run, f) != (size_t)run) {
                    std::cerr << "[tiles] ERROR: short read in " << t->path << "\n";
                    fclose(f);
                    return false;
                }
            }
            fclose(f);
            tiles_touched++;
        }
    }

    std::cout << "[tiles] window " << wcols << "x" << wrows
              << " @ (" << c0 << "," << r0 << ") from "
              << tiles_touched << " tile(s), "
              << ((double)wrows * wcols * sizeof(T)) / 1e6 << " MB\n";
    return true;
}

// Same window read, but from a legacy flat row-major .bin of the full grid.
// Lets the windowed code paths work before/without a tiled rebuild.
template <typename T>
bool read_window_flat(const std::string& path,
                      int grid_rows, int grid_cols,
                      int r0, int c0, int wrows, int wcols,
                      std::vector<T>& out) {
    if (r0 < 0 || c0 < 0 || wrows <= 0 || wcols <= 0 ||
        r0 + wrows > grid_rows || c0 + wcols > grid_cols) {
        std::cerr << "[tiles] ERROR: flat window outside grid\n";
        return false;
    }

    FILE* f = fopen(path.c_str(), "rb");
    if (!f) {
        std::cerr << "[tiles] ERROR: cannot open " << path << "\n";
        return false;
    }

    TG_FSEEK(f, 0, SEEK_END);
    const int64_t file_bytes = TG_FTELL(f);
    const int64_t expect     = (int64_t)grid_rows * grid_cols * (int64_t)sizeof(T);
    if (file_bytes != expect) {
        std::cerr << "[tiles] ERROR: " << path << " size mismatch — expected "
                  << expect / 1000000 << " MB, got " << file_bytes / 1000000 << " MB\n";
        fclose(f);
        return false;
    }

    out.assign((size_t)wrows * wcols, T{});
    for (int gr = r0; gr < r0 + wrows; gr++) {
        const int64_t off = ((int64_t)gr * grid_cols + c0) * (int64_t)sizeof(T);
        if (TG_FSEEK(f, off, SEEK_SET) != 0) {
            std::cerr << "[tiles] ERROR: seek failed in " << path << "\n";
            fclose(f); return false;
        }
        T* dst = out.data() + (size_t)(gr - r0) * wcols;
        if (fread(dst, sizeof(T), wcols, f) != (size_t)wcols) {
            std::cerr << "[tiles] ERROR: short read in " << path << "\n";
            fclose(f); return false;
        }
    }
    fclose(f);

    std::cout << "[tiles] flat window " << wcols << "x" << wrows
              << " @ (" << c0 << "," << r0 << "), "
              << ((double)wrows * wcols * sizeof(T)) / 1e6 << " MB\n";
    return true;
}

// ── tile writing ──────────────────────────────────────────────────────────
// Streaming writer for one tile-row. Rows within a tile are contiguous in
// the file and producers emit rows top-to-bottom, so this holds that tile
// row's files open and appends each row's slice to the matching one. Peak
// memory is ONE row regardless of block size — never a full band.
//
// Open handles while active = number of tile columns.
template <typename T>
struct TileRowWriter {
    std::vector<FILE*> files;
    std::vector<int>   col0;
    std::vector<int>   width;
    bool               ok = false;

    bool open(const std::string& dir, const char* prefix,
              int tile_row, int row0, int rows_in_band,
              int OUT_COLS, int block) {
        close();
        const int n_tile_cols = (OUT_COLS + block - 1) / block;
        files.assign(n_tile_cols, nullptr);
        col0.assign(n_tile_cols, 0);
        width.assign(n_tile_cols, 0);

        for (int tc = 0; tc < n_tile_cols; tc++) {
            const int c0 = tc * block;
            const int w  = std::min(block, OUT_COLS - c0);
            col0[tc]  = c0;
            width[tc] = w;

            char name[320];
            snprintf(name, sizeof(name),
                     "%s_tr%05d_tc%05d_r%09d_c%09d_h%05d_w%05d.bin",
                     prefix, tile_row, tc, row0, c0, rows_in_band, w);

            const std::string p = dir + name;
            files[tc] = fopen(p.c_str(), "wb");
            if (!files[tc]) {
                std::cerr << "[tiles] ERROR: cannot open " << p << " for writing\n";
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
            const size_t w = (size_t)width[tc];
            if (fwrite(full_row + col0[tc], sizeof(T), w, files[tc]) != w) {
                std::cerr << "[tiles] ERROR: short write on tile column " << tc << "\n";
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

inline bool write_grid_meta(const std::string& dir,
                             int ROWS, int COLS,
                             int block, int item_size, int resol,
                             double min_x, double max_x,
                             double min_y, double max_y) {
    const int n_tile_rows = (ROWS + block - 1) / block;
    const int n_tile_cols = (COLS + block - 1) / block;

    char name[320];
    snprintf(name, sizeof(name),
             "__grid_R%d_C%d_B%d_TR%d_TC%d_I%d_res%d.meta",
             ROWS, COLS, block, n_tile_rows, n_tile_cols, item_size, resol);

    const std::string p = dir + name;
    FILE* f = fopen(p.c_str(), "wb");
    if (!f) {
        std::cerr << "[tiles] ERROR: could not write grid meta " << p << "\n";
        return false;
    }
    fprintf(f, "rows=%d\ncols=%d\nblock=%d\ntile_rows=%d\ntile_cols=%d\n"
               "item_size=%d\nresolution=%d\n"
               "min_x=%.6f\nmax_x=%.6f\nmin_y=%.6f\nmax_y=%.6f\nepsg=7755\n",
            ROWS, COLS, block, n_tile_rows, n_tile_cols,
            item_size, resol, min_x, max_x, min_y, max_y);
    fclose(f);
    std::cout << "[tiles] meta → " << p << "\n";
    return true;
}

// ── sequential full-grid row streaming ────────────────────────────────────
// Holds one tile row's files open and assembles full-width grid rows. Peak
// memory is one row. Used by preprocess, which genuinely needs every cell.
struct RowReader {
    std::vector<FILE*>           files;
    std::vector<const TileInfo*> info;
    int  height = 0;
    bool ok     = false;

    bool open(const Grid& g, int tr) {
        close();
        files.assign(g.tile_cols, nullptr);
        info.assign(g.tile_cols, nullptr);
        height = 0;

        for (int tc = 0; tc < g.tile_cols; tc++) {
            const TileInfo* t = g.at(tr, tc);
            if (!t) {
                std::cerr << "[tiles] ERROR: missing tile tr=" << tr << " tc=" << tc << "\n";
                close(); return false;
            }
            files[tc] = fopen(t->path.c_str(), "rb");
            if (!files[tc]) {
                std::cerr << "[tiles] ERROR: cannot open " << t->path << "\n";
                close(); return false;
            }
            info[tc] = t;
            if (tc == 0) height = t->h;
            else if (t->h != height) {
                std::cerr << "[tiles] ERROR: ragged tile row tr=" << tr << "\n";
                close(); return false;
            }
        }
        ok = true;
        return true;
    }

    bool read_row(void* out, size_t elem_size) {
        if (!ok) return false;
        uint8_t* dst = (uint8_t*)out;
        for (size_t tc = 0; tc < files.size(); tc++) {
            const size_t w = (size_t)info[tc]->w;
            if (fread(dst + (size_t)info[tc]->c0 * elem_size,
                      elem_size, w, files[tc]) != w) {
                std::cerr << "[tiles] ERROR: short read on " << info[tc]->path << "\n";
                ok = false;
                return false;
            }
        }
        return true;
    }

    void close() {
        for (auto* f : files) if (f) fclose(f);
        files.clear(); info.clear();
        height = 0; ok = false;
    }

    ~RowReader() { close(); }
};

} // namespace tilegrid

#endif // TILEGRID_H_
