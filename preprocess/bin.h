#ifndef BIN_H_
#define BIN_H_


#include <functional>
#include <gdal_priv.h>
#include <ogr_spatialref.h>
#include <iostream>
#include <fstream>
#include <optional>
#include <vector>
#include <cstdint>
#include <cmath>
#include <string>
#include "./lib.h"
#include <png.h>
// Global constants
const bool WRITE = 1;
const bool IN_MEM = 0;
const bool USE_MULTI_THREADS = 1;
const bool FULL_MATRIX = 1;
const bool LIMIT_MEMORY = 1;
const bool COMBINE = 0;
const bool WITHALL = 1;

struct BoundingBox {
    double min_lat;
    double max_lat;
    double min_lon;
    double max_lon;
};

#pragma pack(1)
struct ElevCell {
    int16_t elevation;
    float   slope;
};
#pragma pack()

struct FeasibleCell {
    int r;
    int c;
    uint8_t value;
    int overlap =0;
};
struct GridPoint {
    int r;
    int c;
    int typeIndex;
};
enum LoadOptions {
    ALL,
    CONSTRAINTS,
    ELEVATION,
    SLOPE_
};

struct FeasibleGrid {
    std::vector<FeasibleCell> cells;
    int origin_r;    // top-left of AABB in original matrix
    int origin_c;
    int height;      // jAABB dimensions in cells
    int width;
};
class Bin {
    // In-memory data caches
    public: std::vector<uint8_t> vec_matrix;
    public: std::vector<int16_t> elev_matrix;
    public: std::vector<float> slope_matrix;
    public: Config conf;


    // Loading tracking flags
    public: bool vector_loaded = false;
    public: bool elevation_loaded = false;
    public: bool slope_loaded = false;
    public: int VEC_ROWS = 0;
    public: int VEC_COLS = 0;

    public: bool flag_road = false;
    public: bool flag_rail = false;
    public: bool flag_water_line = true;
    public: bool flag_water_area = false;
    public: bool flag_ib = false;
    public: bool flag_slope = false;

    public: BoundingBox mbr;


    OGRCoordinateTransformation* poCT;

public:
    Bin() : poCT(nullptr) {
        GDALAllRegister();
        if(Config::is_available())
            conf = Config::load();
        // Instantiate spatial converter from standard Lat/Lng (WGS84) to your local project grid (EPSG:7755)
        OGRSpatialReference oSourceSRS, oTargetSRS;
        oSourceSRS.importFromEPSG(4326);
        oTargetSRS.importFromEPSG(7755);

        // Enforce traditional spatial mapping (X = Longitude, Y = Latitude) regardless of system GDAL version
        oSourceSRS.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
        oTargetSRS.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

        poCT = OGRCreateCoordinateTransformation(&oSourceSRS, &oTargetSRS);
        if (!poCT) {
            std::cerr << "[init] CRITICAL ERROR: Could not initialize OGRCoordinateTransformation!\n";
        }
    }

    ~Bin() {
        if (poCT) OGRCoordinateTransformation::DestroyCT(poCT);
    }
    public: void init(std::set<LoadOptions> Options = {});

    public: static int build_matrix();

    public: static int build_elevation_matrix(bool split = 1);

    public: static int refresh_road_layer(std::optional<std::string> path);

    public: static int refresh_rail_layer(std::optional<std::string> path);

    public: static int refresh_water_area_layer(std::optional<std::string> path);

    public: static int refresh_water_line_layer(std::optional<std::string> path);

    public: static int refresh_ib_layer(std::optional<std::string> path = std::nullopt);

    public: static bool point_in_polygon(double px, double py, const std::vector<std::pair<double,double>>& proj_poly);

    static void write_heatmap_bin_h5 ();

    bool getGridCoords(double lat, double lng, int& out_r, int& out_c);

    int16_t getElevation(double lat, double lng);

    bool getSlope(double lat, double lng);

    bool getStatusRoad(double lat, double lng);

    bool getStatusRail(double lat, double lng);

    bool getStatusIB(double lat, double lng);

    bool getStatusWaterArea(double lat, double lng);

    bool getStatusWaterLine(double lat, double lng);

    bool latlon_to_projected(double lat, double lon, double& px, double& py);

    std::vector<int> getStatusAll(double lat, double lng);

    BoundingBox getMinimumBoundingBox();

    static void getonlyFeasibleGrid(Bin& b);

    static FeasibleGrid getFeasibleGrid(Bin& b, short op = IN_MEM, short opt = 0, short allCells = 0);

    std::pair<double,double> getMapOrigins();

    int getTerrianResolution();

    bool lineOfSight(int x1, int y1, int x2, int y2);

    void visualize(const std::string& output_path);

    static void get_elevation_color(int16_t alt, uint8_t& r, uint8_t& g, uint8_t& b);

    private:
        // ------------------------------------------------------------------
        // Generic row-streaming PNG writer (RGB8). Calls fillRow(r, rgbRow)
        // once per row; rgbRow points to width*3 bytes to fill in
        // R,G,B,R,G,B,... order. Same memory profile as the old PPM loop —
        // only one row buffer lives at a time, not the whole image.
        // ------------------------------------------------------------------
        static bool writeColorPng(const std::string& path, int width, int height,
                                   const std::function<void(int, uint8_t*)>& fillRow) {
            FILE* fp = fopen(path.c_str(), "wb");
            if (!fp) {
                std::cerr << "[verify] FAIL: could not open " << path << " for writing\n";
                return false;
            }

            png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
            if (!png_ptr) {
                std::cerr << "[verify] FAIL: png_create_write_struct failed\n";
                fclose(fp);
                return false;
            }

            png_infop info_ptr = png_create_info_struct(png_ptr);
            if (!info_ptr) {
                std::cerr << "[verify] FAIL: png_create_info_struct failed\n";
                png_destroy_write_struct(&png_ptr, nullptr);
                fclose(fp);
                return false;
            }

            if (setjmp(png_jmpbuf(png_ptr))) {
                std::cerr << "[verify] FAIL: libpng error while writing " << path << "\n";
                png_destroy_write_struct(&png_ptr, &info_ptr);
                fclose(fp);
                return false;
            }

            png_init_io(png_ptr, fp);
            png_set_IHDR(png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_RGB,
                         PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
            png_set_compression_level(png_ptr, 6);
            png_write_info(png_ptr, info_ptr);

            std::vector<uint8_t> row((size_t)width * 3);
            for (int r = 0; r < height; r++) {
                if (r % 2000 == 0 && r > 0)
                    std::cout << "[verify] Writing row " << r << "/" << height << "\n";
                fillRow(r, row.data());
                png_write_row(png_ptr, row.data());
            }

            png_write_end(png_ptr, info_ptr);
            png_destroy_write_struct(&png_ptr, &info_ptr);
            fclose(fp);
            return true;
        }

    public:
        template <typename T>
        static bool verifyMatrixBin(const std::string& path, int rows, int cols) {
            std::cout << "[verify] Checking: " << path << "\n";
            std::cout << "[verify] Expected: " << rows << " x " << cols
                      << " (" << (size_t)rows * cols << " cells, "
                      << ((size_t)rows * cols * sizeof(T)) / 1e6 << " MB, "
                      << sizeof(T) << " bytes/cell)\n";

            // --- 1. check file exists ---
            std::ifstream f(path, std::ios::binary);
            if (!f.is_open()) {
                std::cerr << "[verify] FAIL: could not open " << path << "\n";
                return false;
            }

            // --- 2. check file size ---
            f.seekg(0, std::ios::end);
            size_t file_bytes     = f.tellg();
            size_t expected_bytes = (size_t)rows * cols * sizeof(T);
            f.seekg(0, std::ios::beg);

            if (file_bytes != expected_bytes) {
                std::cerr << "[verify] FAIL: size mismatch — expected "
                          << expected_bytes / 1e6 << " MB, got "
                          << file_bytes / 1e6 << " MB\n";
                return false;
            }
            std::cout << "[verify] Size OK: " << file_bytes / 1e6 << " MB\n";

            // --- 3. read ---
            std::vector<T> matrix((size_t)rows * cols);
            f.read(reinterpret_cast<char*>(matrix.data()), expected_bytes);
            if ((size_t)f.gcount() != expected_bytes) {
                std::cerr << "[verify] FAIL: read " << f.gcount()
                          << " of " << expected_bytes << " bytes\n";
                return false;
            }
            f.close();

            // derive png path from bin path
            std::string png_path = path.substr(0, path.find_last_of('.')) + "_verify.png";

            // --- 4. stats + png ---
            if constexpr (std::is_same<T, uint8_t>::value) {
                // ---- uint8 (vector/bitmask matrix) ----
                T      min_val = matrix[0];
                T      max_val = matrix[0];
                size_t nonzero = 0;

                for (size_t i = 0; i < matrix.size(); i++) {
                    if (matrix[i] < min_val) min_val = matrix[i];
                    if (matrix[i] > max_val) max_val = matrix[i];
                    if (matrix[i] != 0)      nonzero++;
                }
                double coverage = (nonzero / (double)matrix.size()) * 100.0;
                std::cout << "[verify] Min value : " << (int)min_val  << "\n";
                std::cout << "[verify] Max value : " << (int)max_val  << "\n";
                std::cout << "[verify] Non-zero  : " << nonzero << " cells (" << coverage << "%)\n";
                if (nonzero == 0)
                    std::cerr << "[verify] WARN: matrix is entirely zero — likely a failed build\n";

                // png — same color scheme as proof.ppm
                bool ok = writeColorPng(png_path, cols, rows, [&matrix, cols](int r, uint8_t* rgbRow) {
                    for (int c = 0; c < cols; c++) {
                        uint8_t val = matrix[(size_t)r * cols + c];
                        uint8_t red = 20, green = 24, blue = 38;   // dark background
                        if (val & (1 << 0)) { red = 255; green = 50;  blue = 0;   }  // roads   → orange
                        if (val & (1 << 1)) { red = 255; green = 255; blue = 0;   }  // railway → yellow
                        if (val & (1 << 2)) { red = 0;   green = 100; blue = 255; }  // waterline → blue
                        if (val & (1 << 3)) { red = 0;   green = 180; blue = 255; }  // waterbody → light blue
                        if ((val & (1 << 0)) && (val & (1 << 1))) {
                            red = 0; green = 255; blue = 255;   // road+rail crossing → cyan
                        }
                        rgbRow[c * 3 + 0] = red;
                        rgbRow[c * 3 + 1] = green;
                        rgbRow[c * 3 + 2] = blue;
                    }
                });
                if (!ok) return false;

            } else if constexpr (std::is_same<T, int16_t>::value) {
                // ---- int16 (elevation matrix) ----
                T      min_val = matrix[0];
                T      max_val = matrix[0];
                size_t nonzero = 0;

                for (size_t i = 0; i < matrix.size(); i++) {
                    if (matrix[i] < min_val) min_val = matrix[i];
                    if (matrix[i] > max_val) max_val = matrix[i];
                    if (matrix[i] != 0)      nonzero++;
                }
                double coverage = (nonzero / (double)matrix.size()) * 100.0;
                std::cout << "[verify] Min value : " << min_val << "m\n";
                std::cout << "[verify] Max value : " << max_val << "m\n";
                std::cout << "[verify] Non-zero  : " << nonzero << " cells (" << coverage << "%)\n";
                if (nonzero == 0)
                    std::cerr << "[verify] WARN: matrix is entirely zero — likely a failed build\n";

                // png — terrain color scheme
                bool ok = writeColorPng(png_path, cols, rows, [&matrix, cols](int r, uint8_t* rgbRow) {
                    for (int c = 0; c < cols; c++) {
                        uint8_t val = matrix[(size_t)r * cols + c];
                        uint8_t red = 20, green = 24, blue = 38;

                        if (val & (1 << BIT_ROAD))       { red = 255; green = 50;  blue = 0;   }
                        if (val & (1 << BIT_RAIL))       { red = 255; green = 255; blue = 0;   }
                        if (val & (1 << BIT_WATER_LINE)) { red = 0;   green = 100; blue = 255; }
                        if (val & (1 << BIT_WATER_AREA)) { red = 0;   green = 180; blue = 255; }
                        if (val & (1 << BIT_IB))         { red = 80;  green = 0;   blue = 80;  }

                        // intersections
                        if ((val & (1 << BIT_ROAD)) && (val & (1 << BIT_RAIL)))
                            { red = 0; green = 255; blue = 255; }

                        rgbRow[c * 3 + 0] = red;
                        rgbRow[c * 3 + 1] = green;
                        rgbRow[c * 3 + 2] = blue;
                    }
                });
                if (!ok) return false;

            } else {
                // ---- struct (FeasibleCell etc.) ----
                // treat as a sparse point cloud — plot each cell's (r,c) on a blank grid
                std::cout << "[verify] Cell count: " << matrix.size() << "\n";
                std::cout << "[verify] Cell size : " << sizeof(T) << " bytes (struct)\n";

                int min_r = INT_MAX, min_c = INT_MAX;
                int max_r = 0,       max_c = 0;
                for (auto& cell : matrix) {
                    if (cell.r < min_r) min_r = cell.r;
                    if (cell.c < min_c) min_c = cell.c;
                    if (cell.r > max_r) max_r = cell.r;
                    if (cell.c > max_c) max_c = cell.c;
                }

                int img_rows = max_r - min_r + 1;
                int img_cols = max_c - min_c + 1;
                std::cout << "[verify] Reconstructing " << img_rows << " x " << img_cols << " grid from cell coords\n";

                // blank dark canvas
                std::vector<uint8_t> canvas((size_t)img_rows * img_cols * 3, 0);
                for (auto& cell : matrix) {
                    size_t idx = ((size_t)(cell.r - min_r) * img_cols + (cell.c - min_c)) * 3;
                    canvas[idx + 0] = 255;  // white dot per feasible cell
                    canvas[idx + 1] = 255;
                    canvas[idx + 2] = 255;
                }

                bool ok = writeColorPng(png_path, img_cols, img_rows, [&canvas, img_cols](int r, uint8_t* rgbRow) {
                    std::memcpy(rgbRow, &canvas[(size_t)r * img_cols * 3], (size_t)img_cols * 3);
                });
                if (!ok) return false;
            }

            std::cout << "[verify] PNG written → " << png_path << "\n";
            std::cout << "[verify] PASS: " << path << "\n";
            return true;
        }
        static void overlapToColor(int overlap, uint8_t& red, uint8_t& green, uint8_t& blue) {
                        if (overlap < 0) overlap = 0;
                        switch (overlap) {
                            case 0: red = 255; green = 255; blue = 255; break;  // white — unused directly now, kept for callers
                            case 1: red = 0;   green = 200; blue = 0;   break;  // green
                            case 2: red = 255; green = 255; blue = 0;   break;  // yellow
                            case 3: red = 255; green = 165; blue = 0;   break;  // orange
                            case 4: red = 255; green = 0;   blue = 0;   break;  // red
                            default:                                            // 5 and above
                                red = 139; green = 0; blue = 0;
                                break;
                        }
                    }

        // Terrain colour for a raw cell.value byte — same palette as the vec_matrix verify PNG.
        // Only bits whose constraint flag is actually ON get colored; a disabled constraint
        // (flag_* == false) is treated as visually absent even if the physical feature exists,
        // since it isn't restricting radar placement anymore.
        static void terrainToColor(uint8_t val, uint8_t& red, uint8_t& green, uint8_t& blue,
                                    bool flag_road       = true,
                                    bool flag_rail       = true,
                                    bool flag_water_line = true,
                                    bool flag_water_area = true,
                                    bool flag_ib         = true,
                                    bool flag_slope      = true) {
            red = 255; green = 255; blue = 255; // default: plain AOI, no active feature
            if (flag_slope      && (val & (1 << BIT_SLOPE)))      { red = 139; green = 69;  blue = 19;  } // brown
            if (flag_ib         && (val & (1 << BIT_IB)))         { red = 80;  green = 0;   blue = 80;  } // purple
            if (flag_water_area && (val & (1 << BIT_WATER_AREA))) { red = 0;   green = 180; blue = 255; } // light blue
            if (flag_water_line && (val & (1 << BIT_WATER_LINE))) { red = 0;   green = 100; blue = 255; } // blue
            if (flag_rail       && (val & (1 << BIT_RAIL)))       { red = 255; green = 255; blue = 0;   } // yellow
            if (flag_road       && (val & (1 << BIT_ROAD)))       { red = 255; green = 50;  blue = 0;   } // orange
            if (flag_road && flag_rail &&
                (val & (1 << BIT_ROAD)) && (val & (1 << BIT_RAIL)))
                                                                    { red = 0;   green = 255; blue = 255; } // cyan — road/rail crossing
        }

        public:
            // aoiGrid:      ALL cells in the Area of Operation (Bin::getFeasibleGrid(..., allCells=WITHALL)),
            //               with `overlap` already computed on it — this is what actually gets rendered.
            // feasibleGrid: the filtered subset (road/rail-adjacent, water/IB/slope-free) that the GA is
            //               allowed to place radars on — drawn as a light-blue tint so it stays visible on
            //               top of the terrain colours without hiding them, but never overrides coverage.
            // flag_*:       which constraints are actually active (mirrors Bin's flag_road/flag_rail/etc).
            //               Disabled constraints are rendered as plain AOI even if the physical feature exists.
            static bool writeFeasibleGridOverlapPng(
                const FeasibleGrid&            aoiGrid,
                const std::string&             png_path,
                const std::vector<GridPoint>&  placements   = {},
                const FeasibleGrid&            feasibleGrid = {},
                bool                            flag_road       = true,
                bool                            flag_rail       = true,
                bool                            flag_water_line = true,
                bool                            flag_water_area = true,
                bool                            flag_ib         = true,
                bool                            flag_slope      = true)
            {
                if (aoiGrid.height <= 0 || aoiGrid.width <= 0) {
                    std::cerr << "[overlap] FAIL: invalid grid dimensions "
                              << aoiGrid.width << " x " << aoiGrid.height << "\n";
                    return false;
                }
                std::cout << "[overlap] Rendering " << aoiGrid.height << " x " << aoiGrid.width
                          << " grid (" << aoiGrid.cells.size() << " AOI cells, origin "
                          << aoiGrid.origin_r << "," << aoiGrid.origin_c << ")\n";

                // Dark-grey background — anything never touched (outside AOI / outside AABB)
                std::vector<uint8_t> canvas((size_t)aoiGrid.height * aoiGrid.width * 3, 30);

                // ── 0. Paint the full AOI: terrain colour, or coverage colour where overlap > 0 ──
                // Radar signal doesn't care about water/slope/IB — coverage always wins when present.
                int skipped = 0;
                for (const auto& cell : aoiGrid.cells) {
                    int local_r = cell.r - aoiGrid.origin_r;
                    int local_c = cell.c - aoiGrid.origin_c;
                    if (local_r < 0 || local_r >= aoiGrid.height ||
                        local_c < 0 || local_c >= aoiGrid.width) {
                        ++skipped;
                        continue;
                    }
                    uint8_t red, green, blue;
                    if (cell.overlap > 0) overlapToColor(cell.overlap, red, green, blue);
                    else                  terrainToColor(cell.value, red, green, blue,
                                                          flag_road, flag_rail,
                                                          flag_water_line, flag_water_area,
                                                          flag_ib, flag_slope);

                    size_t idx = ((size_t)local_r * aoiGrid.width + local_c) * 3;
                    canvas[idx + 0] = red;
                    canvas[idx + 1] = green;
                    canvas[idx + 2] = blue;
                }
                if (skipped > 0)
                    std::cerr << "[overlap] WARN: " << skipped
                              << " AOI cells fell outside the "

                              << aoiGrid.width << "x" << aoiGrid.height
                              << " AABB and were skipped\n";

                // ── 1. Tint feasible (radar-placeable) cells light blue, blended with ─────
                //      whatever's already there, so terrain/coverage colour still shows through.
                int feas_painted = 0, feas_skipped = 0;
                constexpr uint8_t TR = 135, TG = 206, TB = 250; // light blue tint
                for (const auto& cell : feasibleGrid.cells) {
                    int local_r = cell.r - aoiGrid.origin_r;
                    int local_c = cell.c - aoiGrid.origin_c;
                    if (local_r < 0 || local_r >= aoiGrid.height ||
                        local_c < 0 || local_c >= aoiGrid.width) {
                        ++feas_skipped;
                        continue;
                    }
                    size_t idx = ((size_t)local_r * aoiGrid.width + local_c) * 3;
                    canvas[idx + 0] = (uint8_t)((canvas[idx + 0] + TR) / 2);
                    canvas[idx + 1] = (uint8_t)((canvas[idx + 1] + TG) / 2);
                    canvas[idx + 2] = (uint8_t)((canvas[idx + 2] + TB) / 2);
                    ++feas_painted;
                }
                if (!feasibleGrid.cells.empty())
                    std::cout << "[overlap] Feasible overlay: " << feas_painted
                              << " cells tinted, " << feas_skipped << " outside AABB\n";

                // ── 2. Draw radar centre markers (drawn last → always on top) ────────
                if (!placements.empty()) {
                    // Safe single-pixel writer — clips to canvas bounds silently
                    auto px = [&](int lr, int lc,
                                  uint8_t r, uint8_t g, uint8_t b) __attribute__((always_inline))
                    {
                        if ((unsigned)lr >= (unsigned)aoiGrid.height ||
                            (unsigned)lc >= (unsigned)aoiGrid.width) return;
                        size_t idx = ((size_t)lr * aoiGrid.width + lc) * 3;
                        canvas[idx + 0] = r;
                        canvas[idx + 1] = g;
                        canvas[idx + 2] = b;
                    };

                    constexpr int     ARM    = 5;
                    constexpr uint8_t MR = 255, MG = 0, MB = 255;   // magenta arms
                    constexpr uint8_t BR = 0,   BG = 0, BB = 0;     // black border / dot

                    for (const auto& p : placements) {
                        const int cr = p.r - aoiGrid.origin_r;
                        const int cc = p.c - aoiGrid.origin_c;
                        // Skip placements whose centre maps outside the canvas
                        if (cr < 0 || cr >= aoiGrid.height ||
                            cc < 0 || cc >= aoiGrid.width) continue;

                        // Pass 1 — black halo (1-px border around every arm pixel)
                        //   Vertical halo: columns cc-1 and cc+1, rows cr±(ARM+1)
                        //   Horizontal halo: rows cr-1 and cr+1, cols cc±(ARM+1)
                        for (int d = -(ARM + 1); d <= (ARM + 1); ++d) {
                            px(cr + d, cc - 1, BR, BG, BB);
                            px(cr + d, cc + 1, BR, BG, BB);
                            px(cr - 1, cc + d, BR, BG, BB);
                            px(cr + 1, cc + d, BR, BG, BB);
                        }
                        // End caps (tips of each arm)
                        px(cr - (ARM + 1), cc, BR, BG, BB);
                        px(cr + (ARM + 1), cc, BR, BG, BB);
                        px(cr, cc - (ARM + 1), BR, BG, BB);
                        px(cr, cc + (ARM + 1), BR, BG, BB);

                        // Pass 2 — magenta crosshair arms
                        for (int d = -ARM; d <= ARM; ++d) {
                            px(cr + d, cc, MR, MG, MB);
                            px(cr, cc + d, MR, MG, MB);
                        }

                        // Pass 3 — 3×3 black centre dot (marks the exact cell)
                        for (int dr = -1; dr <= 1; ++dr)
                            for (int dc = -1; dc <= 1; ++dc)
                                px(cr + dr, cc + dc, BR, BG, BB);
                    }

                    std::cout << "[overlap] " << placements.size()
                              << " radar marker(s) drawn\n";
                }

                // ── 3. Write PNG ─────────────────────────────────────────────────────
                bool ok = writeColorPng(
                    png_path, aoiGrid.width, aoiGrid.height,
                    [&canvas, w = aoiGrid.width](int r, uint8_t* rgbRow) {
                        std::memcpy(rgbRow, &canvas[(size_t)r * w * 3], (size_t)w * 3);
                    });
                if (!ok) return false;

                std::cout << "[overlap] PNG written → " << png_path << "\n";
                return true;
            }

                        static bool writeReverseBreakPointsPng(
                            const std::string&            png_path,
                            int                            originR, int originC,   // top-left of bounding box, grid coords
                            int                            width,   int height,     // bounding box size (pixels == cells)
                            int                            srcRow,  int srcCol,     // absolute grid coords of the source
                            const std::vector<GridPoint>& breakPoints)              // absolute grid coords
                        {
                            if (width <= 0 || height <= 0) {
                                std::cerr << "[reverse-coverage] FAIL: invalid canvas dimensions "
                                          << width << " x " << height << "\n";
                                return false;
                            }

                            std::vector<uint8_t> canvas((size_t)height * width * 3, 30); // dark grey background

                            // Safe single-pixel writer — clips to canvas bounds silently
                            auto px = [&](int lr, int lc, uint8_t r, uint8_t g, uint8_t b) {
                                if ((unsigned)lr >= (unsigned)height || (unsigned)lc >= (unsigned)width) return;
                                size_t idx = ((size_t)lr * width + lc) * 3;
                                canvas[idx + 0] = r;
                                canvas[idx + 1] = g;
                                canvas[idx + 2] = b;
                            };

                            // Bresenham line in local (canvas) coordinates
                            auto drawLine = [&](int r0, int c0, int r1, int c1, uint8_t r, uint8_t g, uint8_t b) {
                                const int adr = std::abs(r1 - r0), adc = std::abs(c1 - c0);
                                const int sr  = (r0 < r1) ? 1 : -1;
                                const int sc  = (c0 < c1) ? 1 : -1;
                                int err = adr - adc;
                                int row = r0, col = c0;
                                while (true) {
                                    px(row, col, r, g, b);
                                    if (row == r1 && col == c1) break;
                                    const int e2 = 2 * err;
                                    if (e2 > -adc) { err -= adc; row += sr; }
                                    if (e2 <  adr) { err += adr; col += sc; }
                                }
                            };

                            const int localSrcR = srcRow - originR;
                            const int localSrcC = srcCol - originC;

                            // ── 1. Covered segments — green line source -> each break point ──────
                            for (const auto& bp : breakPoints) {
                                drawLine(localSrcR, localSrcC,
                                          bp.r - originR, bp.c - originC,
                                          60, 200, 80);
                            }

                            // ── 2. Break-point markers — small red dots, drawn on top of lines ────
                            for (const auto& bp : breakPoints) {
                                const int lr = bp.r - originR;
                                const int lc = bp.c - originC;
                                for (int dr = -1; dr <= 1; ++dr)
                                    for (int dc = -1; dc <= 1; ++dc)
                                        px(lr + dr, lc + dc, 220, 30, 30);
                            }

                            // ── 3. Source marker — magenta crosshair + black halo (drawn last) ────
                            constexpr int     ARM = 5;
                            constexpr uint8_t MR = 255, MG = 0, MB = 255;   // magenta arms
                            constexpr uint8_t BR = 0,   BG = 0, BB = 0;     // black border / dot

                            for (int d = -(ARM + 1); d <= (ARM + 1); ++d) {
                                px(localSrcR + d, localSrcC - 1, BR, BG, BB);
                                px(localSrcR + d, localSrcC + 1, BR, BG, BB);
                                px(localSrcR - 1, localSrcC + d, BR, BG, BB);
                                px(localSrcR + 1, localSrcC + d, BR, BG, BB);
                            }
                            px(localSrcR - (ARM + 1), localSrcC, BR, BG, BB);
                            px(localSrcR + (ARM + 1), localSrcC, BR, BG, BB);
                            px(localSrcR, localSrcC - (ARM + 1), BR, BG, BB);
                            px(localSrcR, localSrcC + (ARM + 1), BR, BG, BB);

                            for (int d = -ARM; d <= ARM; ++d) {
                                px(localSrcR + d, localSrcC, MR, MG, MB);
                                px(localSrcR, localSrcC + d, MR, MG, MB);
                            }

                            for (int dr = -1; dr <= 1; ++dr)
                                for (int dc = -1; dc <= 1; ++dc)
                                    px(localSrcR + dr, localSrcC + dc, BR, BG, BB);

                            // ── 4. Write PNG ───────────────────────────────────────────────────────
                            bool ok = writeColorPng(
                                png_path, width, height,
                                [&canvas, w = width](int r, uint8_t* rgbRow) {
                                    std::memcpy(rgbRow, &canvas[(size_t)r * w * 3], (size_t)w * 3);
                                });
                            if (!ok) return false;

                            std::cout << "[reverse-coverage] PNG written -> " << png_path << "\n";
                            return true;
                        }


};

#endif
