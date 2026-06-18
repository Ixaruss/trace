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
    public: bool flag_water_area = true;
    public: bool flag_ib = true;
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

    static FeasibleGrid getFeasibleGrid(short op = IN_MEM, short opt = 0);

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
        private:
            // overlap -> color ramp: 0=white, 2=yellow, 3=orange, 4=red, 5+=dark red.
            // You didn't specify overlap==1, so it's linearly interpolated between
            // the 0 and 2 stops (pale yellow) rather than left undefined — easy to
            // change if you had something else in mind for it.
            static void overlapToColor(int overlap, uint8_t& red, uint8_t& green, uint8_t& blue) {
                if (overlap < 0) overlap = 0;
                switch (overlap) {
                    case 0: red = 255; green = 255; blue = 255; break;  // white
                    case 1: red = 0;   green = 200; blue = 0;   break;  // green
                    case 2: red = 255; green = 255; blue = 0;   break;  // yellow
                    case 3: red = 255; green = 165; blue = 0;   break;  // orange
                    case 4: red = 255; green = 0;   blue = 0;   break;  // red
                    default:                                            // 5 and above
                        red = 139; green = 0; blue = 0;
                        break;
                }
            }

        public:
            // Renders a FeasibleGrid's overlap values into a PNG sized grid.width x
            // grid.height. Cell coordinates are assumed to be in original-matrix
            // space (per your origin_r/origin_c comment), so each cell is placed at
            // (cell.r - origin_r, cell.c - origin_c) on the canvas.
            static bool writeFeasibleGridOverlapPng(const FeasibleGrid& grid, const std::string& png_path) {
                if (grid.height <= 0 || grid.width <= 0) {
                    std::cerr << "[overlap] FAIL: invalid grid dimensions " << grid.width
                              << " x " << grid.height << "\n";
                    return false;
                }

                std::cout << "[overlap] Rendering " << grid.height << " x " << grid.width
                          << " grid (" << grid.cells.size() << " cells, origin "
                          << grid.origin_r << "," << grid.origin_c << ")\n";

                // background for any cell position with no entry in grid.cells — a
                // distinct dark gray so it's never confused with an actual
                // overlap==0 (white) cell. Change this if your cells vector is
                // already dense over the full AABB and this case can't occur.
                std::vector<uint8_t> canvas((size_t)grid.height * grid.width * 3, 30);

                int skipped = 0;
                for (const auto& cell : grid.cells) {
                    int local_r = cell.r - grid.origin_r;
                    int local_c = cell.c - grid.origin_c;
                    if (local_r < 0 || local_r >= grid.height || local_c < 0 || local_c >= grid.width) {
                        skipped++;
                        continue;
                    }
                    uint8_t red, green, blue;
                    overlapToColor(cell.overlap, red, green, blue);

                    size_t idx = ((size_t)local_r * grid.width + local_c) * 3;
                    canvas[idx + 0] = red;
                    canvas[idx + 1] = green;
                    canvas[idx + 2] = blue;
                }
                if (skipped > 0)
                    std::cerr << "[overlap] WARN: " << skipped << " cells fell outside the "
                              << grid.width << "x" << grid.height << " AABB and were skipped\n";

                bool ok = writeColorPng(png_path, grid.width, grid.height,
                                         [&canvas, w = grid.width](int r, uint8_t* rgbRow) {
                                             std::memcpy(rgbRow, &canvas[(size_t)r * w * 3], (size_t)w * 3);
                                         });
                if (!ok) return false;

                std::cout << "[overlap] PNG written → " << png_path << "\n";
                return true;
            }
            public:
                        // Plots the raw (r,c) of each point as a red dot — presence only,
                        // no coverage/overlap rendering. Canvas is sized to fit whatever
                        // coordinates are passed in (same bounding-box approach as the
                        // FeasibleCell point-cloud branch in verifyMatrixBin), so it works
                        // for a single chromosome's radars or the whole population's
                        // combined radars without you having to specify dimensions.
                        // PointT just needs public `.r`/`.c` int members — GridPoint,
                        // FeasibleCell, etc. all qualify.
                        template <typename PointT>
                        static bool writeRadarScatterPng(const std::vector<PointT>& points, const std::string& png_path) {
                            if (points.empty()) {
                                std::cerr << "[radar-viz] FAIL: no points to plot for " << png_path << "\n";
                                return false;
                            }

                            int min_r = points[0].r, max_r = points[0].r;
                            int min_c = points[0].c, max_c = points[0].c;
                            for (const auto& p : points) {
                                if (p.r < min_r) min_r = p.r;
                                if (p.r > max_r) max_r = p.r;
                                if (p.c < min_c) min_c = p.c;
                                if (p.c > max_c) max_c = p.c;
                            }

                            int img_rows = max_r - min_r + 1;
                            int img_cols = max_c - min_c + 1;

                            std::vector<uint8_t> canvas((size_t)img_rows * img_cols * 3, 0);  // black background
                            for (const auto& p : points) {
                                size_t idx = ((size_t)(p.r - min_r) * img_cols + (p.c - min_c)) * 3;
                                canvas[idx + 0] = 255;  // red — radar present
                                canvas[idx + 1] = 0;
                                canvas[idx + 2] = 0;
                            }

                            bool ok = writeColorPng(png_path, img_cols, img_rows,
                                                     [&canvas, img_cols](int r, uint8_t* rgbRow) {
                                                         std::memcpy(rgbRow, &canvas[(size_t)r * img_cols * 3], (size_t)img_cols * 3);
                                                     });
                            if (!ok) return false;

                            std::cout << "[radar-viz] PNG written → " << png_path << " (" << points.size()
                                      << " points, " << img_cols << "x" << img_rows << ")\n";
                            return true;
                        }

};

#endif
