#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <functional>
#include <cstdint>
#include <random>
#include <algorithm>
#include <utility>
#include <cstdlib>
#include <png.h>

// --- structs as given by the user ---
struct FeasibleCell {
    int r;
    int c;
    uint8_t value;
    int overlap = 0;
};

struct FeasibleGrid {
    std::vector<FeasibleCell> cells;
    int origin_r;    // top-left of AABB in original matrix
    int origin_c;
    int height;      // AABB dimensions in cells
    int width;
};

class Renderer {
private:
    static bool writeColorPng(const std::string& path, int width, int height,
                               const std::function<void(int, uint8_t*)>& fillRow) {
        FILE* fp = fopen(path.c_str(), "wb");
        if (!fp) {
            std::cerr << "[overlap] FAIL: could not open " << path << " for writing\n";
            return false;
        }

        png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
        if (!png_ptr) {
            std::cerr << "[overlap] FAIL: png_create_write_struct failed\n";
            fclose(fp);
            return false;
        }

        png_infop info_ptr = png_create_info_struct(png_ptr);
        if (!info_ptr) {
            std::cerr << "[overlap] FAIL: png_create_info_struct failed\n";
            png_destroy_write_struct(&png_ptr, nullptr);
            fclose(fp);
            return false;
        }

        if (setjmp(png_jmpbuf(png_ptr))) {
            std::cerr << "[overlap] FAIL: libpng error while writing " << path << "\n";
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
                std::cout << "[overlap] Writing row " << r << "/" << height << "\n";
            fillRow(r, row.data());
            png_write_row(png_ptr, row.data());
        }

        png_write_end(png_ptr, info_ptr);
        png_destroy_write_struct(&png_ptr, &info_ptr);
        fclose(fp);
        return true;
    }

    // Same as writeColorPng but RGBA — lets non-feasible positions be fully
    // transparent (alpha=0) instead of needing a solid background fill.
    // fillRow gets width*4 bytes per row, in R,G,B,A order.
    static bool writeColorPngRGBA(const std::string& path, int width, int height,
                                   const std::function<void(int, uint8_t*)>& fillRow) {
        FILE* fp = fopen(path.c_str(), "wb");
        if (!fp) {
            std::cerr << "[overlap] FAIL: could not open " << path << " for writing\n";
            return false;
        }

        png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
        if (!png_ptr) {
            std::cerr << "[overlap] FAIL: png_create_write_struct failed\n";
            fclose(fp);
            return false;
        }

        png_infop info_ptr = png_create_info_struct(png_ptr);
        if (!info_ptr) {
            std::cerr << "[overlap] FAIL: png_create_info_struct failed\n";
            png_destroy_write_struct(&png_ptr, nullptr);
            fclose(fp);
            return false;
        }

        if (setjmp(png_jmpbuf(png_ptr))) {
            std::cerr << "[overlap] FAIL: libpng error while writing " << path << "\n";
            png_destroy_write_struct(&png_ptr, &info_ptr);
            fclose(fp);
            return false;
        }

        png_init_io(png_ptr, fp);
        png_set_IHDR(png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_RGB_ALPHA,
                     PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
        png_set_compression_level(png_ptr, 6);
        png_write_info(png_ptr, info_ptr);

        std::vector<uint8_t> row((size_t)width * 4);
        for (int r = 0; r < height; r++) {
            if (r % 2000 == 0 && r > 0)
                std::cout << "[overlap] Writing row " << r << "/" << height << "\n";
            fillRow(r, row.data());
            png_write_row(png_ptr, row.data());
        }

        png_write_end(png_ptr, info_ptr);
        png_destroy_write_struct(&png_ptr, &info_ptr);
        fclose(fp);
        return true;
    }
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
    static bool writeFeasibleGridOverlapPng(const FeasibleGrid& grid, const std::string& png_path) {
        if (grid.height <= 0 || grid.width <= 0) {
            std::cerr << "[overlap] FAIL: invalid grid dimensions " << grid.width
                      << " x " << grid.height << "\n";
            return false;
        }

        std::cout << "[overlap] Rendering " << grid.height << " x " << grid.width
                  << " grid (" << grid.cells.size() << " cells, origin "
                  << grid.origin_r << "," << grid.origin_c << ")\n";

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

    // --- new: read raw FeasibleCell array straight off disk ---
    static bool readFeasibleCellsBin(const std::string& bin_path, std::vector<FeasibleCell>& out_cells) {
        std::ifstream f(bin_path, std::ios::binary);
        if (!f.is_open()) {
            std::cerr << "[radar] FAIL: could not open " << bin_path << "\n";
            return false;
        }

        f.seekg(0, std::ios::end);
        size_t file_bytes = f.tellg();
        f.seekg(0, std::ios::beg);

        if (file_bytes % sizeof(FeasibleCell) != 0) {
            std::cerr << "[radar] FAIL: " << bin_path << " size (" << file_bytes
                      << " bytes) isn't a multiple of sizeof(FeasibleCell) ("
                      << sizeof(FeasibleCell) << ")\n";
            return false;
        }

        size_t count = file_bytes / sizeof(FeasibleCell);
        out_cells.resize(count);
        f.read(reinterpret_cast<char*>(out_cells.data()), file_bytes);
        if ((size_t)f.gcount() != file_bytes) {
            std::cerr << "[radar] FAIL: read " << f.gcount() << " of " << file_bytes << " bytes\n";
            return false;
        }

        std::cout << "[radar] Read " << count << " feasible cells from " << bin_path << "\n";
        return true;
    }

    // --- new: simulate radar deployment, derive overlap counts ---
    // Reads the feasible-cell bin at bin_path, picks num_radars distinct
    // feasible cells as stand-in radar sites (since no real coordinates were
    // given), and for every feasible cell counts how many radars are within
    // radar_range of it (Euclidean, in cell units). That count becomes the
    // cell's overlap value. Result lands in out_grid, ready to hand to
    // writeFeasibleGridOverlapPng().
    static bool simulateRadarDeployment(const std::string& bin_path,
                                         int num_radars,
                                         double radar_range,
                                         FeasibleGrid& out_grid,
                                         unsigned seed = 42) {
        std::vector<FeasibleCell> cells;
        if (!readFeasibleCellsBin(bin_path, cells)) return false;

        if (cells.empty()) {
            std::cerr << "[radar] FAIL: no feasible cells to deploy radars on\n";
            return false;
        }
        if (num_radars > (int)cells.size()) {
            std::cerr << "[radar] WARN: requested " << num_radars << " radars but only "
                      << cells.size() << " feasible cells exist — clamping\n";
            num_radars = (int)cells.size();
        }

        // reset — we're computing overlap fresh for this deployment
        for (auto& cell : cells) cell.overlap = 0;

        // bounding box, same convention writeFeasibleGridOverlapPng expects
        int min_r = cells[0].r, max_r = cells[0].r;
        int min_c = cells[0].c, max_c = cells[0].c;
        for (auto& cell : cells) {
            min_r = std::min(min_r, cell.r);
            max_r = std::max(max_r, cell.r);
            min_c = std::min(min_c, cell.c);
            max_c = std::max(max_c, cell.c);
        }

        // pick `num_radars` distinct feasible cells as radar sites — swap
        // this block for an explicit list of real coordinates if you have them
        std::mt19937 rng(seed);
        std::vector<int> indices(cells.size());
        for (size_t i = 0; i < indices.size(); i++) indices[i] = (int)i;
        std::shuffle(indices.begin(), indices.end(), rng);

        std::vector<std::pair<int, int>> radar_sites;
        for (int i = 0; i < num_radars; i++) {
            const auto& site = cells[indices[i]];
            radar_sites.push_back({site.r, site.c});
        }

        std::cout << "[radar] Deploying " << num_radars << " radars, range "
                  << radar_range << " cells\n";

        double range_sq = radar_range * radar_range;
        for (auto& site : radar_sites) {
            int radar_r = site.first, radar_c = site.second;
            for (auto& cell : cells) {
                double dr = cell.r - radar_r;
                double dc = cell.c - radar_c;
                if (dr * dr + dc * dc <= range_sq) {
                    cell.overlap++;
                }
            }
        }

        out_grid.cells    = std::move(cells);
        out_grid.origin_r = min_r;
        out_grid.origin_c = min_c;
        out_grid.height   = max_r - min_r + 1;
        out_grid.width    = max_c - min_c + 1;

        std::cout << "[radar] Done — overlap computed across "
                  << out_grid.cells.size() << " cells\n";
        return true;
    }
};

int main(int argc, char** argv) {
    // optional CLI overrides: ./test_radar [bin_path] [num_radars] [range] [seed]
    std::string bin_path   = (argc > 1) ? argv[1] : "/tmp/feasible_test.bin";
    int         num_radars = (argc > 2) ? std::atoi(argv[2]) : 4;
    double      radar_range = (argc > 3) ? std::atof(argv[3]) : 4.0;
    unsigned    seed        = (argc > 4) ? (unsigned)std::atoi(argv[4]) : 42;

    // if no real feasible.bin path was given, synthesize one so this is
    // runnable standalone — a 20x20 block of feasible cells
    if (argc <= 1) {
        std::cout << "[test] No bin path given — writing a synthetic 20x20 feasible.bin to "
                  << bin_path << "\n";
        std::vector<FeasibleCell> cells;
        for (int r = 0; r < 20; r++)
            for (int c = 0; c < 20; c++)
                cells.push_back({r, c, /*value=*/1, /*overlap=*/0});

        std::ofstream f(bin_path, std::ios::binary);
        if (!f) {
            std::cerr << "[test] FAIL: could not write synthetic bin to " << bin_path << "\n";
            return 1;
        }
        f.write(reinterpret_cast<char*>(cells.data()), cells.size() * sizeof(FeasibleCell));
    }

    // --- test 1: readFeasibleCellsBin (called internally) + simulateRadarDeployment ---
    FeasibleGrid grid;
    if (!Renderer::simulateRadarDeployment(bin_path, num_radars, radar_range, grid, seed)) {
        std::cerr << "[test] simulateRadarDeployment FAILED\n";
        return 1;
    }

    // quick sanity readout: how many cells fall in each overlap bucket
    int buckets[6] = {0, 0, 0, 0, 0, 0};   // 0,1,2,3,4,5+
    for (const auto& cell : grid.cells) {
        int b = std::min(cell.overlap, 5);
        buckets[b]++;
    }
    std::cout << "[test] Overlap distribution across " << grid.cells.size() << " cells:\n";
    for (int i = 0; i < 6; i++)
        std::cout << "  overlap " << i << (i == 5 ? "+" : "") << " : " << buckets[i] << " cells\n";

    // --- test 2: writeFeasibleGridOverlapPng ---
    std::string png_path = "/tmp/radar_overlap_test.png";
    if (!Renderer::writeFeasibleGridOverlapPng(grid, png_path)) {
        std::cerr << "[test] writeFeasibleGridOverlapPng FAILED\n";
        return 1;
    }

    std::cout << "[test] PASS — open " << png_path << " to inspect the result\n";
    return 0;
}
