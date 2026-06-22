// ================================================================================
// STEP 1: GENERATE LOGS DIRECTLY INTO HDF5 CONTAINER
// ================================================================================
//
// Converted from Python (h5py / NumPy) to C++17 + HDF5 C++ API (H5Cpp).
//
// Build (GCC / Clang):
//   g++ -std=c++17 -O2 generate_hdf5.cpp -o generate_hdf5 -lhdf5_cpp -lhdf5
//
// Install HDF5 dev headers:
//   Ubuntu/Debian : sudo apt  install libhdf5-dev
//   Fedora/RHEL   : sudo dnf  install hdf5-devel
//   macOS (brew)  : brew install hdf5
//
// CMake one-liner:
//   find_package(HDF5 REQUIRED COMPONENTS CXX)
//   target_link_libraries(generate_hdf5 PRIVATE HDF5::HDF5)
//
// NOTE ON RNG VALUES
//   Both Python (np.random.seed) and C++ (std::mt19937) use Mersenne Twister
//   internally, but their distribution-sampling algorithms differ, so the
//   specific floating-point draws will not be identical — the statistical
//   distributions ARE equivalent.
// ================================================================================

#include <H5Cpp.h>

#include <algorithm>    // std::clamp
#include <array>
#include <chrono>
#include <cmath>        // std::sin
#include <cstdio>       // std::snprintf
#include <filesystem>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

// ─── Global constants — adjust to match your project ─────────────────────────
constexpr int          NUM_RUNS           = 100;
const     std::string  HDF5_DATABASE_FILE = "radar_database.h5";
constexpr float        PAN_INDIA_GRID_MAX = 10'000.0f;

struct GlobalLimits {
    double max_range_px      = 5'000.0;
    double max_antenna_hgt_m =   100.0;
    double max_threat_hgt_m  = 10'000.0;
    double max_elevation_deg =    60.0;
} GLOBAL_LIMITS;
// ─────────────────────────────────────────────────────────────────────────────

// ─── Overloaded helpers: write scalar HDF5 attributes ────────────────────────
static void write_attr(H5::Group& g, const char* name, int value)
{
    H5::DataSpace sp(H5S_SCALAR);
    H5::Attribute a = g.createAttribute(name, H5::PredType::NATIVE_INT, sp);
    a.write(H5::PredType::NATIVE_INT, &value);
}

static void write_attr(H5::Group& g, const char* name, double value)
{
    H5::DataSpace sp(H5S_SCALAR);
    H5::Attribute a = g.createAttribute(name, H5::PredType::NATIVE_DOUBLE, sp);
    a.write(H5::PredType::NATIVE_DOUBLE, &value);
}
// ─────────────────────────────────────────────────────────────────────────────

int main()
{
    // ── STEP 1: GENERATE LOGS DIRECTLY INTO HDF5 CONTAINER ────────────────────
    std::cout << "--- STEP 1: PURGING EXISTING DATABASE AND GENERATING "
              << NUM_RUNS << " HDF5 GROUPS ---\n";

    // ── Clean, safe single-file removal logic ─────────────────────────────────
    if (std::filesystem::exists(HDF5_DATABASE_FILE))
    {
        bool removed = false;
        for (int attempt = 0; attempt < 2 && !removed; ++attempt)
        {
            try {
                std::filesystem::remove(HDF5_DATABASE_FILE);
                std::cout << "Successfully deleted existing database: "
                          << HDF5_DATABASE_FILE << "\n";
                removed = true;
            }
            catch (const std::filesystem::filesystem_error&) {
                std::cout << "Warning: File is locked. "
                             "Attempting retry after brief pause...\n";
                std::this_thread::sleep_for(std::chrono::milliseconds(500));
            }
        }
        if (!removed) {
            std::cerr << "Fatal: Unable to remove locked database file.\n";
            return 1;
        }
    }

    // ── RNG — seeded deterministically to mirror np.random.seed(42) ───────────
    std::mt19937 rng(42);

    // np.random.randint(2, 21) → [2, 21) → [2, 20] inclusive
    std::uniform_int_distribution<int>     radar_count_dist(2, 20);
    std::uniform_real_distribution<double> range_dist     (500.0, GLOBAL_LIMITS.max_range_px);
    std::uniform_real_distribution<double> ant_hgt_dist   (15.0,  GLOBAL_LIMITS.max_antenna_hgt_m);
    std::uniform_real_distribution<double> threat_hgt_dist(100.0, GLOBAL_LIMITS.max_threat_hgt_m);
    std::uniform_real_distribution<double> elev_dist      (15.0,  GLOBAL_LIMITS.max_elevation_deg);
    std::uniform_int_distribution<int>     az_idx_dist    (0, 3);      // index into HORIZ_AZ_CHOICES
    std::uniform_real_distribution<float>  coord_dist     (2000.0f, PAN_INDIA_GRID_MAX - 2000.0f);
    std::uniform_real_distribution<double> noise_dist     (-0.20, 0.20);

    // np.random.choice([90.0, 120.0, 180.0, 360.0])
    constexpr std::array<double, 4> HORIZ_AZ_CHOICES = {90.0, 120.0, 180.0, 360.0};
    constexpr double PI = 3.14159265358979323846;

    // ── Open unified HDF5 container (H5F_ACC_TRUNC = "w") ─────────────────────
    try
    {
        H5::H5File h5_file(HDF5_DATABASE_FILE, H5F_ACC_TRUNC);

        for (int run_id = 1; run_id <= NUM_RUNS; ++run_id)
        {
            const int    fixed_num_radars = radar_count_dist(rng);
            const double radar_range      = range_dist(rng);
            const double antenna_hgt      = ant_hgt_dist(rng);
            const double threat_hgt       = threat_hgt_dist(rng);
            const double horiz_az         = HORIZ_AZ_CHOICES[az_idx_dist(rng)];
            const double vert_el          = elev_dist(rng);

            // Construct specific Group namespace inside HDF5 → "run_0001", etc.
            char group_name[16];
            std::snprintf(group_name, sizeof(group_name), "run_%04d", run_id);
            H5::Group run_group = h5_file.createGroup(group_name);

            // Save structural run parameters as high-performance attributes
            write_attr(run_group, "fixed_radar_deployment_count", fixed_num_radars);
            write_attr(run_group, "antenna_height_m",             antenna_hgt);
            write_attr(run_group, "max_range_px",                 radar_range);
            write_attr(run_group, "horizontal_beam_azimuth_deg",  horiz_az);
            write_attr(run_group, "vertical_beam_elevation_deg",  vert_el);
            write_attr(run_group, "target_threat_height_m",       threat_hgt);

            // Generate asset spatial matrix — shape (fixed_num_radars × 3) float32
            // Columns: [ site_y,  site_x,  site_fitness ]
            // Mirrors: np.zeros((fixed_num_radars, 3), dtype=np.float32)
            const std::size_t n = static_cast<std::size_t>(fixed_num_radars);
            std::vector<float> coordinates_and_fitness(n * 3, 0.0f);

            // elevation_factor does not change per-radar; compute once outside loop
            const double elevation_factor = std::sin(vert_el * PI / 180.0);

            for (std::size_t idx = 0; idx < n; ++idx)
            {
                const float  site_y = coord_dist(rng);
                const float  site_x = coord_dist(rng);

                const double base_fitness =
                    (radar_range / GLOBAL_LIMITS.max_range_px) * 0.50 +
                    elevation_factor                            * 0.20;

                const double noise       = noise_dist(rng);
                const float  site_fitness =
                    static_cast<float>(std::clamp(base_fitness + noise, 0.0, 1.0));

                coordinates_and_fitness[idx * 3 + 0] = site_y;
                coordinates_and_fitness[idx * 3 + 1] = site_x;
                coordinates_and_fitness[idx * 3 + 2] = site_fitness;
            }

            // Save array metrics as an internal high-performance Dataset (float32)
            const hsize_t dims[2] = { static_cast<hsize_t>(fixed_num_radars), 3u };
            H5::DataSpace dataspace(2, dims);
            H5::DataSet   dataset = run_group.createDataSet(
                "assets", H5::PredType::NATIVE_FLOAT, dataspace);
            dataset.write(coordinates_and_fitness.data(), H5::PredType::NATIVE_FLOAT);
        }

        std::cout << "HDF5 Binary Serialization Complete.\n";
    }
    catch (const H5::Exception& e) {
        std::cerr << "HDF5 Error: " << e.getCDetailMsg() << "\n";
        return 1;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
