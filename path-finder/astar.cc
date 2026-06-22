#include <iostream>
#include <vector>
#include <queue>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

// =====================================================================
// GLOBAL CONFIGURATIONS & TYPE DEFINITIONS
// =====================================================================

// Dimensions of the Pan-India rasterized map grid.
// Each raster cell represents a 100-meter x 100-meter physical region.
const int ROWS = 30000;
const int COLS = 30000;

// Allocation footprint: 30,000 * 30,000 bytes = ~900 MB flat array memory footprint.
const size_t FILE_SIZE = static_cast<size_t>(ROWS) * COLS;

// Discrete point represented as an integer matrix pair: {row_index, column_index}
using Point = std::pair<int, int>;

/**
 * @brief Custom Hash Function Object (Functor) for Point structures.
 * Necessary for storing coordinates within standard unordered containers.
 * Uses bit-shifting to mix row and column values to minimize hash collisions.
 */
struct PointHash {
    std::size_t operator()(const Point& p) const {
        return std::hash<int>{}(p.first) ^ (std::hash<int>{}(p.second) << 1);
    }
};

/**
 * @brief Chebyshev Distance Heuristic Estimator.
 * Calculates the absolute minimum number of cell steps needed to bridge two coordinates.
 * This heuristic is mathematically accurate for an 8-directional routing system where
 * diagonal and straight moves are equally permitted as basic grid steps.
 */
int heuristic(const Point& a, const Point& b) {
    return std::max(std::abs(a.first - b.first), std::abs(a.second - b.second));
}

/**
 * @brief Evaluation Node Structure utilized by the priority queue min-heap.
 */
struct Node {
    Point coord;  // Discrete grid coordinate matrix index
    int f_score;  // Estimated total path cost ($f(n) = g(n) + h(n)$)

    // Overload the greater-than operator to force the std::priority_queue
    // to function as a min-heap (smallest total f_score extracted first).
    bool operator>(const Node& other) const { return f_score > other.f_score; }
};

/**
 * @brief Inline multi-dimensional array row-major flattening index translator.
 * Maps 2D array matrix notation [r][c] cleanly into a 1D pointer offset sequence.
 */
inline size_t get_index(int r, int c) {
    return static_cast<size_t>(r) * COLS + c;
}

// =====================================================================
// SVG ROUTE PLOTTER (DATA VISUALIZATION SUBSYSTEM)
// =====================================================================

/**
 * @brief Generates a responsive vector graphic document (.svg) of the computed route.
 * Dynamically constructs a cropped bounding box viewport encompassing the trip, avoiding
 * rendering the entire 900M cell array to keep file payload sizes light.
 */
void save_route_as_svg(const unsigned char* matrix, const std::vector<Point>& path, const std::string& filename, int padding = 20) {
    if (path.empty()) return;

    Point start = path.front();
    Point end = path.back();

    // Step 1: Establish a tight bounding box around the generated route path coordinates.
    int min_r = ROWS, max_r = 0, min_c = COLS, max_c = 0;
    for (const auto& p : path) {
        min_r = std::min(min_r, p.first);
        max_r = std::max(max_r, p.first);
        min_c = std::min(min_c, p.second);
        max_c = std::max(max_c, p.second);
    }

    // Step 2: Incorporate visibility padding margins safely bounded by absolute map edge limits.
    min_r = std::max(0, min_r - padding);
    max_r = std::min(ROWS - 1, max_r + padding);
    min_c = std::max(0, min_c - padding);
    max_c = std::min(COLS - 1, max_c + padding);

    // Step 3: Implement dynamic downsampling for long-distance trips.
    // Prevents multi-megabyte file sizes by skipping micro-dots on macro views.
    int step = 1;
    if ((max_r - min_r) > 1000 || (max_c - min_c) > 1000) {
        step = 10; // Evaluate and draw every 10th pixel cell context
    }

    // Calculate canvas size metrics based on normalized coordinate bounds
    int width = ((max_c - min_c) / step + 1) * 6;
    int height = ((max_r - min_r) / step + 1) * 6;

    std::ofstream svg(filename);
    if (!svg.is_open()) {
        std::cerr << "Error: Could not create plot file " << filename << "\n";
        return;
    }

    // Initialize SVG canvas header with custom dark developer theme styling
    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width << "\" height=\"" << height
        << "\" viewBox=\"0 0 " << width << " " << height << "\" style=\"background-color:#1e1e24;\">\n";

    // Sub-render A: Draw the static background valid road grid points in light grey
    for (int r = min_r; r <= max_r; r += step) {
        for (int c = min_c; c <= max_c; c += step) {
            if (matrix[get_index(r, c)] == 1) {
                int x = ((c - min_c) / step) * 6 + 3;
                int y = ((r - min_r) / step) * 6 + 3;
                svg << "  <circle cx=\"" << x << "\" cy=\"" << y << "\" r=\"1.5\" fill=\"#5c5c6d\" opacity=\"0.5\"/>\n";
            }
        }
    }

    // Sub-render B: Draw the calculated continuous path vector as a distinct bold green line string
    svg << "  <path d=\"M ";
    for (size_t i = 0; i < path.size(); ++i) {
        double x = ((double)(path[i].second - min_c) / step) * 6.0 + 3.0;
        double y = ((double)(path[i].first - min_r) / step) * 6.0 + 3.0;
        svg << x << " " << y << (i == path.size() - 1 ? "" : " L ");
    }
    svg << "\" fill=\"none\" stroke=\"#38b000\" stroke-width=\"4\" stroke-linecap=\"round\" stroke-linejoin=\"round\"/>\n";

    // Sub-render C: Draw Passenger Pickup target pin geometry (Bright Blue)
    double start_x = ((double)(start.second - min_c) / step) * 6.0 + 3.0;
    double start_y = ((double)(start.first - min_r) / step) * 6.0 + 3.0;
    svg << "  <circle cx=\"" << start_x << "\" cy=\"" << start_y << "\" r=\"7\" fill=\"#00b4d8\"/>\n";
    svg << "  <circle cx=\"" << start_x << "\" cy=\"" << start_y << "\" r=\"2.5\" fill=\"white\"/>\n";

    // Sub-render D: Draw Passenger Dropoff target pin geometry (Vibrant Red)
    double end_x = ((double)(end.second - min_c) / step) * 6.0 + 3.0;
    double end_y = ((double)(end.first - min_r) / step) * 6.0 + 3.0;
    svg << "  <circle cx=\"" << end_x << "\" cy=\"" << end_y << "\" r=\"7\" fill=\"#ff0054\"/>\n";
    svg << "  <circle cx=\"" << end_x << "\" cy=\"" << end_y << "\" r=\"2.5\" fill=\"white\"/>\n";

    svg << "</svg>\n";
    svg.close();
    std::cout << "[Visual] Saved scaled route map plot render to: " << filename << "\n";
}

// =====================================================================
// MAP-MATCHING: SNAP USER TO NEAREST FEASIBLE ROAD CELL
// =====================================================================

/**
 * @brief Evaluates raw geolocation inputs and snaps them to a valid road coordinate.
 * Uses an unweighted Breadth-First Search (BFS) starting from the user's current cell,
 * scanning outward in expanding concentric circles to find the closest road pixel (`1`).
 */
Point snap_to_nearest_road(const unsigned char* matrix, Point target, int max_radius = 200) {
    // Optimization: If the passenger is already standing precisely on a valid road track, return immediately
    if (matrix[get_index(target.first, target.second)] == 1) return target;

    std::queue<Point> q;
    std::unordered_set<Point, PointHash> visited;

    q.push(target);
    visited.insert(target);

    // Compass lookup matrix mapping out surrounding neighbors in an 8-way grid spacing layout
    const std::vector<Point> directions = {
        {-1, 0}, {1, 0}, {0, -1}, {0, 1}, {-1, -1}, {-1, 1}, {1, -1}, {1, 1}
    };

    while (!q.empty()) {
        Point current = q.front();
        q.pop();

        // Enforce maximum scan radius boundaries to prevent boundless searches
        int current_radius = std::max(std::abs(current.first - target.first), std::abs(current.second - target.second));
        if (current_radius > max_radius) break;

        // Process surrounding adjacency structures
        for (const auto& dir : directions) {
            int nr = current.first + dir.first;
            int nc = current.second + dir.second;

            // Verify index coordinates stay safely within spatial matrix bounds
            if (nr >= 0 && nr < ROWS && nc >= 0 && nc < COLS) {
                Point neighbor = {nr, nc};
                if (visited.find(neighbor) == visited.end()) {
                    // Return the coordinate immediately once a valid road pixel is identified
                    if (matrix[get_index(nr, nc)] == 1) return neighbor;

                    visited.insert(neighbor);
                    q.push(neighbor);
                }
            }
        }
    }
    return {-1, -1}; // Returns a sentinel error value if no valid roads are located within range
}

// =====================================================================
// STRICT ROAD-ONLY 8-DIRECTIONAL A* ROUTING ENGINE
// =====================================================================

/**
 * @brief Computes the shortest alternative path linking two nodes across the grid.
 * Operates with strict validation: any pixel evaluation matching a zero value (`0`)
 * is treated as an impassable barrier, completely preventing off-road path creation.
 */
std::vector<Point> astar_mmap(const unsigned char* matrix, Point start, Point end) {
    // Open discovery priority queue sorted by minimum estimated path costs
    std::priority_queue<Node, std::vector<Node>, std::greater<Node>> open_set;

    // Scale tracking costs up by a factor of 10 to manage integers safely without using floats
    open_set.push({start, heuristic(start, end) * 10});

    // Hash maps tracking cost-to-reach scores ($g\_score$) and step linkage lineages ($came\_from$)
    std::unordered_map<Point, int, PointHash> g_score;
    std::unordered_map<Point, Point, PointHash> came_from;
    g_score[start] = 0;

    // Movement mappings for 8-way traversal directions
    const std::vector<Point> directions = {
        {-1, 0}, {1, 0}, {0, -1}, {0, 1}, {-1, -1}, {-1, 1}, {1, -1}, {1, 1}
    };

    // Calculate a tracking corridor box to stop the engine from wandering off-course
    int min_allowed_row = std::min(start.first, end.first) - 500;
    int max_allowed_row = std::max(start.first, end.first) + 500;
    int min_allowed_col = std::min(start.second, end.second) - 500;
    int max_allowed_col = std::max(start.second, end.second) + 500;

    int iterations = 0;

    while (!open_set.empty()) {
        Point current = open_set.top().coord;
        open_set.pop();

        // Print progress indicators every 500,000 steps to monitor long-distance performance
        if (++iterations % 500000 == 0) {
            std::cout << "[Processing] Evaluated " << iterations << " nodes... Current Cell: {"
                      << current.first << ", " << current.second << "}\n";
        }

        // Target Reached Condition: Backtrack and reconstruct the generated routing path
        if (current == end) {
            std::vector<Point> path;
            while (came_from.find(current) != came_from.end()) {
                path.push_back(current);
                current = came_from[current];
            }
            path.push_back(start);
            std::reverse(path.begin(), path.end()); // Revert path into sequential forward order
            return path;
        }

        // Evaluate valid connections for neighboring cells
        for (const auto& dir : directions) {
            int nr = current.first + dir.first;
            int nc = current.second + dir.second;

            if (nr >= 0 && nr < ROWS && nc >= 0 && nc < COLS) {
                // Reject evaluation if the cell drifts outside the tracking corridor bounding box
                if (nr < min_allowed_row || nr > max_allowed_row || nc < min_allowed_col || nc > max_allowed_col) {
                    continue;
                }

                size_t idx = get_index(nr, nc);

                // ==========================================
                // STRICT CONSTRAINT BLOCK
                // ==========================================
                // If the target cell contains a 0 value (representing missing data or off-road spaces),
                // skip it immediately. Bridging gaps or jumping across empty cells is strictly forbidden.
                if (matrix[idx] != 1) {
                    continue;
                }

                // Base cost weights: Straight moves cost 10 units, diagonal moves cost ~14 units ($\sqrt{2} \times 10$).
                int base_move_weight = (dir.first != 0 && dir.second != 0) ? 14 : 10;
                Point neighbor = {nr, nc};
                int tentative_g_score = g_score[current] + base_move_weight;

                // If this path to the neighbor is shorter than any previous attempt, update records
                if (g_score.find(neighbor) == g_score.end() || tentative_g_score < g_score[neighbor]) {
                    came_from[neighbor] = current;
                    g_score[neighbor] = tentative_g_score;

                    // Compute the total estimated cost ($f(n) = g(n) + h(n)$) and insert into the open set
                    int f_score = tentative_g_score + (heuristic(neighbor, end) * 10);
                    open_set.push({neighbor, f_score});
                }
            }
        }
    }
    return {}; // Return an empty vector if no continuous road link exists to reach the destination
}

// =====================================================================
// MAIN EXECUTABLE ENVIRONMENT PIPELINE
// =====================================================================
int main() {
    const char* filepath = "bin/road.bin";

    // Step 1: Open the binary road map data asset
    int fd = open(filepath, O_RDONLY);
    if (fd == -1) {
        std::cerr << "Error: Could not open binary file: " << filepath << "\n";
        return 1;
    }

    // Step 2: Establish Direct Kernel Memory Mapping (`mmap`).
    // Maps the binary payload directly into the process's virtual memory addresses.
    // This bypasses traditional file-read overhead and optimizes page cache performance.
    unsigned char* matrix = (unsigned char*)mmap(NULL, FILE_SIZE, PROT_READ, MAP_SHARED, fd, 0);
    if (matrix == MAP_FAILED) {
        std::cerr << "Error: Memory mapping failed.\n";
        close(fd);
        return 1;
    }

    std::cout << "Successfully mapped " << filepath << " directly into address space.\n";

    // Define coordinate location inputs
    Point passenger_pickup = {2001, 9320};
    Point passenger_dropoff = {3012, 10345};

    std::cout << "\nTaxi Request Received:\n";
    std::cout << "-> Pickup Query : {" << passenger_pickup.first << ", " << passenger_pickup.second << "}\n";
    std::cout << "-> Dropoff Query: {" << passenger_dropoff.first << ", " << passenger_dropoff.second << "}\n";

    // Step 3: Run map-matching to snap raw target selections to valid roads
    std::cout << "\nSnapping points to the nearest road tracks...\n";
    Point start = snap_to_nearest_road(matrix, passenger_pickup, 200);
    Point end = snap_to_nearest_road(matrix, passenger_dropoff, 200);

    if (start.first == -1 || end.first == -1) {
        std::cerr << "Critical Error: Snapping failed. Coordinates are entirely off-grid.\n";
        munmap(matrix, FILE_SIZE);
        close(fd);
        return 1;
    }

    std::cout << "-> Snapped Pickup Point  : {" << start.first << ", " << start.second << "}\n";
    std::cout << "-> Snapped Dropoff Point : {" << end.first << ", " << end.second << "}\n\n";

    // Step 4: Run the strict pathfinder engine across the mapped infrastructure
    std::cout << "Calculating optimal corridor route (Strict Road Tracking Only)...\n";
    std::vector<Point> path = astar_mmap(matrix, start, end);

    // Step 5: Evaluate pathfinding results and present status messages
    if (!path.empty()) {
        std::cout << "\n================================================\n";
        std::cout << "SUCCESS: Route generated successfully!\n";
        std::cout << "Total Path Coordinates : " << path.size() << "\n";
        std::cout << "Estimated Distance     : " << (double)(path.size() - 1) * 100 / 1000.0 << " km.\n";
        std::cout << "================================================\n\n";

        // Save a visual plot of the computed route
        save_route_as_svg(matrix, path, "route_plot.svg", 40);
    } else {
        // Strict failure notice output when road paths contain unresolvable structural gaps
        std::cout << "\n============================================================\n";
        std::cout << "ERROR: Route does not exist.\n";
        std::cout << "Reason: Continuous road path could not be strictly established\n";
        std::cout << "        due to unlinked fragments or gaps in the grid data.\n";
        std::cout << "============================================================\n";
    }

    // Step 6: Clean up system resource handles before exiting cleanly
    munmap(matrix, FILE_SIZE);
    close(fd);
    return 0;
}
