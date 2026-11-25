// #include "bin.h"
// using namespace std;
// int main() {
//     Bin b;
//     FeasibleGrid grid = b.getFeasibleGrid(WRITE, FULL_MATRIX);
//     // find cells with any bit set
//     int found = 0;
//     for (auto& cell : grid.cells) {
//         if (cell.value > 0 && found < 10) {
//             cout << "cell r=" << cell.r << " c=" << cell.c
//                  << " road="  << ((cell.value >> BIT_ROAD)       & 1)
//                  << " rail="  << ((cell.value >> BIT_RAIL)       & 1)
//                  << " wline=" << ((cell.value >> BIT_WATER_LINE) & 1)
//                  << " warea=" << ((cell.value >> BIT_WATER_AREA) & 1)
//                  << " ib="    << ((cell.value >> BIT_IB)         & 1)
//                  << " raw="   << (int)cell.value
//                  << "\n";
//             found++;
//         }
//     }
//     cout << "total non-zero cells: " <<
//         count_if(grid.cells.begin(), grid.cells.end(),
//                  [](const FeasibleCell& c){ return c.value > 0; }) << "\n";
//     Bin::verifyMatrixBin<FeasibleCell>("feasable.bin", 1, grid.cells.size());
//     return 0;
// }

// #include "bin.h"
// #include <chrono>

// using namespace lib;
// int main() {
//     Bin b;
//     b.visualize("./image/combined.png");
//     return 0;
// }
// #include "bin.h"

// using namespace std;

// int main() {
//     auto k = Bin::getFeasibleGrid(WRITE,FULL_MATRIX);
//     Bin::verifyMatrixBin<FeasibleCell>("./feasable.bin",1,k.cells.size());
//     return 0;
// }

#include "bin.h"

#include "../operations.h"
using namespace std;

int main() {
    Global g;
    g.bin.init();
    cout << "elev at lat: 33.987 lng: 74.7742 elev: " << g.bin.getElevation(33.9871, 74.7742) << "\n";
    // int x,y;
    // g.bin.getGridCoords(33.9871, 74.7742,x,y);
    // Bin::getonlyFeasibleGrid(b); // output: feasible.bin in root folder based on constraints defined for PSO
    // Bin::write_heatmap_bin_h5(); // withall option for all cells without feasibility check for heatmap : output aoi_cells.h5 in root folder
    return 0;
}

// #include "bin.h"

// using namespace std;

// int main() {
//     Bin b;
//     Bin::build_elevation_matrix();
//     Bin::verifyMatrixBin<int16_t>("./bin/elev.bin", b.conf.grid.rows, b.conf.grid.cols);
//     return 0;
// }
