#include "radar.h"

int main() {
    Global g;
    int x,y;
    g.bin.getGridCoords(33.9871, 74.7742,x,y);
    std::vector<GridPoint> placements = {{x,y}};
    auto f = Bin::getFeasibleGrid(g.bin);
    calculate_radar_coverage_overlap(g, f, placements,5,"main.png");
}
