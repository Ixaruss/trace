#include <cstddef>
#include <vector>
#include <cmath>
#include <chrono>
#include "../common.h"
using namespace std;

/******************************************************************************
 * BRESENHAM LINE
 *
 * Returns grid cells along the line from (x0,y0) to (x1,y1).
 *
 * Classic integer line algorithm.
 * No floating point, no gaps, no duplicates.
 *
 ******************************************************************************/

std::vector<std::pair<int,int>>
Global::bresenham(int x0, int y0, int x1, int y1)
{
    std::vector<std::pair<int,int>> cells;

    int dx = std::abs(x1 - x0);
    int dy = std::abs(y1 - y0);
    int sx = (x0 < x1) ? 1 : -1;
    int sy = (y0 < y1) ? 1 : -1;
    int err = dx - dy;

    while (true)
    {
        cells.push_back({x0, y0});

        if (x0 == x1 && y0 == y1)
            break;

        int e2 = 2 * err;

        if (e2 > -dy) { err -= dy; x0 += sx; }
        if (e2 <  dx) { err += dx; y0 += sy; }
    }

    return cells;
}

