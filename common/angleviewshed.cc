/**
 * @file Viewshed.cc
 *
 * Correct LOS algorithm: max-angle horizon sweep.
 *
 * For each ray from source outward:
 *   - Track running maximum slope seen so far
 *   - Cell is visible if its slope >= maxSlope
 *   - Update maxSlope when a higher slope is found
 *
 * This is physically correct. A cell is visible if the line from
 * the radar to that cell is not obstructed by any closer terrain.
 */

#include "../operations.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

    // ------------------------------------------------------------------
    // CORE RAY — max-angle horizon sweep (LOS visibility)
    //
    // Walk outward from (srcRow, srcCol) toward (endRow, endCol).
    // Track maxSlope = max(dh/dist) seen so far.
    // Cell visible <=> its slope >= maxSlope.
    //
    // dist is accumulated correctly:
    //   axis step     -> RESOLUTION
    //   diagonal step -> RESOLUTION * sqrt(2)
    //
    // Earth curvature: subtract d²/2R from terrain elevation.
    // At 300km this is ~7km — significant for radar.
    //
    // Free function, same shape as RadarEngine::castRadarRayDiag/Fast in
    // coverage.cc, rather than a ViewshedEngine member — the loop only
    // needs the grid geometry and DEM/mask buffers, not any other engine
    // state, so it's kept out of the class like radar's ray casters.
    // ------------------------------------------------------------------
    static void castViewshedRay(
        int srcRow, int srcCol,
        int endRow, int endCol,
        double srcElev,
        int maxRadiusCells,
        const short* __restrict__ dem,
        uint8_t*     __restrict__ mask,
        int COLS, int ROWS, bool use_curvature)
    {
        int adr = std::abs(endRow - srcRow), adc = std::abs(endCol - srcCol);
        int sr  = (srcRow < endRow) ? 1 : -1, sc = (srcCol < endCol) ? 1 : -1;
        int err = adr - adc;

        // layout: elev_matrix[row * COLS + col]
        const std::ptrdiff_t strideR = static_cast<std::ptrdiff_t>(sr) * COLS;
        const std::ptrdiff_t strideC = sc;
        std::ptrdiff_t off = static_cast<std::ptrdiff_t>(srcRow) * COLS + srcCol;

        const short* demCell  = dem  + off;
        uint8_t*     maskCell = mask + off;

        int row = srcRow, col = srcCol;
        double dist     = 0.0;
        double maxSlope = -1e18;

        // incremental curvature: drop = dist² / 2R
        // updated as: drop(d+step) = drop(d) + (2d*step + step²) / 2R
        double curvDrop = 0.0;

        const double maxDist = static_cast<double>(maxRadiusCells) * RESOLUTION;

        while (true)
        {
            if (row == endRow && col == endCol) break;

            int e2 = 2 * err;
            bool mr = (e2 > -adc);
            bool mc = (e2 <  adr);

            if (mr) { err -= adc; row += sr; demCell += strideR; maskCell += strideR; }
            if (mc) { err += adr; col += sc; demCell += strideC; maskCell += strideC; }

            // exact step length
            double step = (mr && mc) ? RESOLUTION * 1.41421356237 : RESOLUTION;

            // incremental curvature update
            if (use_curvature)
                curvDrop += (2.0 * dist * step + step * step) * INV_2EARTHR;

            dist += step;

            if (dist > maxDist) break;

            short rawElev = *demCell;
            if (rawElev < -9000) break;

            // corrected terrain elevation
            double terrainElev = static_cast<double>(rawElev) - curvDrop;

            // slope = tan(elevation angle) — no trig needed, ordering preserved
            double slope = (terrainElev - srcElev) / dist;

            if (slope >= maxSlope)
            {
                *maskCell = 1;
                if (slope > maxSlope)
                    maxSlope = slope;
            }
        }
    }

    // ------------------------------------------------------------------
    // Reuses RadarEngine::circleBoundary (same midpoint-circle algorithm,
    // 8-way symmetry, sort+dedupe) instead of duplicating it here. Both
    // engines share the same row/col grid convention (centerX = row,
    // centerY = col <-> srcRow, srcCol), so the points line up directly.
    // ------------------------------------------------------------------
    void ViewshedEngine::precomputeCirclePerimeter()
    {
        circlePerimeter = RadarEngine::circleBoundary(centerX, centerY, maxRadiusCells);
    }

    // ------------------------------------------------------------------
    // PREFILL
    //
    // Reset viewshedMask to 0 (not visible) — but only across the active
    // bounding box (rowMin..rowMax, colMin..colMax), not the whole grid.
    // Mirrors RadarEngine::prefill() / the per-row memset reset in
    // RadarEngine::computeFast(): the mask here is already binary, so
    // there's no separate "outside sector" value to write — a sector
    // that's smaller than 360° just leaves more of the reset area at 0,
    // which castViewshedRay never touches anyway. The win is the same as
    // radar's: repeated computeViewshed() calls on a small radius stay
    // cheap instead of paying for a full-grid std::fill every time.
    // ------------------------------------------------------------------
    void ViewshedEngine::prefill()
    {
        #pragma omp parallel for schedule(static, 32)
        for (int r = rowMin; r < rowMax; r++) {
            uint8_t* p = viewshedMask.data() + (size_t)r * gridCols + colMin;
            std::memset(p, 0, colMax - colMin);
        }
    }



    // ------------------------------------------------------------------
    // SECTOR SUPPORT
    //
    // A sector is an arbitrary arc given as two azimuths (degrees,
    // North = 0°, clockwise positive):
    //
    //   leftDeg  = start of the arc (left edge, looking outward from source)
    //   rightDeg = end of the arc   (right edge)
    //
    // The scanned region is everything swept clockwise from leftDeg to
    // rightDeg. If leftDeg > rightDeg the arc wraps through North
    // (e.g. leftDeg=340°, rightDeg=20° scans 340->360->0->20, an
    // 40°-wide arc straddling North).
    //
    // leftDeg < 0 || rightDeg < 0 -> full 360°.
    //
    // Reuses RadarEngine::cellAzimuth / inArc / normAz rather than
    // duplicating them, since ViewshedEngine uses the same row/col grid
    // layout (centerX = row, centerY = col) as RadarEngine.
    // ------------------------------------------------------------------
    std::vector<Point2D> ViewshedEngine::sectorPerimeter(double leftDeg, double rightDeg) const
    {
        if (leftDeg < 0.0 || rightDeg < 0.0) return circlePerimeter;  // full 360

        std::vector<Point2D> out;
        double width = rightDeg - leftDeg;
        if (width < 0.0) width += 360.0;
        out.reserve(static_cast<size_t>(circlePerimeter.size() * (width / 360.0)) + 8);

        for (const auto& pt : circlePerimeter) {
            double az = RadarEngine::cellAzimuth(centerX, centerY, pt.x, pt.y);
            if (RadarEngine::inArc(az, leftDeg, rightDeg))
                out.push_back(pt);
        }
        return out;
    }

    // ------------------------------------------------------------------
    // leftDeg/rightDeg: arc swept clockwise from leftDeg to rightDeg,
    //                    North = 0°. leftDeg < 0 || rightDeg < 0 -> full 360.
    //                    leftDeg > rightDeg wraps through North.
    // ------------------------------------------------------------------
    void ViewshedEngine::computeViewshed(double leftDeg, double rightDeg)
    {
        auto t0 = std::chrono::high_resolution_clock::now();

        bool fullCircle = (leftDeg < 0.0 || rightDeg < 0.0);
        auto activePerimeter = sectorPerimeter(leftDeg, rightDeg);

        std::cout << "[INFO]   Sector  : ";
        if (fullCircle) {
            std::cout << "360° (full)\n";
        } else {
            std::cout << RadarEngine::normAz(leftDeg) << "° -> "
                      << RadarEngine::normAz(rightDeg)
                      << "° (clockwise)  rays=" << activePerimeter.size() << "\n";
        }

        // Reset mask so repeated calls (e.g. different sectors on the same
        // engine) don't leak stale visible cells from a previous call.
        // Bounding-box prefill instead of a full-grid std::fill — see
        // ViewshedEngine::prefill().
        prefill();

        // centerX = row, centerY = col
        size_t srcIdx =
            static_cast<size_t>(centerX) * gridCols + centerY;

        short srcTerrain = g.bin.elev_matrix[srcIdx];
        if (srcTerrain < -9000) {
            std::cerr << "[ERROR] No-data elevation at radar cell\n";
            return;
        }

        double srcElev =
            static_cast<double>(srcTerrain) + radarHeightMeters;

        // source cell always visible
        viewshedMask[srcIdx] = 1;

        const size_t N       = activePerimeter.size();
        const short* demPtr  = g.bin.elev_matrix.data();
        uint8_t*     maskPtr = viewshedMask.data();
        const int    COLS    = gridCols;
        const int    ROWS    = gridRows;
        const int    MRC     = maxRadiusCells;
        const bool   UC      = useEarthCurvature;

        // rays are independent — safe to parallelize
        // each ray only writes cells it visits (no two rays share a cell
        // on the same Bresenham path), and both write value=1 if they do
        #pragma omp parallel for schedule(dynamic, 64) default(none) \
            shared(demPtr, maskPtr, activePerimeter) \
            firstprivate(N, COLS, ROWS, srcElev, MRC, UC, centerX, centerY)
        for (size_t i = 0; i < N; ++i)
        {
            // activePerimeter stores (row, col) in (.x, .y) — see
            // precomputeCirclePerimeter / RadarEngine::circleBoundary.
            int er = activePerimeter[i].x;
            int ec = activePerimeter[i].y;
            if (er >= 0 && er < ROWS && ec >= 0 && ec < COLS)
                castViewshedRay(centerX, centerY, er, ec, srcElev, MRC,
                                 demPtr, maskPtr, COLS, ROWS, UC);
        }

        auto t1 = std::chrono::high_resolution_clock::now();
        std::cout << "[METRIC] Viewshed compute             : "
                  << std::chrono::duration<double,std::milli>(t1-t0).count()
                  << " ms\n";

        // stats
        int vis = 0;
        for (auto v : viewshedMask) vis += v;
        std::cout << "[INFO]   Visible cells: " << vis << "\n";
    }





/******************************************************************************
 * MAIN
 ******************************************************************************/

// int main()
// {
//     auto t0 = std::chrono::high_resolution_clock::now();

//     Global g;

//     // Srinagar airport
//     auto engine = std::make_unique<ViewshedEngine>(
//         g,
//         33.9871, 74.7742,  // lon
//         300.0f,       // range km
//         30.0f,        // antenna height m
//         true          // earth curvature
//     );

//     // full 360
//     engine->computeViewshed();
//     engine->writeViewshedToSVG("viewshed_timed.svg");

//     // // 60°-wide arc centred on NE (20° -> 80°)
//     // engine->computeViewshed(20, 80);
//     // engine->writeViewshedToSVG("viewshed_sector_ne.svg");

//     // // 40°-wide arc straddling North (340° -> 20°)
//     // engine->computeViewshed(340, 20);
//     // engine->writeViewshedToSVG("viewshed_sector_north.svg");

//     auto t1 = std::chrono::high_resolution_clock::now();
//     std::cout << "\n=======================================================\n"
//               << ">> TOTAL: "
//               << std::chrono::duration<double,std::milli>(t1-t0).count()/1000.0
//               << " s\n=======================================================\n";
// }
