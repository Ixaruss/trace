
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <cstdint>
#include <cstring>

#if defined(_OPENMP)
#  include <omp.h>
#endif

#include "../common.h"   // Global, FeasibleGrid, GridPoint, RESOLUTION,
#include "../ga/radar.h"

using namespace rc_detail;
// ═══════════════════════════════════════════════════════════════════════════
// ADDED: reverse ray cast (target -> source)
// ═══════════════════════════════════════════════════════════════════════════
//
// Whether a cell is visible from the source depends on every closer cell on
// the ray (the running shadow slope), so that can only be computed walking
// source -> target, exactly like castRayFast does. What this function does
// "in reverse" is: run that exact forward visibility test once (identical
// formula/order to castRayFast, so results agree with it cell-for-cell),
// then scan the result back-to-front (target -> source) to find the
// FARTHEST visible cell — that is the "break point". Coverage is then
// written for the whole segment from the source up to and including that
// break point in one pass.
//
// This keeps the outer boundary of the filled region aligned with the real
// LOS coverage boundary (the farthest-visible-cell-per-ray is exactly what
// defines that boundary), instead of just filling to the R-circle regardless
// of terrain.
//
// This does not touch or replace castRayFast / computeLocalMask — it is an
// additive, independent code path.
inline void castRayReverse(
    int srcRow, int srcCol,
    int endRow, int endCol,
    double srcElev,    double targetHeight,
    double minSlope,   double maxSlope,
    double horizonM,   int maxRadiusCells,
    const short* __restrict__ dem,
    uint8_t* __restrict__ statBase,
    int localCols,
    int gridCols)
{
    // Walk the Bresenham path from source -> target once, recording the
    // ordered list of visited cells (same stepping rule as castRayFast).
    struct RayCell { int row, col; };
    std::vector<RayCell> path;
    path.reserve(static_cast<size_t>(maxRadiusCells) + 2);

    {
        const int adr = std::abs(endRow - srcRow);
        const int adc = std::abs(endCol - srcCol);
        const int sr  = (srcRow < endRow) ?  1 : -1;
        const int sc  = (srcCol < endCol) ?  1 : -1;
        int err = adr - adc;
        int row = srcRow, col = srcCol;

        path.push_back({row, col});   // source cell itself

        int steps = 0;
        while (!(row == endRow && col == endCol)) {
            const int  e2 = 2 * err;
            const bool mr = (e2 > -adc);
            const bool mc = (e2 <  adr);

            if (mr) { err -= adc; row += sr; }
            if (mc) { err += adr; col += sc; }

            path.push_back({row, col});

            if (++steps > maxRadiusCells) break;
        }
    }

    // ── Forward pass: identical test/order to castRayFast, so "visible[i]"
    // agrees with what castRayFast would have marked at that cell. Stops
    // permanently at the first horizon-exceeded / no-data cell, exactly as
    // castRayFast's own loop does (nothing past that point is ever reached).
    std::vector<uint8_t> visible(path.size(), 0);
    double maxTerrainSlope = -1e18;
    int    lastValid        = 0;

    for (size_t i = 1; i < path.size(); ++i) {
        const int row = path[i].row;
        const int col = path[i].col;

        const double drow = static_cast<double>(row - srcRow);
        const double dcol = static_cast<double>(col - srcCol);
        const double dist = std::sqrt(drow * drow + dcol * dcol) * RESOLUTION;

        if (dist > horizonM) break;   // matches castRayFast: walk stops here

        const short raw = dem[static_cast<std::ptrdiff_t>(row) * gridCols + col];
        if (raw < -9000) break;       // matches castRayFast: walk stops here

        const double curvDrop     = (dist * dist) * INV_2_KEFFR;
        const double terrainElev  = static_cast<double>(raw) - curvDrop;
        const double terrainSlope = (terrainElev - srcElev) / dist;
        const double targetSlope  = terrainSlope + targetHeight / dist;

        if (maxTerrainSlope <= targetSlope &&
            targetSlope     >= minSlope    &&
            targetSlope     <= maxSlope)
            visible[i] = 1;

        if (terrainSlope > maxTerrainSlope) maxTerrainSlope = terrainSlope;

        lastValid = static_cast<int>(i);
    }

    // ── Backward search: scan target -> source over the visibility we just
    // computed, and stop at the first (i.e. farthest-from-source) visible
    // cell. That is the break point.
    int breakIdx = -1;
    for (int i = lastValid; i >= 1; --i) {
        if (visible[i]) { breakIdx = i; break; }
    }

    if (breakIdx < 0) return;   // nothing visible anywhere along this ray

    // Fill source (index 0) through breakIdx inclusive as covered.
    for (int i = 0; i <= breakIdx; ++i) {
        uint8_t* cell = statBase +
            static_cast<std::ptrdiff_t>(path[i].row) * localCols + path[i].col;
        *cell = 1;
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// ADDED: reverse ray cast that returns the break point instead of writing a mask
// ═══════════════════════════════════════════════════════════════════════════
//
// Same corrected approach as castRayReverse: visibility genuinely depends on
// every closer cell, so it's computed with one forward pass identical to
// castRayFast's own test; the result is then scanned back-to-front to find
// the farthest visible cell, which is reported as the break point.
//
// found is set to true/false depending on whether any visible cell was
// found along this ray; the return value is only meaningful when found==true.
 Pt findBreakPointReverse(
    int srcRow, int srcCol,
    int endRow, int endCol,
    double srcElev,    double targetHeight,
    double minSlope,   double maxSlope,
    double horizonM,   int maxRadiusCells,
    const short* __restrict__ dem,
    int gridCols,
    bool& found)
{
    found = false;

    struct RayCell { int row, col; };
    std::vector<RayCell> path;
    path.reserve(static_cast<size_t>(maxRadiusCells) + 2);

    {
        const int adr = std::abs(endRow - srcRow);
        const int adc = std::abs(endCol - srcCol);
        const int sr  = (srcRow < endRow) ?  1 : -1;
        const int sc  = (srcCol < endCol) ?  1 : -1;
        int err = adr - adc;
        int row = srcRow, col = srcCol;

        path.push_back({row, col});

        int steps = 0;
        while (!(row == endRow && col == endCol)) {
            const int  e2 = 2 * err;
            const bool mr = (e2 > -adc);
            const bool mc = (e2 <  adr);

            if (mr) { err -= adc; row += sr; }
            if (mc) { err += adr; col += sc; }

            path.push_back({row, col});

            if (++steps > maxRadiusCells) break;
        }
    }

    // ── Forward pass: identical test/order to castRayFast ─────────────────
    std::vector<uint8_t> visible(path.size(), 0);
    double maxTerrainSlope = -1e18;
    int    lastValid        = 0;

    for (size_t i = 1; i < path.size(); ++i) {
        const int row = path[i].row;
        const int col = path[i].col;

        const double drow = static_cast<double>(row - srcRow);
        const double dcol = static_cast<double>(col - srcCol);
        const double dist = std::sqrt(drow * drow + dcol * dcol) * RESOLUTION;

        if (dist > horizonM) break;   // matches castRayFast: walk stops here

        const short raw = dem[static_cast<std::ptrdiff_t>(row) * gridCols + col];
        if (raw < -9000) break;       // matches castRayFast: walk stops here

        const double curvDrop     = (dist * dist) * INV_2_KEFFR;
        const double terrainElev  = static_cast<double>(raw) - curvDrop;
        const double terrainSlope = (terrainElev - srcElev) / dist;
        const double targetSlope  = terrainSlope + targetHeight / dist;

        if (maxTerrainSlope <= targetSlope &&
            targetSlope     >= minSlope    &&
            targetSlope     <= maxSlope)
            visible[i] = 1;

        if (terrainSlope > maxTerrainSlope) maxTerrainSlope = terrainSlope;

        lastValid = static_cast<int>(i);
    }

    // ── Backward search: farthest visible cell = the break point ──────────
    for (int i = lastValid; i >= 1; --i) {
        if (visible[i]) {
            found = true;
            return Pt{path[i].row, path[i].col};
        }
    }

    return Pt{0, 0};   // unused when found == false
}





// ═══════════════════════════════════════════════════════════════════════════
// ADDED: build a local mask using the reverse (target -> source) ray cast
// ═══════════════════════════════════════════════════════════════════════════
// Identical setup/allocation to computeLocalMask — only the per-ray function
// differs (castRayReverse instead of castRayFast). Left as a separate
// function so computeLocalMask itself is untouched.
LocalMask computeLocalMaskReverse(
    Global& g,
    int     srcRow,     int srcCol,
    double  antennaH,   double targetH,
    double  minElevDeg, double maxElevDeg,
    double  maxRangeM)
{
    const int gridCols = g.conf.grid.cols;
    const int gridRows = g.conf.grid.rows;

    int R = static_cast<int>(std::ceil(maxRangeM / RESOLUTION));
    R = std::min(R, std::min(
            std::min(srcRow,         gridRows - srcRow - 1),
            std::min(srcCol,         gridCols - srcCol - 1)));

    const double horizonM =
        std::min(maxRangeM,
                 std::sqrt(2.0 * K_EFF_R * antennaH) +
                 std::sqrt(2.0 * K_EFF_R * targetH));

    const double minSlope = std::tan(minElevDeg * DEG2RAD);
    const double maxSlope = std::tan(maxElevDeg * DEG2RAD);

    const double srcElev =
        static_cast<double>(g.bin.elev_matrix[
            static_cast<size_t>(srcRow) * gridCols + srcCol]) + antennaH;

    LocalMask m;
    m.rowMin    = std::max(0,        srcRow - R - 1);
    m.rowMax    = std::min(gridRows, srcRow + R + 2);
    m.colMin    = std::max(0,        srcCol - R - 1);
    m.colMax    = std::min(gridCols, srcCol + R + 2);
    m.localCols = m.colMax - m.colMin;
    const int localRows = m.rowMax - m.rowMin;

    m.buf.assign(static_cast<size_t>(localRows) * m.localCols, 0);

    m.buf[static_cast<size_t>(srcRow - m.rowMin) * m.localCols +
          (srcCol - m.colMin)] = 1;

    uint8_t* statBase =
        m.buf.data()
        - static_cast<std::ptrdiff_t>(m.rowMin) * m.localCols
        - m.colMin;

    const short*      demPtr   = g.bin.elev_matrix.data();
    const auto        boundary = circleBoundary(srcRow, srcCol, R);
    const int         N        = static_cast<int>(boundary.size());

#if defined(_OPENMP)
    #pragma omp parallel for schedule(dynamic, 64) if(!omp_in_parallel())
#endif
    for (int i = 0; i < N; ++i) {
        const int er = boundary[i].x, ec = boundary[i].y;
        if (er >= 0 && ec >= 0 && er < gridRows && ec < gridCols)
            castRayReverse(
                srcRow, srcCol, er, ec,
                srcElev, targetH,
                minSlope, maxSlope,
                horizonM, R,
                demPtr, statBase,
                m.localCols, gridCols);
    }

    return m;  // moves the buf vector — no copy
}










/**
 * ═══════════════════════════════════════════════════════════════════════════
 * ADDED: get_coverage_reverse
 * ═══════════════════════════════════════════════════════════════════════════
 *
 * For every ray cast from the source out to the R-circle boundary, scans
 * that ray from the TARGET/boundary end back toward the SOURCE
 * (rc_detail::findBreakPointReverse) and finds the first cell (walking from
 * far -> near) that is visible. That single cell — the "break point" — is
 * collected into the returned list. Rays along which nothing is visible
 * contribute no entry.
 *
 * Returns: one GridPoint per boundary ray that had a break point (not an
 * actual/ideal ratio). No local mask buffer is built — this is a direct,
 * per-ray scan.
 *
 * This function does not modify or call get_individual_coverage — it is a
 * fully independent, additive code path built on castRayReverse machinery.
 */
 std::vector<GridPoint> RadarEngine::get_coverage_reverse(
    Global&              g,
    int                  srcRow,
    int                  srcCol,
    double               antennaHeight,
    double               targetH,
    double               minAngle,
    double               maxAngle,
    double               maxRangeM,
    const std::string&   pngPath       )   // ADDED — optional PNG of source + break points
{
    using namespace rc_detail;

    const int gridRows = g.conf.grid.rows;
    const int gridCols = g.conf.grid.cols;

    int R = static_cast<int>(std::ceil(maxRangeM / RESOLUTION));
    R = std::min(R, std::min(
            std::min(srcRow,         gridRows - srcRow - 1),
            std::min(srcCol,         gridCols - srcCol - 1)));

    const double horizonM =
        std::min(maxRangeM,
                 std::sqrt(2.0 * K_EFF_R * antennaHeight) +
                 std::sqrt(2.0 * K_EFF_R * targetH));

    const double minSlope = std::tan(minAngle * DEG2RAD);
    const double maxSlope = std::tan(maxAngle * DEG2RAD);

    const double srcElev =
        static_cast<double>(g.bin.elev_matrix[
            static_cast<size_t>(srcRow) * gridCols + srcCol]) + antennaHeight;

    const short* demPtr   = g.bin.elev_matrix.data();
    const auto   boundary = circleBoundary(srcRow, srcCol, R);
    const int    N        = static_cast<int>(boundary.size());

    std::vector<GridPoint> breakPoints;
    breakPoints.reserve(N);

    for (int i = 0; i < N; ++i) {
        const int er = boundary[i].x, ec = boundary[i].y;
        if (er < 0 || ec < 0 || er >= gridRows || ec >= gridCols) continue;

        bool found = false;
        const Pt bp = findBreakPointReverse(
            srcRow, srcCol, er, ec,
            srcElev, targetH,
            minSlope, maxSlope,
            horizonM, R,
            demPtr, gridCols,
            found);

        if (found)
            breakPoints.push_back(GridPoint{bp.x, bp.y});
    }

    // ── Optional PNG render ────────────────────────────────────────────────
    // Bounding box mirrors the one computeLocalMaskReverse would allocate
    // (srcRow/srcCol ± R + pad), so the picture frames exactly the area rays
    // were cast across. Requires no FeasibleGrid / terrain-value grid —
    // just the source and the break points already computed above.
    if (!pngPath.empty()) {
        const int rowMin = std::max(0,        srcRow - R - 1);
        const int rowMax = std::min(gridRows, srcRow + R + 2);
        const int colMin = std::max(0,        srcCol - R - 1);
        const int colMax = std::min(gridCols, srcCol + R + 2);

        Bin::writeReverseBreakPointsPng(
            pngPath,
            rowMin, colMin,
            colMax - colMin, rowMax - rowMin,
            srcRow, srcCol,
            breakPoints);
    }

    return breakPoints;
}
