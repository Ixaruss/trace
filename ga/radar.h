/**
 * radar_overlap.h
 *
 * Merged, zero-full-grid-allocation radar coverage overlap.
 * Replaces: operations.h + coverage.cc + compute_local_coverage()
 *
 * Memory per placement:
 *   OLD: 2 × gridRows×gridCols bytes  (statusBuf + fastBuf)
 *   NEW: (2R+3)² bytes                (local bounding-box only)
 *
 * For a 10 000×10 000 DEM at 300 km range (R≈1 200 cells):
 *   OLD ≈ 200 MB held in RadarEngine lifetime
 *   NEW ≈   5.8 MB, allocated / freed per placement
 */
#pragma once

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
                      // Bin::writeFeasibleGridOverlapPng

// ═══════════════════════════════════════════════════════════════════════════
// Implementation detail — keep out of global namespace
// ═══════════════════════════════════════════════════════════════════════════
namespace rc_detail {

// ── Physical constants ────────────────────────────────────────────────────
constexpr double EARTH_R     = 6'371'000.0;
constexpr double K_EFF_R     = (4.0 / 3.0) * EARTH_R;   // 4/3 earth radius
constexpr double INV_2_KEFFR = 1.0 / (2.0 * K_EFF_R);
constexpr double DEG2RAD     = 3.14159265358979323846 / 180.0;

// ── Minimal 2-D point ────────────────────────────────────────────────────
struct Pt { int x, y; };

// ── Local bounding-box mask ───────────────────────────────────────────────
// buf is indexed as  buf[(r - rowMin)*localCols + (c - colMin)]
// 1 = covered by this radar, 0 = not
struct LocalMask {
    int rowMin, rowMax;   // [rowMin, rowMax)
    int colMin, colMax;   // [colMin, colMax)
    int localCols;
    std::vector<uint8_t> buf;
};

// ── Bresenham circle perimeter ────────────────────────────────────────────
inline std::vector<Pt> circleBoundary(int cr, int cc, int R)
{
    std::vector<Pt> pts;
    pts.reserve(8 * R);

    auto push8 = [&](int px, int py) {
        pts.push_back({cr+px, cc+py}); pts.push_back({cr-px, cc+py});
        pts.push_back({cr+px, cc-py}); pts.push_back({cr-px, cc-py});
        pts.push_back({cr+py, cc+px}); pts.push_back({cr-py, cc+px});
        pts.push_back({cr+py, cc-px}); pts.push_back({cr-py, cc-px});
    };

    int x = 0, y = R, d = 3 - 2 * R;
    while (y >= x) {
        push8(x, y); ++x;
        if (d > 0) { --y; d += 4*(x - y) + 10; } else { d += 4*x + 6; }
    }

    std::sort(pts.begin(), pts.end(), [](const Pt& a, const Pt& b){
        return a.x < b.x || (a.x == b.x && a.y < b.y); });
    pts.erase(std::unique(pts.begin(), pts.end(), [](const Pt& a, const Pt& b){
        return a.x == b.x && a.y == b.y; }), pts.end());
    return pts;
}

// ── Binary ray cast into local bounding-box buffer ───────────────────────
//
// statBase is a "virtual origin" pointer such that
//   statBase[row * localCols + col]
// resolves to
//   buf[(row - rowMin) * localCols + (col - colMin)]
//
// The DEM keeps its own global stride (gridCols); only the stat pointer
// uses localCols — this is the key difference from the old castRadarRayFast
// that shared COLS for both.
//
inline void castRayFast(
    int srcRow, int srcCol,
    int endRow, int endCol,
    double srcElev,    double targetHeight,
    double minSlope,   double maxSlope,
    double horizonM,   int maxRadiusCells,
    const short* __restrict__ dem,
    uint8_t*     __restrict__ statBase,
    int localCols,
    int gridCols)
{
    const int adr = std::abs(endRow - srcRow);
    const int adc = std::abs(endCol - srcCol);
    const int sr  = (srcRow < endRow) ?  1 : -1;
    const int sc  = (srcCol < endCol) ?  1 : -1;
    int err = adr - adc;

    // Row strides differ between the global DEM and the local stat buffer
    const auto demStrideR  = static_cast<std::ptrdiff_t>(sr) * gridCols;
    const auto statStrideR = static_cast<std::ptrdiff_t>(sr) * localCols;

    const short* demCell  = dem      + static_cast<std::ptrdiff_t>(srcRow) * gridCols  + srcCol;
    uint8_t*     statCell = statBase + static_cast<std::ptrdiff_t>(srcRow) * localCols + srcCol;

    int    row = srcRow, col = srcCol;
    double dist = 0.0, curvDrop = 0.0;
    double maxTerrainSlope = -1e18;
    int    steps = 0;

    while (true)
    {
        if (row == endRow && col == endCol) break;

        const int  e2 = 2 * err;
        const bool mr = (e2 > -adc);
        const bool mc = (e2 <  adr);

        if (mr) { err -= adc; row += sr; demCell += demStrideR;  statCell += statStrideR; }
        if (mc) { err += adr; col += sc; demCell += sc;           statCell += sc;           }

        const double step = (mr && mc) ? RESOLUTION * 1.41421356237 : RESOLUTION;
        curvDrop += (2.0 * dist * step + step * step) * INV_2_KEFFR;
        dist     += step;
        ++steps;

        if (steps > maxRadiusCells) break;
        if (dist  > horizonM)       break;

        const short raw = *demCell;
        if (raw < -9000) break;   // no-data sentinel

        const double terrainElev  = static_cast<double>(raw) - curvDrop;
        const double terrainSlope = (terrainElev - srcElev) / dist;
        const double targetSlope  = terrainSlope + targetHeight / dist;

        // Covered: not shadowed, within elevation band
        if (maxTerrainSlope <= targetSlope &&
            targetSlope     >= minSlope    &&
            targetSlope     <= maxSlope)
            *statCell = 1;

        if (terrainSlope > maxTerrainSlope) maxTerrainSlope = terrainSlope;
    }
}

// ── Build a local mask for one radar placement ────────────────────────────
// Peak allocation: (2R+3)² bytes — nothing grid-sized.
inline LocalMask computeLocalMask(
    Global& g,
    int     srcRow,     int srcCol,
    double  antennaH,   double targetH,
    double  minElevDeg, double maxElevDeg,
    double  maxRangeM)
{
    const int gridCols = g.conf.grid.cols;
    const int gridRows = g.conf.grid.rows;

    // Clamp radius to grid bounds around source
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

    // ── Allocate bounding box only ────────────────────────────────────────
    LocalMask m;
    m.rowMin    = std::max(0,        srcRow - R - 1);
    m.rowMax    = std::min(gridRows, srcRow + R + 2);
    m.colMin    = std::max(0,        srcCol - R - 1);
    m.colMax    = std::min(gridCols, srcCol + R + 2);
    m.localCols = m.colMax - m.colMin;
    const int localRows = m.rowMax - m.rowMin;

    m.buf.assign(static_cast<size_t>(localRows) * m.localCols, 0);

    // Source cell is always covered
    m.buf[static_cast<size_t>(srcRow - m.rowMin) * m.localCols +
          (srcCol - m.colMin)] = 1;

    // "Virtual origin" base so castRayFast can use global row/col directly
    uint8_t* statBase =
        m.buf.data()
        - static_cast<std::ptrdiff_t>(m.rowMin) * m.localCols
        - m.colMin;

    const short*      demPtr   = g.bin.elev_matrix.data();
    const auto        boundary = circleBoundary(srcRow, srcCol, R);
    const int         N        = static_cast<int>(boundary.size());

    // Rays are independent — safe to parallelise (binary write, no race hazard)
    // Guard against nested parallelism if this function is already called from
    // a parallel region (e.g. per-chromosome GA threads).
#if defined(_OPENMP)
    #pragma omp parallel for schedule(dynamic, 64) if(!omp_in_parallel())
#endif
    for (int i = 0; i < N; ++i) {
        const int er = boundary[i].x, ec = boundary[i].y;
        if (er >= 0 && ec >= 0 && er < gridRows && ec < gridCols)
            castRayFast(
                srcRow, srcCol, er, ec,
                srcElev, targetH,
                minSlope, maxSlope,
                horizonM, R,
                demPtr, statBase,
                m.localCols, gridCols);
    }

    return m;  // moves the buf vector — no copy
}

}  // namespace rc_detail


// ═══════════════════════════════════════════════════════════════════════════
// PUBLIC API
// ═══════════════════════════════════════════════════════════════════════════

/**
 * calculate_radar_coverage_overlap
 *
 * For each placement, casts all radar rays into a local bounding-box mask
 * (≈ (2R+3)² bytes), then increments per-cell overlap counters.
 * No grid-sized buffer is ever allocated.
 *
 * Returns result[4]:
 *   [0]  number of placement pairs closer than minSpacingCells
 *   [1]  feasible cells covered by exactly 1 radar
 *   [2]  feasible cells covered by exactly 2 radars
 *   [3]  feasible cells covered by 3+ radars
 *
 * Radar spec (fixed): 50 m AGL antenna, 1 000 m target, 10–30° elevation, 300 km range.
 */
inline std::vector<int> calculate_radar_coverage_overlap(
    Global&                       g,
    FeasibleGrid                  grid,               // BY VALUE — caller's private copy
    const std::vector<GridPoint>& placements,
    double                        minSpacingCells = 5.0,
    const std::string&            pngPath         = "")
{
    std::vector<int> result(4, 0);
    if (placements.empty() || grid.cells.empty()) return result;

    // ── [0] Under-separation penalty ─────────────────────────────────────
    for (size_t i = 0; i < placements.size(); ++i)
        for (size_t j = i + 1; j < placements.size(); ++j) {
            const double dr = placements[i].r - placements[j].r;
            const double dc = placements[i].c - placements[j].c;
            if (std::sqrt(dr*dr + dc*dc) < minSpacingCells) ++result[0];
        }

    // ── Per-placement: build local mask, tally overlap ───────────────────
    for (const auto& p : placements) {
        const rc_detail::LocalMask mask =
            rc_detail::computeLocalMask(
                g, p.r, p.c,
                /*antennaH*/  50.0,
                /*targetH*/  1'000.0,
                /*minElev*/   0.5,
                /*maxElev*/   30.0,
                /*rangeM*/  300'000.0);

        for (auto& cell : grid.cells) {
            // Quick bounds check before indexing
            if (cell.r < mask.rowMin || cell.r >= mask.rowMax ||
                cell.c < mask.colMin || cell.c >= mask.colMax) continue;

            const size_t flat =
                static_cast<size_t>(cell.r - mask.rowMin) * mask.localCols +
                (cell.c - mask.colMin);

            if (mask.buf[flat]) ++cell.overlap;
        }
        // mask goes out of scope here — bounding-box memory released immediately
    }

    // ── Bucket final per-cell overlap counts ─────────────────────────────
    for (const auto& cell : grid.cells) {
        if      (cell.overlap == 1) ++result[1];
        else if (cell.overlap == 2) ++result[2];
        else if (cell.overlap >= 3) ++result[3];
    }

    // ── Optional PNG render ───────────────────────────────────────────────
    if (!pngPath.empty())
        Bin::writeFeasibleGridOverlapPng(grid, pngPath, placements);

    return result;
}
