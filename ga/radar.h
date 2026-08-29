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
 *
 * Sector coverage (threat_sector_flag):
 *   RadarEngine has been removed entirely — it required its own
 *   grid-sized statusBuf/fastBuf allocations, which is exactly what
 *   this file exists to avoid. Sector-covered-cell counting is now
 *   done directly on the LocalMask already built for overlap scoring
 *   (see rc_detail::countSectorCovered): the same (2R+3)² buffer is
 *   reused, azimuth-filtered against [loa, roa], no extra ray casts,
 *   no extra allocation.
 */
#pragma once

#include "config.h"
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <cstdint>
#include <cstring>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "../common.h"   // Global, FeasibleGrid, GridPoint, RESOLUTION,
#include "../operations.h"                      // Bin::writeFeasibleGridOverlapPng, OutputCrs, BandLayout, writeOverlapGeoTIFF
// ═══════════════════════════════════════════════════════════════════════════
// Implementation detail — keep out of global namespace
// ═══════════════════════════════════════════════════════════════════════════
namespace rc_detail {

// ── Physical constants ────────────────────────────────────────────────────

// ── Minimal 2-D point ────────────────────────────────────────────────────
struct Pt { int x, y; };

// ── Local bounding-box mask ───────────────────────────────────────────────
// buf is indexed as  buf[(r - rowMin)*localCols + (c - colMin)]
// 1 = covered by this radar, 0 = not
struct LocalMask {
    int rowMin, rowMax;   // [rowMin, rowMax)
    int colMin, colMax;   // [colMin, colMax)
    int localCols;
    int R;                // clamped radius (cells) used to build this mask
    std::vector<uint8_t> buf;
};

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
        push8(x, y);

        // FIX: Calculate the decision variable using current x,y
        // values BEFORE x is incremented.
        if (d > 0) {
            d += 4 * (x - y) + 10;
            --y;
        } else {
            d += 4 * x + 6;
        }
        ++x; // Increment x safely at the very end of the step iteration
    }

    // Keep your exact sorting & deduplication pipeline intact
    std::sort(pts.begin(), pts.end(), [](const Pt& a, const Pt& b){
        return a.x < b.x || (a.x == b.x && a.y < b.y); });

    pts.erase(std::unique(pts.begin(), pts.end(), [](const Pt& a, const Pt& b){
        return a.x == b.x && a.y == b.y; }), pts.end());

    return pts;
}

// ── Azimuth helpers (ported from RadarEngine, operations.h/coverage.cc) ──
//
// North-based clockwise azimuth from source to (row,col).
// North = 0°, East = 90°, South = 180°, West = 270°
inline double normAz(double deg)
{
    double a = std::fmod(deg, 360.0);
    if (a < 0.0) a += 360.0;
    return a;
}

inline double cellAzimuth(int srcRow, int srcCol, int r, int c)
{
    const double drow = r - srcRow;
    const double dcol = c - srcCol;
    // atan2(east, north) = clockwise from north
    double az = std::atan2(dcol, -drow) * RAD2DEG;
    if (az < 0.0) az += 360.0;
    return az;
}

// Returns true if azimuth `az` lies on the arc swept clockwise from
// `leftDeg` to `rightDeg` (both inclusive). North = 0°, clockwise positive.
// Handles the 0°/360° wraparound: if leftDeg > rightDeg the arc is taken
// to cross North.
inline bool inArc(double az, double leftDeg, double rightDeg)
{
    az       = normAz(az);
    leftDeg  = normAz(leftDeg);
    rightDeg = normAz(rightDeg);

    if (leftDeg <= rightDeg)
        return az >= leftDeg && az <= rightDeg;
    // arc wraps through 0°/360°
    return az >= leftDeg || az <= rightDeg;
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
    uint8_t* __restrict__ statBase,
    int localCols,
    int gridCols, bool use_curv)
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
    uint8_t* statCell = statBase + static_cast<std::ptrdiff_t>(srcRow) * localCols + srcCol;

    int    row = srcRow, col = srcCol;
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

        ++steps;
        if (steps > maxRadiusCells) break;

        // True straight-line ground distance from source
        const double drow = static_cast<double>(row - srcRow);
        const double dcol = static_cast<double>(col - srcCol);
        const double dist = std::sqrt(drow * drow + dcol * dcol) * RESOLUTION;
        const double curvDrop = use_curv ? (dist * dist) * INV_2_KEFFR : 0;

        if(use_curv) if (dist > horizonM) break;

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
    double  maxRangeM, bool use_curv)
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
    m.R         = R;
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
                m.localCols, gridCols, use_curv);
    }

    return m;  // moves the buf vector — no copy
}

// ── Count covered cells within a sector, reusing an existing LocalMask ───
//
// leftDeg/rightDeg: arc swept clockwise from leftDeg to rightDeg,
//                    North = 0°. leftDeg < 0 || rightDeg < 0 → full 360°.
//                    leftDeg > rightDeg wraps through North.
//
// Walks only the mask's local bounding box (already (2R+3)² at most —
// no grid-sized buffer, no extra ray casting). This is the direct
// replacement for RadarEngine::compute(...).countCovered.
inline int countSectorCovered(
    const LocalMask& m,
    int srcRow, int srcCol,
    double leftDeg, double rightDeg)
{
    const bool fullCircle = (leftDeg < 0.0 || rightDeg < 0.0);
    const long long rSq = static_cast<long long>(m.R) * m.R;

    int count = 0;
    for (int r = m.rowMin; r < m.rowMax; ++r) {
        const long long dr = r - srcRow;
        const uint8_t* rowBuf = m.buf.data() +
            static_cast<size_t>(r - m.rowMin) * m.localCols;

        for (int c = m.colMin; c < m.colMax; ++c) {
            const long long dc = c - srcCol;
            if (dr * dr + dc * dc > rSq) continue;
            if (!rowBuf[c - m.colMin]) continue;

            if (!fullCircle &&
                !inArc(cellAzimuth(srcRow, srcCol, r, c), leftDeg, rightDeg))
                continue;

            ++count;
        }
    }
    return count;
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
 * If threat_sector_flag is true, also counts (per placement, summed into
 * result[4]) how many of that placement's covered cells fall within the
 * [loa, roa] azimuth sector — computed directly off the same local mask,
 * with no additional allocation or ray casting.
 *
 * Returns result[5]:
 *   [0]  number of placement pairs closer than minSpacingCells
 *   [1]  feasible cells covered by exactly 1 radar
 *   [2]  feasible cells covered by exactly 2 radars
 *   [3]  feasible cells covered by 3+ radars
 *   [4]  sector-covered cell count, summed across placements
 *        (only populated when threat_sector_flag is true; 0 otherwise)
 *
 * Radar spec (fixed): 50 m AGL antenna, 1 000 m target, 10–30° elevation, 300 km range.
 */
inline std::vector<int> calculate_radar_coverage_overlap(
    Global&                       g,
    FeasibleGrid                  grid,               // BY VALUE — feasible cells only; drives result[] stats (GA fitness signal)
    const std::vector<GridPoint>& placements,
    double                        minSpacingCells = 5.0,
    const std::string&            tifPath         = "",
    OutputCrs opt = OutputCrs::Wgs84, BandLayout band = BandLayout::Rgba,
    double                        antennaHeight   = 10.0,    // metres AGL
    double                        minAngle        = -5,     // degrees
    double                        target          = 1000.0,  // metres
    double                        maxAngle        = 90.0,    // degrees
    double                        rangeM          = 25000.0, // metres
    const std::vector<RadarConf>& radarRanges = {},
    FeasibleGrid aoiGrid = {} ,
    bool write_png = false,
    int reduction = 50,
    bool use_curv = true,
    bool threat_sector_flag = false,
    double loa = 0.5,
    double roa = 60
)
{
    std::vector<int> result(5, 0);
    int sectorcount = 0;
    if (placements.empty() || grid.cells.empty()) return result;
    for (size_t i = 0; i < placements.size(); ++i)
        for (size_t j = i + 1; j < placements.size(); ++j) {
            const double dr = placements[i].r - placements[j].r;
            const double dc = placements[i].c - placements[j].c;
            if (std::sqrt(dr*dr + dc*dc) < minSpacingCells) ++result[0];
        }

    // Small helper: apply one radar's mask to a grid's cells (bump overlap on hit).
    auto applyMask = [](const rc_detail::LocalMask& m, std::vector<FeasibleCell>& cells) {
        for (auto& cell : cells) {
            if (cell.r < m.rowMin || cell.r >= m.rowMax ||
                cell.c < m.colMin || cell.c >= m.colMax) continue;
            const size_t local_r = cell.r - m.rowMin;
            const size_t local_c = cell.c - m.colMin;
            const size_t flat = local_r * m.localCols + local_c;
            if (flat < m.buf.size() && m.buf[flat]) ++cell.overlap;
        }
    };

    if (!radarRanges.empty()) {
        // Build the flattened radar list ONCE.
        std::vector<RadarConf> radaran;
        for (const auto& x : radarRanges)
            for (int k = 0; k < x.quantity; ++k)
                radaran.push_back(x);
        if (radaran.empty()) {
            // Nothing configured — bail out instead of UB out-of-bounds reads.
            return result;
        }

        for (size_t i = 0; i < placements.size(); ++i) {
            const auto& p = placements[i];
            const auto& r = radaran[i % radaran.size()]; // cycle if fewer radars than placements

            const rc_detail::LocalMask cov =
                rc_detail::computeLocalMask(
                    g, p.r, p.c,
                    r.antennaHeight,
                    r.target,
                    r.minAngle,
                    r.maxAngle,
                    r.range * 1000,use_curv);

            applyMask(cov, grid.cells);
            applyMask(cov, aoiGrid.cells);

            // Sector-covered-cell count — reuses `cov`, no RadarEngine,
            // no grid-sized buffers, no extra ray casting.
            if (threat_sector_flag)
                sectorcount += rc_detail::countSectorCovered(cov, p.r, p.c, loa, roa);
        }
    } else {
        // Homogenous path — all radars share the same physical spec.
        // rangeM is in metres (default 300 km); previously this was
        // `int range = 25` which passed 25 *metres* and produced no coverage.
        for (const auto& p : placements) {
            const rc_detail::LocalMask mask =
                rc_detail::computeLocalMask(
                    g, p.r, p.c,
                    antennaHeight,
                    target,
                    minAngle,
                    maxAngle,
                    rangeM,use_curv);   // ← was `range` (int 25 = 25 m), now correct metres

            applyMask(mask, grid.cells);
            applyMask(mask, aoiGrid.cells); // same mask, no extra ray-casting cost
            // mask goes out of scope here — bounding-box memory released immediately

            // Sector-covered-cell count — reuses `mask`, no RadarEngine,
            // no grid-sized buffers, no extra ray casting.

            if (threat_sector_flag)
                sectorcount += rc_detail::countSectorCovered(mask, p.r, p.c, loa, roa);
        }
    }


    // ── Bucket final per-cell overlap counts ─────────────────────────────
    // Prefer the full AOI grid when given, so result[]/fitness reflects coverage
    // over the whole Area of Operation, not just the road/rail-feasible subset.
    // (Placement itself is still gated to feasible cells — only the scoring changes.)

    const std::vector<FeasibleCell>& statsCells =
        !aoiGrid.cells.empty() ? aoiGrid.cells : grid.cells;
    for (const auto& cell : statsCells) {
        if      (cell.overlap == 1) ++result[1];
        else if (cell.overlap == 2) ++result[2];
        else if (cell.overlap >= 3) ++result[3];
        // std::cout<< cell.overlap << std::endl;
    }
    result[4] = sectorcount;
    // ── Optional PNG render ───────────────────────────────────────────────
    // aoiGrid (all AOI cells, with overlap now computed on it) drives the picture;
    // grid (feasible-only) is passed through purely to show roads/rails/water/etc.
    if (!tifPath.empty()){
        auto idx = tifPath.find_last_of("\\/");
        std::string p;
        if(std::string::npos != idx){
            p = tifPath.substr(0,idx+1);
        }
        writeOverlapGeoTIFF(g, aoiGrid, tifPath, placements, reduction, opt, band);
        if(write_png){
            Bin::writeFeasibleGridOverlapPng(aoiGrid, p+"result.png", placements, grid,
                    g.bin.flag_road, g.bin.flag_rail,
                    g.bin.flag_water_line, g.bin.flag_water_area,
                    g.bin.flag_ib, g.bin.flag_slope);
        }
    }

    return result;
}

/**
 * get_individual_coverage
 *
 * actual = feasible cells the radar can see (LOS mask buf == 1)
 * ideal  = feasible cells within radius R of the placement
 *          (what a perfectly unobstructed radar would cover)
 *
 * Returns  actual / ideal  in [0.0, 1.0].
 * Returns  0.0             if no feasible cells fall within range.
 */
inline double get_individual_coverage(
    Global&              g,
    const FeasibleGrid&  grid,
    int                  srcRow,
    int                  srcCol,
    double               antennaHeight = 50.0,
    double               targetH       = 1000.0,
    double               minAngle      = 0.5,
    double               maxAngle      = 30.0,
    double               maxRangeM     = 25000.0, bool use_curv = true)
{
    using namespace rc_detail;

    const int gridRows = g.conf.grid.rows;
    const int gridCols = g.conf.grid.cols;

    // ── Build local LOS mask ──────────────────────────────────────────────
    const LocalMask mask = computeLocalMask(
        g, srcRow, srcCol,
        antennaHeight, targetH,
        minAngle, maxAngle,
        maxRangeM, use_curv);

    // Reproduce the same R clamping computeLocalMask applies so the
    // radius we test feasible cells against is identical to the one
    // used when casting rays.
    int R = static_cast<int>(std::ceil(maxRangeM / RESOLUTION));
    R = std::min(R, std::min(
            std::min(srcRow,         gridRows - srcRow - 1),
            std::min(srcCol,         gridCols - srcCol - 1)));

    const double R2 = static_cast<double>(R) * R;

    // ── Count over feasible cells only ────────────────────────────────────
    int actual = 0;
    int ideal  = 0;

    for (const auto& cell : grid.cells) {
        // Must fall within the radar's radius
        const double dr = cell.r - srcRow;
        const double dc = cell.c - srcCol;
        if (dr * dr + dc * dc > R2) continue;

        ++ideal;

        // Check LOS mask — also guards against cells outside bounding box
        if (cell.r < mask.rowMin || cell.r >= mask.rowMax ||
            cell.c < mask.colMin || cell.c >= mask.colMax) continue;

        const size_t flat =
            static_cast<size_t>(cell.r - mask.rowMin) * mask.localCols +
            static_cast<size_t>(cell.c - mask.colMin);

        if (flat < mask.buf.size() && mask.buf[flat] == 1) ++actual;
    }

    if (ideal == 0) return 0.0;
    return static_cast<double>(actual) / static_cast<double>(ideal);
}
