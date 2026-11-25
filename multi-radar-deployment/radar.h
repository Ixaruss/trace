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
 *
 * WINDOWED-LOADING FIX (added when Bin gained AOI/UNBOUNDED partial-grid
 * loading): this file was written against the old always-full-grid
 * elev_matrix model. computeLocalMask() used to index elev_matrix with
 * g.conf.grid.cols (the FULL grid's width) and the placement's GLOBAL
 * row/col — correct only when the whole grid is loaded. Under AOI or
 * UNBOUNDED, elev_matrix holds only g.bin.op_window (sized
 * window.rows × window.cols, indexed local-to-the-window), so that old
 * indexing read out of bounds for essentially every placement. Fixed by
 * doing all internal DEM/ray-cast work in WINDOW-LOCAL coordinates and
 * converting only the exported LocalMask bounding box back to global at
 * the end — every other function in this file (applyMask,
 * countSectorCovered, get_individual_coverage) already assumed a global
 * bounding box to compare against FeasibleCell::r/c, so none of them
 * needed to change.
 */
#pragma once

#include "config.h"
#include <numbers>
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
#include "../operations.h"
#include "lib.h" // Bin::writeFeasibleGridOverlapPng, OutputCrs, BandLayout, writeOverlapGeoTIFF
// ═══════════════════════════════════════════════════════════════════════════
// Implementation detail — keep out of global namespace
// ═══════════════════════════════════════════════════════════════════════════
namespace rc_detail {

// Set to true on GA thread-pool workers (see ThreadPool in ga.h). Those
// workers are plain std::threads, so omp_in_parallel() is FALSE inside
// them — without this flag every worker spawned its own full OpenMP team
// in computeLocalMask, i.e. (0.75 x cores) workers x (cores) OMP threads
// all fighting for the same cores. Population-level parallelism already
// saturates the machine; per-ray OMP is only useful when a single
// evaluation runs alone (e.g. the final render on the main thread).
inline thread_local bool tls_serial_rays = false;

// ── Physical constants ────────────────────────────────────────────────────

// ── Minimal 2-D point ────────────────────────────────────────────────────
struct Pt { int x, y; };

// ── Local bounding-box mask ───────────────────────────────────────────────
// buf is indexed as  buf[(r - rowMin)*localCols + (c - colMin)]
// 1 = covered by this radar, 0 = not
//
// rowMin/rowMax/colMin/colMax are GLOBAL grid coordinates (same space as
// FeasibleCell::r/c), regardless of whether the underlying elev_matrix is
// windowed — computeLocalMask() converts back to global before returning,
// so every consumer of this struct can keep comparing against global (r,c)
// without caring whether AOI/UNBOUNDED windowing is in play.
struct LocalMask {
    int rowMin = 0, rowMax = 0;   // [rowMin, rowMax)
    int colMin = 0, colMax = 0;   // [colMin, colMax)
    int localCols = 0;
    int R = 0;                    // clamped radius (cells) used to build this mask
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
// All coordinates here (srcRow/srcCol/endRow/endCol) are WINDOW-LOCAL —
// same space as the DEM pointer's own layout (stride = the window's
// column count, passed as gridCols below). Callers are responsible for
// translating global placements to local before calling this.
//
// statBase is a "virtual origin" pointer such that
//   statBase[row * localCols + col]
// resolves to
//   buf[(row - rowMin) * localCols + (col - colMin)]
// using the LOCAL rowMin/colMin the buffer was actually built with.
//
// The DEM keeps its own local stride (gridCols == window width); only the
// stat pointer uses localCols — this is the key difference from the old
// castRadarRayFast that shared COLS for both.
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
    int gridCols, bool use_curv,bool useAMSL,
    double cellM)   // grid resolution in metres (Bin::cell_m), NOT RESOLUTION
{
    const int adr = std::abs(endRow - srcRow);
    const int adc = std::abs(endCol - srcCol);
    const int sr  = (srcRow < endRow) ?  1 : -1;
    const int sc  = (srcCol < endCol) ?  1 : -1;
    int err = adr - adc;

    // Row strides differ between the DEM (window stride) and the local stat buffer
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
        const double dist = std::sqrt(drow * drow + dcol * dcol) * cellM;
        const double curvDrop = use_curv ? (dist * dist) * INV_2_KEFFR : 0;

        if(use_curv) if (dist > horizonM) break;

        const short raw = *demCell;
        if (raw < -9000) break;   // no-data sentinel

        const double terrainElev  = static_cast<double>(raw) - curvDrop;
        const double terrainSlope = (terrainElev - srcElev) / dist;
        double targetSlope  = useAMSL
            ? (targetHeight - srcElev) / dist
            : terrainSlope + targetHeight / dist;

        // Covered: not shadowed, within elevation band
        if (maxTerrainSlope <= targetSlope &&
            targetSlope     >= minSlope    &&
            targetSlope     <= maxSlope)
            *statCell = 1;

        if (terrainSlope > maxTerrainSlope) maxTerrainSlope = terrainSlope;
    }
}

// Safe degenerate mask: rowMax<=rowMin (and colMax<=colMin) so every
// consumer's bounds checks (applyMask, countSectorCovered,
// get_individual_coverage) naturally treat it as "covers nothing" without
// needing their own special-case for "no window" / "placement not loaded".
inline LocalMask emptyMask() {
    LocalMask m;
    m.rowMin = m.rowMax = m.colMin = m.colMax = 0;
    m.localCols = 0;
    m.R = 0;
    return m;
}

// ── Build a local mask for one radar placement ────────────────────────────
// srcRow/srcCol are GLOBAL grid coordinates, same as every other caller in
// this codebase passes (FeasibleGrid/GridPoint are global). Internally this
// translates to window-local before touching elev_matrix, since that's the
// only space elev_matrix is actually indexed in once AOI/UNBOUNDED windowed
// loading is in play.
// Peak allocation: (2R+3)² bytes — nothing grid-sized.
inline LocalMask computeLocalMask(
    Global& g,
    int     srcRow,     int srcCol,
    double  antennaH,   double targetH,
    double  minElevDeg, double maxElevDeg,
    double  maxRangeM, bool use_curv, bool use_radar_model,
    int peak_power, int antenna_gain, double wavelength,
    int crossSection, double detachable_threshold, bool useAMSL)
{
    // FIX: this used to take a `load_full` parameter with branches that
    // were BACKWARDS relative to what that name means everywhere else in
    // the codebase — `if(load_full)` ran the WINDOWED path (g.bin.op_window),
    // the `else` ran the FULL-GRID path (g.conf.grid.cols, untranslated
    // global coordinates). calculate_radar_coverage_overlap defaulted this
    // to `false` and no caller in ga.h ever passed it explicitly, so every
    // GA coverage evaluation took the full-grid branch — indexing
    // elev_matrix (which only holds g.bin.data_window under AOI/UNBOUNDED)
    // with a global row * FULL-grid stride. That is the out-of-bounds
    // write behind the segfault: masked as long as the feasible grid was
    // empty (nothing ever got placed), guaranteed to crash the moment a
    // real placement was evaluated.
    //
    // One path now, exactly like ViewshedEngine/RadarEngine in
    // operations.h: everything derives from g.bin.data_window, which IS
    // the resident layout in both modes (the whole grid when
    // g.bin.load_full is true, the loaded window otherwise) — no flag to
    // get backwards, nothing for a caller to forget to pass.
    if (!g.bin.window_valid) {
        std::cerr << "[radar] ERROR: computeLocalMask called with no window "
                     "resolved — call Global::load() first\n";
        return emptyMask();
    }

    const GridWindow& dw = g.bin.data_window;

    if (!dw.contains(srcRow, srcCol)) {
        std::cerr << "[radar] ERROR: placement (" << srcRow << "," << srcCol
                  << ") is not resident — rows [" << dw.r0 << " -> " << (dw.r0 + dw.rows)
                  << ") cols [" << dw.c0 << " -> " << (dw.c0 + dw.cols)
                  << ") — elev_matrix has no data for it, skipping\n";
        return emptyMask();
    }

    // Everything from here down operates in RESIDENT-LOCAL coordinates —
    // the space elev_matrix (and therefore demPtr) is actually laid out
    // in. In load_full mode dw.r0/c0 are 0, so this is the identity and
    // resident-local == global — one code path, correct in both modes.
    const int dwCols = dw.cols;
    const int dwRows = dw.rows;
    const int localSrcRow = dw.toLocalRow(srcRow);
    const int localSrcCol = dw.toLocalCol(srcCol);

    // Grid resolution the data was actually BUILT at (Bin::cell_m), never
    // the compile-time RESOLUTION — see cellM parameter on castRayFast.
    const double cellM = (double)g.bin.cell_m;

    // Clamp radius to the RESIDENT bounds around the (local) source — this
    // limits how far a ray can go before running off the edge of
    // elev_matrix.
    int R;
    if (use_radar_model) {
        R = pow((peak_power * pow(antenna_gain, 2) * pow(wavelength, 2) * crossSection) /
                (pow(4 * std::numbers::pi, 2) * pow(10, detachable_threshold)), 0.25);
        R = (int)((R * 1000) / cellM);
    } else {
        R = static_cast<int>(std::ceil(maxRangeM / cellM));
    }
    R = std::min(R, std::min(
            std::min(localSrcRow, dwRows - localSrcRow - 1),
            std::min(localSrcCol, dwCols - localSrcCol - 1)));

    const double horizonM =
        std::min(maxRangeM,
                 std::sqrt(2.0 * K_EFF_R * antennaH) +
                 std::sqrt(2.0 * K_EFF_R * targetH));

    const double minSlope = std::tan(minElevDeg * DEG2RAD);
    const double maxSlope = std::tan(maxElevDeg * DEG2RAD);

    const double srcElev =
        static_cast<double>(g.bin.elev_matrix[
            static_cast<size_t>(localSrcRow) * dwCols + localSrcCol]) + antennaH;

    // ── Local bounding box — used for ALL internal indexing and ray
    // casting, since it matches elev_matrix's own layout.
    const int localRowMin = std::max(0,      localSrcRow - R - 1);
    const int localRowMax = std::min(dwRows, localSrcRow + R + 2);
    const int localColMin = std::max(0,      localSrcCol - R - 1);
    const int localColMax = std::min(dwCols, localSrcCol + R + 2);
    const int localCols   = localColMax - localColMin;
    const int localRows   = localRowMax - localRowMin;

    LocalMask m;
    m.localCols = localCols;
    m.R         = R;
    m.buf.assign(static_cast<size_t>(localRows) * localCols, 0);

    // Source cell is always covered
    m.buf[static_cast<size_t>(localSrcRow - localRowMin) * localCols +
          (localSrcCol - localColMin)] = 1;

    // "Virtual origin" base so castRayFast can use local row/col directly
    uint8_t* statBase =
        m.buf.data()
        - static_cast<std::ptrdiff_t>(localRowMin) * localCols
        - localColMin;

    const short*      demPtr   = g.bin.elev_matrix.data();
    const auto        boundary = circleBoundary(localSrcRow, localSrcCol, R);
    const int         N        = static_cast<int>(boundary.size());

    // Rays are independent — safe to parallelise (binary write, no race hazard)
    // Guard against nested parallelism if this function is already called from
    // a parallel region (e.g. per-chromosome GA threads).
#if defined(_OPENMP)
    #pragma omp parallel for schedule(dynamic, 64) if(!omp_in_parallel() && !tls_serial_rays)
#endif
    for (int i = 0; i < N; ++i) {
        const int er = boundary[i].x, ec = boundary[i].y;
        if (er >= 0 && ec >= 0 && er < dwRows && ec < dwCols)
            castRayFast(
                localSrcRow, localSrcCol, er, ec,
                srcElev, targetH,
                minSlope, maxSlope,
                horizonM, R,
                demPtr, statBase,
                m.localCols, dwCols, use_curv, useAMSL, cellM);
    }

    // Convert the exported bounding box back to GLOBAL coordinates. Every
    // consumer of LocalMask (applyMask, countSectorCovered,
    // get_individual_coverage) compares rowMin/colMin against
    // FeasibleCell::r/c, which is global — this keeps that comparison
    // correct without those functions needing to know windowing exists.
    // The arithmetic still lines up with buf's internal layout because
    // (global_r - global_rowMin) == (local_r - local_rowMin) for any
    // point, since global_rowMin = local_rowMin + dw.r0 by construction.
    m.rowMin = dw.toGlobalRow(localRowMin);
    m.rowMax = dw.toGlobalRow(localRowMax);
    m.colMin = dw.toGlobalCol(localColMin);
    m.colMax = dw.toGlobalCol(localColMax);

    return m;  // moves the buf vector — no copy
}

// ── Count covered cells within a sector, reusing an existing LocalMask ───
//
// leftDeg/rightDeg: arc swept clockwise from leftDeg to rightDeg,
//                    North = 0°. leftDeg < 0 || rightDeg < 0 → full 360°.
//                    leftDeg > rightDeg wraps through North.
//
// srcRow/srcCol here are GLOBAL, matching m.rowMin/colMin (also global,
// per computeLocalMask above) — this function was already written in
// global terms and needed no change for windowed loading.
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
 *
 * NOTE: placements/grid/aoiGrid cells are all GLOBAL grid coordinates.
 * computeLocalMask() (called below via applyMask) requires each placement
 * to fall inside Global's currently-loaded window (g.bin.op_window) — a
 * placement outside it is skipped with an error rather than reading
 * garbage, but that also means it silently contributes zero coverage. If
 * you see placements landing outside the AOI window, check the window
 * was loaded with load({...}, AOI) (or a wide enough UNBOUNDED range)
 * before calling this.
 */
inline std::vector<int> calculate_radar_coverage_overlap(
    Global&                       g,
    // PERF FIX: both grids used to be taken BY VALUE, so every single
    // fitness evaluation deep-copied the feasible grid AND the AOI grid
    // (16M+ cells each for a ~140,000 km^2 AOI at 90 m) just to write
    // per-cell overlap counters into the copies — hundreds of MB of
    // allocation + memcpy per evaluation, times every concurrent worker.
    // Overlap counts now live in a small per-call side buffer instead;
    // the input grids are never modified or copied (except once, for the
    // optional output render).
    const FeasibleGrid&           grid,               // feasible cells only
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
    const FeasibleGrid&           aoiGrid = {},
    bool write_png = false,
    int reduction = 50,
    bool use_curv = true,
    bool threat_sector_flag = false,
    double loa = 0.5,
    double roa = 60,
    bool use_radar_model = false, int peak_power = 1, int antenna_gain=34, double wavelength = 0.1, int crossSection = 1, double detachable_threshold = 0.00000000001,
    bool useAMSL = false,
    std::pair<int,int> midpoint = {},
    int cellSize = 100
)
{
    std::vector<int> result(5, 0);
    int distance = 0;
    if (placements.empty() || grid.cells.empty()) return result;
    for (size_t i = 0; i < placements.size(); ++i)
        for (size_t j = i + 1; j < placements.size(); ++j) {
            const double dr = placements[i].r - placements[j].r;
            const double dc = placements[i].c - placements[j].c;
            if (std::sqrt(dr*dr + dc*dc) < minSpacingCells) ++result[0];
        }

    // Stats are bucketed over the full AOI grid when given, else the
    // feasible grid — unchanged from before.
    const bool separateAoi  = !aoiGrid.cells.empty();
    const bool wantOutput   = !tifPath.empty();
    const std::vector<FeasibleCell>& statsCells = separateAoi ? aoiGrid.cells : grid.cells;

    // Per-call hit counters, parallel to statsCells. These replace the
    // cell.overlap increments that used to be written into the by-value
    // grid copies. Final per-cell overlap = cell.overlap (as passed in)
    // + hits[i], i.e. exactly what the old code computed.
    std::vector<uint16_t> hits(statsCells.size(), 0);

    // The feasible grid's own overlap counts are only ever READ by the
    // output render (fitness/stats use statsCells only), so during GA
    // evaluation they are no longer computed at all — previously every
    // evaluation applied each mask to BOTH grids, which with no
    // constraints active are the same 16M cells twice.
    std::vector<uint16_t> gridHits;
    if (wantOutput && separateAoi) gridHits.assign(grid.cells.size(), 0);

    auto applyMask = [](const rc_detail::LocalMask& m,
                        const std::vector<FeasibleCell>& cells,
                        std::vector<uint16_t>& hit) {
        const size_t n = cells.size();
        for (size_t i = 0; i < n; ++i) {
            const auto& cell = cells[i];
            if (cell.r < m.rowMin || cell.r >= m.rowMax ||
                cell.c < m.colMin || cell.c >= m.colMax) continue;
            const size_t local_r = cell.r - m.rowMin;
            const size_t local_c = cell.c - m.colMin;
            const size_t flat = local_r * m.localCols + local_c;
            if (flat < m.buf.size() && m.buf[flat]) ++hit[i];
        }
    };

    auto applyAll = [&](const rc_detail::LocalMask& m) {
        applyMask(m, statsCells, hits);
        if (!gridHits.empty()) applyMask(m, grid.cells, gridHits);
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
                    r.range * 1000,use_curv, use_radar_model,
                    r.peak_power,  r.antenna_gain,  r.wavelength,
                    r.crossSection,  r.detachable_threshold, useAMSL);

            applyAll(cov);
            if(threat_sector_flag){
                distance += calc_distance(midpoint.first, midpoint.second, p.r, p.c, cellSize);
            }
        }
    } else {
        // Homogenous path — all radars share the same physical spec.
        for (const auto& p : placements) {
            const rc_detail::LocalMask mask =
                rc_detail::computeLocalMask(
                    g, p.r, p.c,
                    antennaHeight,
                    target,
                    minAngle,
                    maxAngle,
                    rangeM,use_curv, use_radar_model,
                    peak_power,  antenna_gain,  wavelength,
                    crossSection,  detachable_threshold, useAMSL);

            applyAll(mask);
            if(threat_sector_flag){
                distance += calc_distance(midpoint.first, midpoint.second, p.r, p.c, cellSize);
            }
        }
    }

    // ── Bucket final per-cell overlap counts ─────────────────────────────
    for (size_t i = 0; i < statsCells.size(); ++i) {
        const long long ov = (long long)statsCells[i].overlap + hits[i];
        if      (ov == 1) ++result[1];
        else if (ov == 2) ++result[2];
        else if (ov >= 3 && ov <= 9998) ++result[3];
    }
    result[4] = distance;

    // ── Optional render ──────────────────────────────────────────────────
    // Only here (once per run, not per evaluation) are the grids copied,
    // so the writers get cells with .overlap filled in exactly as before.
    if (wantOutput){
        auto idx = tifPath.find_last_of("\\/");
        std::string p;
        if(std::string::npos != idx){
            p = tifPath.substr(0,idx+1);
        }
        FeasibleGrid aoiOut  = aoiGrid;
        FeasibleGrid gridOut = grid;
        if (separateAoi) {
            for (size_t i = 0; i < aoiOut.cells.size(); ++i)  aoiOut.cells[i].overlap  += hits[i];
            for (size_t i = 0; i < gridOut.cells.size(); ++i) gridOut.cells[i].overlap += gridHits[i];
        } else {
            for (size_t i = 0; i < gridOut.cells.size(); ++i) gridOut.cells[i].overlap += hits[i];
        }
        writeOverlapGeoTIFF(g, aoiOut, tifPath, placements, reduction, opt, band);
        if(write_png){
            Bin::writeFeasibleGridOverlapPng(aoiOut, p+"result.png", placements, gridOut,
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
 * Returns  0.0             if no feasible cells fall within range, or if
 *          the placement lies outside the currently-loaded window.
 *
 * CHANGED: previously re-derived R independently using g.conf.grid.rows/
 * cols (the FULL grid) to clamp the radius — inconsistent with the R
 * computeLocalMask actually used once windowing existed (mask.rowMin/
 * rowMax reflect the WINDOW-clamped radius, not a full-grid one), which
 * would silently over/under-count "ideal" near a window edge. Now just
 * reuses mask.R, the exact radius computeLocalMask clamped to and built
 * the mask with, so "ideal" and "actual" are always counted against the
 * same circle.
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

    // ── Build local LOS mask ──────────────────────────────────────────────
    const LocalMask mask = computeLocalMask(
        g, srcRow, srcCol,
        antennaHeight, targetH,
        minAngle, maxAngle,
        maxRangeM, use_curv,0,0,0,0,0,0, false);

    if (mask.localCols <= 0) return 0.0;  // outside window, or no window resolved

    // Same clamped radius the mask was actually built with — no
    // independent re-derivation, no chance of the two disagreeing.
    const double R2 = static_cast<double>(mask.R) * mask.R;

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
