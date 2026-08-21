#ifndef COVERAGE_IAD_
#define COVERAGE_IAD_

/**
 * coverage.h
 *
 * Radar coverage overlap, computed directly on the full DEM.
 *
 * The feasible-grid machinery is gone. It existed because the GA needed a
 * filtered candidate set (road/rail-adjacent, water/IB/slope-free) to place
 * radars on, and that filtering lives in Bin::vec_matrix. Line-of-sight needs
 * elevation and nothing else, so this path now reads g.bin.elev_matrix only --
 * vec_matrix does not have to be built or loaded at all.
 *
 * Memory:
 *   per placement : (2R+3)^2 bytes   -- LocalMask, allocated/freed per radar
 *   accumulated   : union bbox bytes -- OverlapAccum, 1 byte per cell
 *
 * The union bbox is the tight envelope of every radar's reach, NOT the whole
 * DEM, so a couple of 10 km sites cost a few hundred KB rather than the
 * rows*cols the old AOI grid implied.
 */

#include "config.h"
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iostream>

#if defined(_OPENMP)
#  include <omp.h>
#endif

#include "../common.h"
#include "../operations.h"

struct Pt { int x, y; };

// Coverage bitmap for one radar, over its own bounding box.
// Rows [rowMin, rowMax), cols [colMin, colMax); buf is localRows*localCols.
struct LocalMask {
    int rowMin = 0, rowMax = 0;
    int colMin = 0, colMax = 0;
    int localCols = 0;
    std::vector<uint8_t> buf;

    bool empty() const { return buf.empty(); }
    int  rows()  const { return rowMax - rowMin; }
};

// Overlap counts over the union of every mask, in full-grid coordinates.
struct OverlapAccum {
    int rowMin = 0, rowMax = 0;
    int colMin = 0, colMax = 0;
    int cols = 0;
    std::vector<uint8_t> counts;

    bool empty() const { return counts.empty(); }
    int  rows()  const { return rowMax - rowMin; }
};


class Iad {

// ======================================================================
// Ray endpoints
// ======================================================================

static inline std::vector<Pt> circleBoundary(int cr, int cc, int R)
{
    std::vector<Pt> pts;
    if (R <= 0) return pts;          // was: returned 8 copies of the centre,
    pts.reserve(8 * R);              // which made every ray a no-op

    auto push8 = [&](int px, int py) {
        pts.push_back({cr+px, cc+py}); pts.push_back({cr-px, cc+py});
        pts.push_back({cr+px, cc-py}); pts.push_back({cr-px, cc-py});
        pts.push_back({cr+py, cc+px}); pts.push_back({cr-py, cc+px});
        pts.push_back({cr+py, cc-px}); pts.push_back({cr-py, cc-px});
    };

    int x = 0, y = R, d = 3 - 2 * R;
    while (y >= x) {
        push8(x, y);
        if (d > 0) { d += 4 * (x - y) + 10; --y; }
        else       { d += 4 * x + 6; }
        ++x;
    }

    std::sort(pts.begin(), pts.end(), [](const Pt& a, const Pt& b){
        return a.x < b.x || (a.x == b.x && a.y < b.y); });

    pts.erase(std::unique(pts.begin(), pts.end(), [](const Pt& a, const Pt& b){
        return a.x == b.x && a.y == b.y; }), pts.end());

    return pts;
}

// Alternative endpoint generator: uniform angular sampling rather than a
// rasterized circle. Kept for A/B-ing ray density against circleBoundary.
static inline std::vector<Pt> angularBoundary(int cr, int cc, int R)
{
    std::vector<Pt> pts;
    if (R <= 0) return pts;

    const double angStep = 1.0 / (double)R;
    const int    N       = (int)std::ceil(2.0 * M_PI / angStep) + 1;
    pts.reserve(N);

    for (int i = 0; i < N; ++i) {
        const double angle = i * angStep;
        const int er = cr + (int)std::round((double)R * std::cos(angle));
        const int ec = cc + (int)std::round((double)R * std::sin(angle));
        if (!pts.empty() && pts.back().x == er && pts.back().y == ec) continue;
        pts.push_back({er, ec});
    }

    if (pts.size() > 1 && pts.front().x == pts.back().x
                       && pts.front().y == pts.back().y)
        pts.pop_back();

    return pts;
}

// ======================================================================
// Single ray
// ======================================================================

static inline void castRayFast(
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
    uint8_t*     statCell = statBase + static_cast<std::ptrdiff_t>(srcRow) * localCols + srcCol;

    int    row = srcRow, col = srcCol;
    double maxTerrainSlope = -1e18;
    int    steps = 0;

    while (true)
    {
        if (row == endRow && col == endCol) break;

        const int  e2 = 2 * err;
        const bool mr = (e2 > -adc);
        const bool mc = (e2 <  adr);

        if (mr) { err -= adc; row += sr; demCell += demStrideR; statCell += statStrideR; }
        if (mc) { err += adr; col += sc; demCell += sc;         statCell += sc;          }

        ++steps;
        if (steps > maxRadiusCells) break;

        // True straight-line ground distance from source
        const double drow = static_cast<double>(row - srcRow);
        const double dcol = static_cast<double>(col - srcCol);
        const double dist = std::sqrt(drow * drow + dcol * dcol) * RESOLUTION;
        //
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

// ======================================================================
// One radar's mask
// ======================================================================

static inline LocalMask computeLocalMask(
    Global& g,
    int     srcRow,     int srcCol,
    double  antennaH,   double targetH,
    double  minElevDeg, double maxElevDeg,
    double  maxRangeM,
    bool    verbose = false, bool use_curv = true)
{
    const int gridCols = g.conf.grid.cols;
    const int gridRows = g.conf.grid.rows;

    LocalMask m;

    if (srcRow < 0 || srcCol < 0 || srcRow >= gridRows || srcCol >= gridCols) {
        std::cerr << "[mask] source (" << srcRow << "," << srcCol
                  << ") is outside the " << gridRows << "x" << gridCols
                  << " grid -- no coverage produced\n";
        return m;
    }
    if (g.bin.elev_matrix.size() <
        static_cast<size_t>(gridRows) * static_cast<size_t>(gridCols)) {
        std::cerr << "[mask] elevation matrix holds " << g.bin.elev_matrix.size()
                  << " cells, grid needs "
                  << static_cast<size_t>(gridRows) * gridCols
                  << " -- call Bin::init({ELEVATION}) first\n";
        return m;
    }

    // Radius is NOT clamped to the distance-to-nearest-edge any more. That
    // clamp collapsed R to 0 for any site at or near the border, which
    // silently produced a one-cell mask. Instead the window is clipped to the
    // grid and out-of-grid ray endpoints are skipped -- Bresenham between two
    // points inside a box stays inside that box, so this is still safe.
    const int R = static_cast<int>(std::ceil(maxRangeM / RESOLUTION));
    if (R <= 0) {
        std::cerr << "[mask] range " << maxRangeM << " m is under one cell ("
                  << RESOLUTION << " m) -- no coverage produced\n";
        return m;
    }

    const double horizonM =
        std::min(maxRangeM,
                 std::sqrt(2.0 * K_EFF_R * antennaH) +
                 std::sqrt(2.0 * K_EFF_R * targetH));

    const double minSlope = std::tan(minElevDeg * DEG2RAD);
    const double maxSlope = std::tan(maxElevDeg * DEG2RAD);

    const double srcElev =
        static_cast<double>(g.bin.elev_matrix[
            static_cast<size_t>(srcRow) * gridCols + srcCol]) + antennaH;

    // Bounding box, clipped to the grid
    m.rowMin    = std::max(0,        srcRow - R - 1);
    m.rowMax    = std::min(gridRows, srcRow + R + 2);
    m.colMin    = std::max(0,        srcCol - R - 1);
    m.colMax    = std::min(gridCols, srcCol + R + 2);
    m.localCols = m.colMax - m.colMin;
    const int localRows = m.rowMax - m.rowMin;

    if (localRows <= 0 || m.localCols <= 0) {
        std::cerr << "[mask] degenerate window for source ("
                  << srcRow << "," << srcCol << ")\n";
        m = LocalMask{};
        return m;
    }

    m.buf.assign(static_cast<size_t>(localRows) * m.localCols, 0);

    // Source cell is always covered
    m.buf[static_cast<size_t>(srcRow - m.rowMin) * m.localCols +
          (srcCol - m.colMin)] = 1;

    uint8_t* statBase =
        m.buf.data()
        - static_cast<std::ptrdiff_t>(m.rowMin) * m.localCols
        - m.colMin;

    const short* demPtr   = g.bin.elev_matrix.data();
    const auto   boundary = circleBoundary(srcRow, srcCol, R);
    const int    N        = static_cast<int>(boundary.size());

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

    if (verbose) {
        const size_t set = std::count(m.buf.begin(), m.buf.end(), uint8_t(1));
        std::cout << "[mask] src=(" << srcRow << "," << srcCol << ") R=" << R
                  << " win=[" << m.rowMin << "," << m.rowMax << ")x["
                  << m.colMin << "," << m.colMax << ")"
                  << " rays=" << N
                  << " horizon=" << horizonM << "m"
                  << " covered=" << set << " cells\n";
        if (set <= 1)
            std::cerr << "[mask] WARN: only the source cell is covered -- the "
                         "marker will hide it and the raster will look empty\n";
    }

    return m;  // moves buf -- no copy
}

// ======================================================================
// Accumulation across radars
// ======================================================================

// Union envelope of every placement's reach, clipped to the grid.
static inline OverlapAccum makeAccum(
    Global& g,
    const std::vector<GridPoint>&    placements,
    const std::vector<IadRadarConf>& radars)
{
    OverlapAccum a;
    const int gridRows = g.conf.grid.rows;
    const int gridCols = g.conf.grid.cols;

    int rMin = gridRows, rMax = -1, cMin = gridCols, cMax = -1;

    for (size_t i = 0; i < placements.size(); ++i) {
        const auto& p = placements[i];
        if (p.r < 0 || p.c < 0 || p.r >= gridRows || p.c >= gridCols) continue;

        const double rangeM = static_cast<double>(radars[i].range) * 1000.0;
        const int    R      = static_cast<int>(std::ceil(rangeM / RESOLUTION));

        rMin = std::min(rMin, std::max(0,            p.r - R - 1));
        rMax = std::max(rMax, std::min(gridRows - 1, p.r + R + 1));
        cMin = std::min(cMin, std::max(0,            p.c - R - 1));
        cMax = std::max(cMax, std::min(gridCols - 1, p.c + R + 1));
    }

    if (rMax < rMin || cMax < cMin) return a;   // nothing usable

    a.rowMin = rMin; a.rowMax = rMax + 1;
    a.colMin = cMin; a.colMax = cMax + 1;
    a.cols   = a.colMax - a.colMin;
    a.counts.assign(static_cast<size_t>(a.rows()) * a.cols, 0);
    return a;
}

static inline void accumulate(OverlapAccum& a, const LocalMask& m)
{
    if (a.empty() || m.empty()) return;

    const int r0 = std::max(a.rowMin, m.rowMin), r1 = std::min(a.rowMax, m.rowMax);
    const int c0 = std::max(a.colMin, m.colMin), c1 = std::min(a.colMax, m.colMax);
    if (r0 >= r1 || c0 >= c1) return;

    for (int r = r0; r < r1; ++r) {
        const uint8_t* src = m.buf.data()
                           + static_cast<size_t>(r - m.rowMin) * m.localCols
                           - m.colMin;
        uint8_t* dst = a.counts.data()
                     + static_cast<size_t>(r - a.rowMin) * a.cols
                     - a.colMin;
        for (int c = c0; c < c1; ++c)
            if (src[c] && dst[c] < 255) ++dst[c];
    }
}

// Thin adapter: writeOverlapGeoTIFF still speaks FeasibleGrid, so hand it one
// built straight from the accumulator. Only covered cells are emitted -- the
// writer zero-fills the rest, and zero is NoData in the output palette.
static inline FeasibleGrid toRenderGrid(const OverlapAccum& a)
{
    FeasibleGrid fg;
    fg.origin_r = a.rowMin;
    fg.origin_c = a.colMin;
    fg.height   = a.rows();
    fg.width    = a.cols;
    if (a.empty()) return fg;

    size_t hits = 0;
    for (uint8_t v : a.counts) if (v) ++hits;
    fg.cells.reserve(hits);

    for (int r = a.rowMin; r < a.rowMax; ++r) {
        const uint8_t* row = a.counts.data()
                           + static_cast<size_t>(r - a.rowMin) * a.cols;
        for (int c = a.colMin; c < a.colMax; ++c) {
            const uint8_t v = row[c - a.colMin];
            if (v) fg.cells.push_back({r, c, uint8_t(1), static_cast<int>(v)});
        }
    }
    return fg;
}

// ======================================================================
// Driver
// ======================================================================
//
// result[0] unused, result[1..3] = cells covered by exactly 1, exactly 2,
// and 3-or-more radars. Same contract as before, counted off the dense
// accumulator instead of a cell list.

static inline std::vector<int> calculate_radar_coverage_multi(
    Global&                          g,
    const std::vector<GridPoint>&    placements,
    const std::vector<IadRadarConf>& radars,
    OutputCrs opt, BandLayout band,
    int                              target    = 2000,
    const std::string&               tifPath   = "",
    int                              reduction = 0,
    bool                             verbose   = true, bool use_curv = true)
{
    std::vector<int> result(4, 0);

    if (placements.empty()) {
        std::cerr << "[coverage] no placements\n";
        return result;
    }
    if (radars.size() != placements.size()) {
        std::cerr << "[coverage] " << radars.size() << " radar config(s) for "
                  << placements.size() << " placement(s) -- they must match "
                     "1:1, the caller must filter both together\n";
        return result;
    }

    OverlapAccum accum = makeAccum(g, placements, radars);
    if (accum.empty()) {
        std::cerr << "[coverage] every placement is outside the "
                  << g.conf.grid.rows << "x" << g.conf.grid.cols
                  << " grid -- nothing to compute\n";
        return result;
    }

    if (verbose)
        std::cout << "[coverage] envelope rows[" << accum.rowMin << ","
                  << accum.rowMax << ") cols[" << accum.colMin << ","
                  << accum.colMax << ") = " << accum.rows() << "x" << accum.cols
                  << " cells\n";

    for (size_t i = 0; i < placements.size(); ++i) {
        const auto& p = placements[i];
        const auto& r = radars[i];

        if (verbose)
            std::cout << "[coverage] radar " << i
                      << (r.label.empty() ? "" : " (" + r.label + ")")
                      << " at (" << p.r << "," << p.c << ")"
                      << " range=" << r.range << "km"
                      << " antenna=" << r.antennaHeight << "m"
                      << " elev=[" << r.minAngle << "," << r.maxAngle << "]deg\n";

        const LocalMask cov = computeLocalMask(
            g, p.r, p.c,
            r.antennaHeight,
            target,
            r.minAngle,
            r.maxAngle,
            static_cast<double>(r.range) * 1000.0,
            verbose, use_curv);

        accumulate(accum, cov);
    }

    for (uint8_t v : accum.counts) {
        if      (v == 1) ++result[1];
        else if (v == 2) ++result[2];
        else if (v >= 3) ++result[3];
    }

    if (verbose)
        std::cout << "[coverage] cells by overlap: 1=" << result[1]
                  << " 2=" << result[2] << " 3+=" << result[3] << "\n";

    if (!tifPath.empty()) {
        if (result[1] + result[2] + result[3] == 0) {
            std::cerr << "[coverage] nothing covered -- skipping " << tifPath
                      << " rather than writing an empty raster\n";
            return result;
        }
        const FeasibleGrid renderGrid = toRenderGrid(accum);
        writeOverlapGeoTIFF(g, renderGrid, tifPath, placements,
                            reduction, opt, band);
    }

    return result;
}


public:
    static void get_iad_coverage_tif(Global& g,OutputCrs, BandLayout, std::string output_location = "./output/iad/coverage.tif");
};

#endif
