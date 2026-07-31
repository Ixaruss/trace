#include "radarfdr.h"
#include "../operations.h"
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
    double                  srcLat,
    double                  srcLon,
    double               antennaHeight,
    double               targetH,
    double               minAngle,
    double               maxAngle,
    double               maxRangeM,
    const std::string&   tifPath,
    const std::string&   pngPath, int reduction)   // ADDED — optional PNG of source + break points
{
    using namespace rc_detail;

    const int gridRows = g.conf.grid.rows;
    const int gridCols = g.conf.grid.cols;

    int srcRow, srcCol;
    g.bin.getGridCoords(srcLat, srcLon, srcRow, srcCol);

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
    if (!tifPath.empty()) {
        const int rowMin = std::max(0,        srcRow - R - 1);
        const int rowMax = std::min(gridRows, srcRow + R + 2);
        const int colMin = std::max(0,        srcCol - R - 1);
        const int colMax = std::min(gridCols, srcCol + R + 2);
        writeReverseBreakPointsGeoTIFF(
            g, tifPath,
            rowMin, colMin, colMax - colMin, rowMax - rowMin,
            srcRow, srcCol, breakPoints, reduction);
    }
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
