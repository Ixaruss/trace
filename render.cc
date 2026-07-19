/**
 * @file geotiff_export.cc
 *
 * GeoTIFF export for ViewshedEngine and RadarEngine outputs.
 *
 * CRS: EPSG:7755 — WGS 84 / India NSF LCC
 *   Lambert Conformal Conic (2SP), lat_0=24, lon_0=80,
 *   std_parallel_1=12.472955, std_parallel_2=35.1728044444444,
 *   false_easting=4,000,000 m, false_northing=4,000,000 m.
 *   This is the CRS the grid's x/y (elevation matrix, origin_x/origin_y
 *   in Config) are already expressed in — no reprojection happens here,
 *   we only tag the raster with it.
 *
 * Geotransform: derived directly from Config's Grid_Conf
 *   (origin_x, origin_y, resolution) — origin_x/origin_y is the
 *   top-left corner of cell (row=0, col=0), matching the convention
 *   used by Bin::getGridCoords (col grows east with +x, row grows
 *   south as y decreases), so pixelHeight is negative (north-up raster).
 *
 * TRANSPARENCY: written as a real 4-band RGBA raster (not a single-band
 * indexed palette with per-entry alpha). Palette alpha is only reliably
 * honored by GIS-aware tools; a genuine alpha band is respected by
 * QGIS, GDAL, browsers, and plain image viewers alike. "Not covered" /
 * "not visible" cells get alpha=0 so only the actual green footprint
 * shows against whatever the raster is composited over.
 *
 * Rows are streamed one at a time rather than building full-grid RGBA
 * buffers up front — this grid is ~1.07 billion cells, so 4 full-size
 * planes would mean ~4 GB of extra pixel-data memory alone.
 */

#include "operations.h"

#include <gdal_priv.h>
#include <ogr_spatialref.h>
#include <cpl_string.h>

namespace {

// ----------------------------------------------------------------------
// Hardcoded WKT-1 fallback for EPSG:7755, used only if the GDAL/PROJ
// database installed on the build machine predates this code (7755 was
// added to EPSG relatively recently — some older PROJ installs won't
// have it). Sourced from the authoritative EPSG registry definition.
// ----------------------------------------------------------------------
constexpr const char* EPSG_7755_WKT_FALLBACK =
    "PROJCS[\"WGS 84 / India NSF LCC\","
    "GEOGCS[\"WGS 84\","
    "DATUM[\"WGS_1984\","
    "SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],"
    "AUTHORITY[\"EPSG\",\"6326\"]],"
    "PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],"
    "UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],"
    "AUTHORITY[\"EPSG\",\"4326\"]],"
    "PROJECTION[\"Lambert_Conformal_Conic_2SP\"],"
    "PARAMETER[\"latitude_of_origin\",24],"
    "PARAMETER[\"central_meridian\",80],"
    "PARAMETER[\"standard_parallel_1\",12.472955],"
    "PARAMETER[\"standard_parallel_2\",35.1728044444444],"
    "PARAMETER[\"false_easting\",4000000],"
    "PARAMETER[\"false_northing\",4000000],"
    "UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]],"
    "AXIS[\"Easting\",EAST],"
    "AXIS[\"Northing\",NORTH],"
    "AUTHORITY[\"EPSG\",\"7755\"]]";

std::string epsg7755Wkt()
{
    OGRSpatialReference srs;
    srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

    if (srs.importFromEPSG(7755) == OGRERR_NONE) {
        char* wkt = nullptr;
        srs.exportToWkt(&wkt);
        std::string result(wkt ? wkt : "");
        CPLFree(wkt);
        if (!result.empty()) return result;
    }

    std::cerr << "[WARN] PROJ database missing EPSG:7755, "
                 "using bundled WKT fallback\n";
    return EPSG_7755_WKT_FALLBACK;
}

struct RGBA { uint8_t r, g, b, a; };

RGBA overlapColorAlpha(int overlap, uint8_t alpha)
{
    if (overlap <= 0) return {0, 0, 0, 0};
    switch (overlap) {
        case 1:  return {0,   200, 0,   alpha};  // green
        case 2:  return {255, 255, 0,   alpha};  // yellow
        case 3:  return {255, 165, 0,   alpha};  // orange
        case 4:  return {255, 0,   0,   alpha};  // red
        default: return {139, 0,   0,   alpha};  // dark red, 5+
    }
}
// ----------------------------------------------------------------------
// Core writer: 4-band (R,G,B,A) Byte GeoTIFF, LZW-compressed.
//
// `src` is a single-band classification buffer (cols*rows bytes,
// row-major, row 0 = north/top). Each cell's color is chosen by
// `coveredColor` if src[i] == coveredValue, else `otherColor` — this
// keeps the caller's "just show me the covered footprint" semantics
// with no per-status enumeration needed.
// ----------------------------------------------------------------------
bool writeRGBAGeoTIFF(
    const std::string& filename,
    const uint8_t*      src,
    uint8_t              coveredValue,
    RGBA                  coveredColor,
    RGBA                  otherColor,
    int                    cols,
    int                    rows,
    double                 originX,
    double                 originY,
    double                 resolution,
    int                    markerRow = -1,       // -1 = no marker
    int                    markerCol = -1,
    int                    markerRadiusCells = 8,
    RGBA                   markerColor = {255, 0, 0, 255},
    RGBA                   markerRingColor = {255, 255, 255, 255},
    int                    markerRingThicknessCells = 2)
{
    GDALAllRegister();

    GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("GTiff");
    if (!driver) {
        std::cerr << "[ERROR] GDAL GTiff driver not available\n";
        return false;
    }

    char** opts = nullptr;
    opts = CSLSetNameValue(opts, "COMPRESS",   "LZW");
    opts = CSLSetNameValue(opts, "TILED",      "YES");
    opts = CSLSetNameValue(opts, "BIGTIFF",    "IF_SAFER");
    opts = CSLSetNameValue(opts, "ALPHA",      "YES");
    opts = CSLSetNameValue(opts, "PHOTOMETRIC","RGB");

    GDALDataset* ds = driver->Create(
        filename.c_str(), cols, rows, /*bands=*/4, GDT_Byte, opts);
    CSLDestroy(opts);

    if (!ds) {
        std::cerr << "[ERROR] Failed to create GeoTIFF: " << filename << "\n";
        return false;
    }

    double gt[6] = {
        originX, resolution, 0.0,
        originY, 0.0,        -resolution
    };
    ds->SetGeoTransform(gt);
    ds->SetProjection(epsg7755Wkt().c_str());

    GDALRasterBand* rBand = ds->GetRasterBand(1);
    GDALRasterBand* gBand = ds->GetRasterBand(2);
    GDALRasterBand* bBand = ds->GetRasterBand(3);
    GDALRasterBand* aBand = ds->GetRasterBand(4);
    rBand->SetColorInterpretation(GCI_RedBand);
    gBand->SetColorInterpretation(GCI_GreenBand);
    bBand->SetColorInterpretation(GCI_BlueBand);
    aBand->SetColorInterpretation(GCI_AlphaBand);

    const bool hasMarker = (markerRow >= 0 && markerCol >= 0);
    const long long outerRSq = (long long)markerRadiusCells * markerRadiusCells;
    const long long innerRSq = (long long)(markerRadiusCells - markerRingThicknessCells)
                              * (markerRadiusCells - markerRingThicknessCells);

    // Bounded to O(cols) per iteration -- one scanline of each band at a
    // time, not the full W*H*4 image at once.
    std::vector<uint8_t> rowR(cols), rowG(cols), rowB(cols), rowA(cols);

    bool ok = true;
    for (int y = 0; y < rows && ok; ++y) {
        const uint8_t* srcRow = src + static_cast<size_t>(y) * cols;

        // Skip the marker distance-check entirely for rows far from it --
        // this loop still runs once per row regardless, so keep the
        // common (non-marker) case cheap.
        const bool rowNearMarker = hasMarker &&
            std::abs(y - markerRow) <= markerRadiusCells;

        for (int x = 0; x < cols; ++x) {
            RGBA c = (srcRow[x] == coveredValue) ? coveredColor : otherColor;

            if (rowNearMarker) {
                long long dx = x - markerCol, dy = y - markerRow;
                long long distSq = dx * dx + dy * dy;
                if (distSq <= outerRSq)
                    c = (distSq >= innerRSq) ? markerRingColor : markerColor;
            }

            rowR[x] = c.r; rowG[x] = c.g; rowB[x] = c.b; rowA[x] = c.a;
        }
        ok = ok && rBand->RasterIO(GF_Write, 0, y, cols, 1, rowR.data(), cols, 1, GDT_Byte, 0, 0) == CE_None;
        ok = ok && gBand->RasterIO(GF_Write, 0, y, cols, 1, rowG.data(), cols, 1, GDT_Byte, 0, 0) == CE_None;
        ok = ok && bBand->RasterIO(GF_Write, 0, y, cols, 1, rowB.data(), cols, 1, GDT_Byte, 0, 0) == CE_None;
        ok = ok && aBand->RasterIO(GF_Write, 0, y, cols, 1, rowA.data(), cols, 1, GDT_Byte, 0, 0) == CE_None;
    }

    GDALClose(ds);

    if (!ok) {
        std::cerr << "[ERROR] RasterIO write failed for " << filename << "\n";
        return false;
    }

    std::cout << "[INFO]   Wrote GeoTIFF: " << filename
              << " (" << cols << " x " << rows << ", RGBA)\n";
    return true;
}

// Alpha=160/255 (~63% opaque) so the basemap underneath remains visible
// through the coverage overlay. Bump toward 255 for a more solid fill,
// or lower for more see-through.
constexpr RGBA GREEN       { 0, 200, 0, 160 };
constexpr RGBA TRANSPARENT { 0,   0, 0,   0 };

} // namespace

// ----------------------------------------------------------------------
// ViewshedEngine — binary mask (0 = not visible, 1 = visible)
// ----------------------------------------------------------------------
void ViewshedEngine::writeViewshedToGeoTIFF(const std::string& filename)
{
    writeRGBAGeoTIFF(
        filename,
        viewshedMask.data(),
        /*coveredValue=*/1, GREEN, TRANSPARENT,
        gridCols, gridRows,
        static_cast<double>(g.conf.grid.origin_x),
        static_cast<double>(g.conf.grid.origin_y),
        static_cast<double>(RESOLUTION),
        /*markerRow=*/centerX, /*markerCol=*/centerY);
}

// ----------------------------------------------------------------------
// RadarEngine
// ----------------------------------------------------------------------

// Full diagnostic result (per-cell CoverageStatus) — only COVERED cells
// are shown; every other status (masked, out of range/sector, beyond
// horizon, etc.) is fully transparent.
void RadarEngine::writeCoverageToGeoTIFF(
    const RadarCoverageResult& res, const std::string& filename)
{
    static_assert(sizeof(CoverageStatus) == sizeof(uint8_t),
        "CoverageStatus must be byte-sized for the reinterpret_cast below "
        "to be safe -- check the enum's underlying type in operations.h");

    writeRGBAGeoTIFF(
        filename,
        reinterpret_cast<const uint8_t*>(res.status.data()),
        /*coveredValue=*/static_cast<uint8_t>(CoverageStatus::COVERED),
        GREEN, TRANSPARENT,
        res.gridCols, res.gridRows,
        static_cast<double>(g.conf.grid.origin_x),
        static_cast<double>(g.conf.grid.origin_y),
        static_cast<double>(RESOLUTION),
        /*markerRow=*/res.srcRow, /*markerCol=*/res.srcCol);
}

// Binary fast-mode result (0/1 covered mask) from computeFast().
void RadarEngine::writeFastToGeoTIFF(
    const std::vector<uint8_t>& fastMask, const std::string& filename)
{
    writeRGBAGeoTIFF(
        filename,
        fastMask.data(),
        /*coveredValue=*/1, GREEN, TRANSPARENT,
        gridCols, gridRows,
        static_cast<double>(g.conf.grid.origin_x),
        static_cast<double>(g.conf.grid.origin_y),
        static_cast<double>(RESOLUTION),
        /*markerRow=*/srcRow, /*markerCol=*/srcCol);
}

bool writeOverlapGeoTIFF(
    Global&                        g,
    const FeasibleGrid&            aoiGrid,
    const std::string&             filename,
    const std::vector<GridPoint>&  placements)
{
    if (aoiGrid.height <= 0 || aoiGrid.width <= 0) {
        std::cerr << "[overlap-tif] FAIL: invalid grid dimensions "
                  << aoiGrid.width << " x " << aoiGrid.height << "\n";
        return false;
    }

    const int width  = aoiGrid.width;
    const int height = aoiGrid.height;

    // Dense local overlap buffer, scoped to the AOI's own bounding box
    // only -- not the full DEM grid (same allocation philosophy as
    // rc_detail::LocalMask). Default 0 == "no coverage" == transparent.
    std::vector<uint8_t> overlapBuf(static_cast<size_t>(height) * width, 0);

    int skipped = 0;
    for (const auto& cell : aoiGrid.cells) {
        int lr = cell.r - aoiGrid.origin_r;
        int lc = cell.c - aoiGrid.origin_c;
        if (lr < 0 || lr >= height || lc < 0 || lc >= width) { ++skipped; continue; }
        overlapBuf[static_cast<size_t>(lr) * width + lc] =
            static_cast<uint8_t>(std::min(cell.overlap, 255));
    }
    if (skipped > 0)
        std::cerr << "[overlap-tif] WARN: " << skipped
                  << " AOI cells fell outside the " << width << "x" << height
                  << " AABB and were skipped\n";

    // Local geotransform origin = full-grid origin, offset by this AABB's
    // own row/col offset -- same convention as the rest of this file
    // (row 0 = north/top, col grows east).
    const double originX = static_cast<double>(g.conf.grid.origin_x)
                          + static_cast<double>(aoiGrid.origin_c) * RESOLUTION;
    const double originY = static_cast<double>(g.conf.grid.origin_y)
                          - static_cast<double>(aoiGrid.origin_r) * RESOLUTION;

    GDALAllRegister();
    GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("GTiff");
    if (!driver) {
        std::cerr << "[ERROR] GDAL GTiff driver not available\n";
        return false;
    }

    char** opts = nullptr;
    opts = CSLSetNameValue(opts, "COMPRESS",   "LZW");
    opts = CSLSetNameValue(opts, "TILED",      "YES");
    opts = CSLSetNameValue(opts, "ALPHA",      "YES");
    opts = CSLSetNameValue(opts, "PHOTOMETRIC","RGB");

    GDALDataset* ds = driver->Create(
        filename.c_str(), width, height, /*bands=*/4, GDT_Byte, opts);
    CSLDestroy(opts);

    if (!ds) {
        std::cerr << "[ERROR] Failed to create GeoTIFF: " << filename << "\n";
        return false;
    }

    double gt[6] = {
        originX, static_cast<double>(RESOLUTION), 0.0,
        originY, 0.0,                             -static_cast<double>(RESOLUTION)
    };
    ds->SetGeoTransform(gt);
    ds->SetProjection(epsg7755Wkt().c_str());

    GDALRasterBand* rBand = ds->GetRasterBand(1);
    GDALRasterBand* gBand = ds->GetRasterBand(2);
    GDALRasterBand* bBand = ds->GetRasterBand(3);
    GDALRasterBand* aBand = ds->GetRasterBand(4);
    rBand->SetColorInterpretation(GCI_RedBand);
    gBand->SetColorInterpretation(GCI_GreenBand);
    bBand->SetColorInterpretation(GCI_BlueBand);
    aBand->SetColorInterpretation(GCI_AlphaBand);

    // Marker centers in this AOI's local pixel space.
    struct LocalMarker { int lr, lc; };
    std::vector<LocalMarker> markers;
    markers.reserve(placements.size());
    for (const auto& p : placements)
        markers.push_back({p.r - aoiGrid.origin_r, p.c - aoiGrid.origin_c});

    constexpr RGBA MARKER_FILL{ 255, 0, 255, 255 };  // magenta, opaque
    constexpr RGBA MARKER_RING{ 0,   0, 0,   255 };  // black outline
    constexpr int  MARKER_RADIUS    = 8;
    constexpr int  MARKER_RING_THICK = 2;
    constexpr long long OUTER_R_SQ = (long long)MARKER_RADIUS * MARKER_RADIUS;
    constexpr long long INNER_R_SQ =
        (long long)(MARKER_RADIUS - MARKER_RING_THICK) * (MARKER_RADIUS - MARKER_RING_THICK);

    std::vector<uint8_t> rowR(width), rowG(width), rowB(width), rowA(width);

    bool ok = true;
    for (int y = 0; y < height && ok; ++y) {
        const uint8_t* srcRow = overlapBuf.data() + static_cast<size_t>(y) * width;

        for (int x = 0; x < width; ++x) {
            RGBA c = overlapColorAlpha(srcRow[x], /*alpha=*/160);
            rowR[x] = c.r; rowG[x] = c.g; rowB[x] = c.b; rowA[x] = c.a;
        }

        // Stamp markers whose vertical extent touches this row -- bounded
        // to O(placements_near_row * MARKER_RADIUS), never a full-row scan.
        for (const auto& mk : markers) {
            const int dy = y - mk.lr;
            if (std::abs(dy) > MARKER_RADIUS) continue;
            const long long dy2 = static_cast<long long>(dy) * dy;

            const int xLo = std::max(0, mk.lc - MARKER_RADIUS);
            const int xHi = std::min(width - 1, mk.lc + MARKER_RADIUS);
            for (int x = xLo; x <= xHi; ++x) {
                const long long dx = x - mk.lc;
                const long long d2 = dx * dx + dy2;
                if (d2 > OUTER_R_SQ) continue;
                const RGBA& c = (d2 >= INNER_R_SQ) ? MARKER_RING : MARKER_FILL;
                rowR[x] = c.r; rowG[x] = c.g; rowB[x] = c.b; rowA[x] = c.a;
            }
        }

        ok = ok && rBand->RasterIO(GF_Write, 0, y, width, 1, rowR.data(), width, 1, GDT_Byte, 0, 0) == CE_None;
        ok = ok && gBand->RasterIO(GF_Write, 0, y, width, 1, rowG.data(), width, 1, GDT_Byte, 0, 0) == CE_None;
        ok = ok && bBand->RasterIO(GF_Write, 0, y, width, 1, rowB.data(), width, 1, GDT_Byte, 0, 0) == CE_None;
        ok = ok && aBand->RasterIO(GF_Write, 0, y, width, 1, rowA.data(), width, 1, GDT_Byte, 0, 0) == CE_None;
    }

    GDALClose(ds);

    if (!ok) {
        std::cerr << "[ERROR] RasterIO write failed for " << filename << "\n";
        return false;
    }

    std::cout << "[overlap-tif] Wrote GeoTIFF: " << filename
              << " (" << width << " x " << height << ", RGBA, "
              << placements.size() << " marker(s))\n";
    return true;
}
