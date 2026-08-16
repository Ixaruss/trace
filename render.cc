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
 * TRANSPARENCY: written as a real 4-band RGBA raster by default (not a
 * single-band indexed palette with per-entry alpha). Palette alpha is
 * only reliably honored by GIS-aware tools; a genuine alpha band is
 * respected by QGIS, GDAL, browsers, and plain image viewers alike.
 * "Not covered" / "not visible" cells get alpha=0 so only the actual
 * green footprint shows against whatever the raster is composited over.
 *
 * BAND LAYOUT (bands): every writer also takes an optional BandLayout
 * (Rgba / SingleBandPalette). Rgba is the default described above.
 * SingleBandPalette writes one Byte band plus a GDAL color table instead
 * of four Byte bands — about 1/4 the pixel data — but the TIFF ColorMap
 * tag carries RGB only, GDAL does not persist per-entry alpha into it.
 * So single-band transparency is expressed as NoData on index 0 (which
 * is always reserved for fully-transparent) instead of a true alpha
 * channel. QGIS and GDAL-based tools honor that NoData value; a plain
 * image viewer will not, and will show index 0's color-table entry
 * (black) instead of "nothing here." Pick SingleBandPalette only when
 * the consumer is known to be GIS-aware and the size savings matter.
 *
 * Rows are streamed one at a time rather than building full-grid
 * buffers up front — this grid is ~1.07 billion cells, so full-size
 * planes would mean multiple GB of extra pixel-data memory alone.
 *
 * DOWNSAMPLING (reductionPercent): every writer takes an optional
 * `reductionPercent` in [0, 100]. 0 (the default) writes at native
 * resolution; 50 halves the pixel count along EACH axis (so ~1/4 the
 * total pixels); 90 leaves a tenth of the cells per side. The
 * geographic extent is always preserved exactly — the output has
 * fewer, proportionally larger pixels covering the same ground, and
 * the geotransform's pixel size is scaled to match. Nothing shifts,
 * so a reduced raster still lines up with a full-resolution one.
 *
 * WINDOWED RASTERS: writeOverlapGeoTIFF and the two reverse-raycast
 * writers emit only their own sub-window of the DEM grid (the AOI's
 * AABB, or the mask's rowMin/colMin box) rather than the full grid.
 * The geotransform origin is shifted by that window's row/col offset,
 * so these still overlay the full-grid products correctly.
 *
 * OUTPUT CRS (outCrs): every writer also takes an optional OutputCrs
 * (Native7755 / Wgs84 / Both). The raster is always *rasterized* once,
 * in EPSG:7755 -- that's the CRS the grid math is already in, so it's
 * the only path that's free of resampling error. When EPSG:4326 is
 * requested, that native raster is warped (nearest-neighbour, see
 * warpToWgs84 below) into a second file rather than reprojecting the
 * grid coordinates per cell during rasterization. This keeps a single
 * source of truth for "what does this cell look like" and makes the
 * 4326 product a strict, verifiable derivative of the 7755 one. For a
 * SingleBandPalette source, the warp is nearest with NoData=0 on both
 * ends rather than an alpha-aware warp, so the "outside the rotated
 * footprint" corner wedges stay NoData instead of picking up index 0's
 * color-table entry.
 */

#include "operations.h"
#include "common/radarfdr.h"

#include <gdal_priv.h>
#include <ogr_spatialref.h>
#include <cpl_string.h>
#include <gdal_utils.h>
#include <cpl_vsi.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_map>

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

// Builds EPSG:7755 into `srs`. Returns false if neither PROJ nor the
// bundled WKT yields a usable CRS -- the caller must not then pretend
// the raster got georeferenced.
bool epsg7755Srs(OGRSpatialReference& srs)
{
    srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

    if (srs.importFromEPSG(7755) == OGRERR_NONE) return true;

    std::cerr << "[WARN] PROJ could not resolve EPSG:7755 (proj.db missing, "
                 "unset PROJ_LIB/PROJ_DATA, or too old a PROJ) -- "
                 "trying bundled WKT\n";

    srs.Clear();
    srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    if (srs.importFromWkt(EPSG_7755_WKT_FALLBACK) == OGRERR_NONE) return true;

    std::cerr << "[ERROR] Bundled EPSG:7755 WKT also failed to parse -- "
                 "raster will be written WITHOUT a CRS\n";
    return false;
}

// Applies geotransform + CRS, reporting each failure rather than
// letting an untagged raster out the door silently.
void applyGeoreferencing(GDALDataset* ds, const double gt[6],
                         const std::string& filename)
{
    if (ds->SetGeoTransform(const_cast<double*>(gt)) != CE_None)
        std::cerr << "[ERROR] SetGeoTransform failed for " << filename << "\n";

    OGRSpatialReference srs;
    if (epsg7755Srs(srs)) {
        if (ds->SetSpatialRef(&srs) != CE_None)
            std::cerr << "[ERROR] SetSpatialRef failed for " << filename
                      << " -- output will not be georeferenced\n";
    }
}

// ======================================================================
// Reprojection to EPSG:4326
// ======================================================================
//
// We rasterize once in EPSG:7755 and warp that, rather than projecting
// per cell during the write. The grid's x/y ARE 7755 metres, so the
// native product stays byte-exact and the 4326 one is a pure resample
// of it — the two can't disagree about what's covered.
//
// RESAMPLING IS NEAREST, deliberately. These are classification colors
// with a hard alpha edge: "not covered" is (0,0,0,0), so any averaging
// kernel blends RGB 0 into the footprint boundary and leaves a dark
// fringe around every coverage blob, plus half-transparent cells that
// mean nothing. Nearest keeps the palette exact.
//
// INIT_DEST=0 makes destination pixels the warp never touches fully
// transparent (alpha 0) instead of undefined — the reprojected raster
// is a rotated parallelogram inside a north-up 4326 box, so those
// corner wedges exist in every output. For a paletted source there is
// no alpha band to init, so srcnodata/dstnodata=0 does the equivalent
// job: untouched pixels come out as index 0, which is also NoData.

// Insert `suffix` before the extension: "cov.tif" -> "cov_wgs84.tif".
std::string insertSuffix(const std::string& path, const std::string& suffix)
{
    const size_t slash = path.find_last_of("/\\");
    const size_t dot   = path.find_last_of('.');
    if (dot == std::string::npos || (slash != std::string::npos && dot < slash))
        return path + suffix;                       // no extension
    return path.substr(0, dot) + suffix + path.substr(dot);
}

// `paletted` selects how the destination is told to preserve
// transparency: a real alpha band (RGBA source) or NoData=0 (single-band
// palette source, which has no alpha band to carry across the warp).
bool warpToWgs84(const std::string& srcPath, const std::string& dstPath,
                 bool paletted)
{
    GDALDatasetH src = GDALOpen(srcPath.c_str(), GA_ReadOnly);
    if (!src) {
        std::cerr << "[ERROR] Reprojection: could not reopen " << srcPath << "\n";
        return false;
    }

    // If epsg7755Srs() fell all the way through, the source carries no
    // CRS and a warp to 4326 would be a fabrication. Bail loudly and
    // leave the untagged native raster as the only output.
    const char* proj = GDALGetProjectionRef(src);
    if (!proj || proj[0] == '\0') {
        std::cerr << "[ERROR] Reprojection: " << srcPath << " has no CRS "
                     "(EPSG:7755 lookup must have failed) -- cannot warp to "
                     "EPSG:4326\n";
        GDALClose(src);
        return false;
    }

    char** argv = nullptr;
    argv = CSLAddString(argv, "-of");     argv = CSLAddString(argv, "GTiff");
    argv = CSLAddString(argv, "-t_srs");  argv = CSLAddString(argv, "EPSG:4326");
    argv = CSLAddString(argv, "-r");      argv = CSLAddString(argv, "near");
    argv = CSLAddString(argv, "-wo");     argv = CSLAddString(argv, "INIT_DEST=0");
    argv = CSLAddString(argv, "-co");     argv = CSLAddString(argv, "COMPRESS=LZW");
    argv = CSLAddString(argv, "-co");     argv = CSLAddString(argv, "TILED=YES");
    argv = CSLAddString(argv, "-co");     argv = CSLAddString(argv, "BIGTIFF=IF_SAFER");
    if (paletted) {
        // Nearest keeps the palette indices exact (no averaging garbage
        // indices into existence); 0 is NoData on both ends so the
        // untouched corner wedges stay "nothing here" rather than
        // becoming whatever colour table slot 0 happens to hold.
        argv = CSLAddString(argv, "-srcnodata"); argv = CSLAddString(argv, "0");
        argv = CSLAddString(argv, "-dstnodata"); argv = CSLAddString(argv, "0");
    } else {
        argv = CSLAddString(argv, "-co"); argv = CSLAddString(argv, "ALPHA=YES");
        argv = CSLAddString(argv, "-co"); argv = CSLAddString(argv, "PHOTOMETRIC=RGB");
    }
    argv = CSLAddString(argv, "-multi");
    argv = CSLAddString(argv, "-overwrite");

    GDALWarpAppOptions* wopts = GDALWarpAppOptionsNew(argv, nullptr);
    CSLDestroy(argv);
    if (!wopts) {
        std::cerr << "[ERROR] Reprojection: GDALWarpAppOptionsNew failed\n";
        GDALClose(src);
        return false;
    }

    int usageErr = 0;
    GDALDatasetH out = GDALWarp(dstPath.c_str(), nullptr, 1, &src, wopts, &usageErr);

    GDALWarpAppOptionsFree(wopts);
    GDALClose(src);

    if (!out) {
        std::cerr << "[ERROR] Reprojection to EPSG:4326 failed for " << dstPath
                  << " (GDAL: " << CPLGetLastErrorMsg() << ")\n";
        return false;
    }

    const int w = GDALGetRasterXSize(out), h = GDALGetRasterYSize(out);
    GDALClose(out);

    std::cout << "[INFO]   Reprojected to EPSG:4326: " << dstPath
              << " (" << w << " x " << h << ")\n";
    return true;
}

// Path planning shared by both top-level writers. The 7755 raster is
// always physically written first (it's the warp input); `nativePath`
// is a temp file when the caller only asked for 4326.
struct CrsPlan {
    std::string nativePath;   // where the 7755 raster is written
    std::string wgs84Path;    // where the 4326 raster goes (empty if none)
    bool        keepNative;
};

CrsPlan planCrsOutputs(const std::string& filename, OutputCrs outCrs)
{
    switch (outCrs) {
        case OutputCrs::Wgs84:
            return { filename + ".native7755.tmp.tif", filename, false };
        case OutputCrs::Both:
            return { filename, insertSuffix(filename, "_wgs84"), true };
        case OutputCrs::Native7755:
        default:
            return { filename, std::string(), true };
    }
}

// Run after the native raster is closed. Returns false only if a
// requested reprojection actually failed. `paletted` must match the
// band layout the native raster was actually written with.
bool finalizeCrsOutputs(const CrsPlan& plan, bool paletted)
{
    if (plan.wgs84Path.empty()) return true;

    const bool ok = warpToWgs84(plan.nativePath, plan.wgs84Path, paletted);

    if (!plan.keepNative) {
        // Caller asked for 4326 only. Keep the intermediate if the warp
        // failed -- losing the data entirely is worse than an extra file.
        if (ok) VSIUnlink(plan.nativePath.c_str());
        else    std::cerr << "[WARN] Keeping intermediate 7755 raster at "
                          << plan.nativePath << " since the warp failed\n";
    }
    return ok;
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

// ======================================================================
// Downsampling helpers
// ======================================================================

// Validate + clamp the caller's percentage. Out-of-range is a caller
// bug, so it's loud, but we don't fail the export over it.
int clampReduction(int reductionPercent, const char* who)
{
    if (reductionPercent < 0 || reductionPercent > 100) {
        std::cerr << "[WARN] " << who << ": reductionPercent "
                  << reductionPercent
                  << " outside [0,100], clamping\n";
        reductionPercent = std::min(100, std::max(0, reductionPercent));
    }
    return reductionPercent;
}

// Output pixel count along one axis. 0% -> unchanged, 50% -> half,
// 100% -> would be zero, so we floor at 1 px (a 100% reduction is
// degenerate but shouldn't produce an invalid dataset).
int reducedDim(int n, int reductionPercent)
{
    if (reductionPercent <= 0) return n;
    const double scale = (100.0 - reductionPercent) / 100.0;
    int out = static_cast<int>(std::lround(n * scale));
    return std::min(n, std::max(1, out));
}

// Half-open source index range [lo, hi) that output index `outIdx`
// aggregates over. The blocks exactly tile the source extent with no
// gaps or overlap even when outN doesn't divide srcN evenly, so no
// source cell is silently dropped and none is counted twice.
inline void srcSpan(int outIdx, int outN, int srcN, int& lo, int& hi)
{
    lo = static_cast<int>((static_cast<long long>(outIdx)     * srcN) / outN);
    hi = static_cast<int>((static_cast<long long>(outIdx) + 1) * srcN / outN);
    if (hi <= lo) hi = lo + 1;
    if (hi > srcN) hi = srcN;
}

// Scans a covered/mask buffer for the tightest [rowMin,rowMax] x
// [colMin,colMax] box containing every cell == coveredValue. Returns
// false if nothing is covered at all -- callers should skip the write
// (or fall back to a placeholder) rather than ask for a degenerate
// window.
bool computeCoverageBBox(
    const uint8_t* src, int cols, int rows, uint8_t coveredValue,
    int& rowMin, int& rowMax, int& colMin, int& colMax)
{
    rowMin = rows; rowMax = -1;
    colMin = cols; colMax = -1;

    for (int y = 0; y < rows; ++y) {
        const uint8_t* row = src + static_cast<size_t>(y) * cols;
        int lo = -1, hi = -1;
        for (int x = 0; x < cols; ++x) {
            if (row[x] == coveredValue) {
                if (lo < 0) lo = x;
                hi = x;
            }
        }
        if (lo >= 0) {
            rowMin = std::min(rowMin, y);
            rowMax = std::max(rowMax, y);
            colMin = std::min(colMin, lo);
            colMax = std::max(colMax, hi);
        }
    }
    return rowMax >= rowMin;
}

// Map a source row/col onto the reduced grid (used to keep placement
// markers on the right cell after downsampling).
inline int scaleIndex(int idx, int srcN, int outN)
{
    if (idx < 0) return -1;
    int out = static_cast<int>((static_cast<long long>(idx) * outN) / srcN);
    return std::min(outN - 1, std::max(0, out));
}

// ======================================================================
// Single-band paletted output
// ======================================================================
//
// The color decision per pixel is IDENTICAL to the RGBA path -- we just
// emit an index into a color table instead of four samples. So the two
// layouts are guaranteed to agree pixel-for-pixel on what colour a cell
// got; only the encoding differs.
//
// Index 0 is reserved for fully-transparent and tagged as NoData,
// because the TIFF ColorMap tag has no alpha channel: GDAL will write
// the RGB triplets and silently drop the alpha we set on each entry.
// NoData is the only transparency signal that survives the round trip
// (and it's also what INIT_DEST=0 / srcnodata=0 leaves behind at the
// warp step, so the two mechanisms line up).
//
// The whole palette is registered BEFORE the first RasterIO -- the GTiff
// driver requires SetColorTable() to precede pixel writes, so a lazily
// grown table would be dropped on the floor. Callers therefore pre-seed
// this with every colour that can possibly appear (markers first, then
// the fixed classification/ramp colours) rather than discovering colours
// as pixels are visited.

struct RgbaIndexer {
    std::vector<RGBA>                     lut;
    std::unordered_map<uint32_t, uint8_t> byColor;
    size_t                                overflow = 0;

    RgbaIndexer() { (void)index(RGBA{0, 0, 0, 0}); }  // slot 0 == NoData

    static uint32_t pack(RGBA c) {
        return (uint32_t(c.r) << 24) | (uint32_t(c.g) << 16) |
               (uint32_t(c.b) <<  8) |  uint32_t(c.a);
    }

    uint8_t index(RGBA c) {
        const uint32_t key = pack(c);
        auto it = byColor.find(key);
        if (it != byColor.end()) return it->second;
        if (lut.size() >= 256) { ++overflow; return nearest(c); }
        const uint8_t idx = static_cast<uint8_t>(lut.size());
        lut.push_back(c);
        byColor.emplace(key, idx);
        return idx;
    }

    // Only reachable with a >256-colour source palette. Alpha is
    // weighted hard so an opaque colour never collapses onto a
    // transparent entry (that would punch holes in the footprint,
    // which is far worse than a slightly-off hue).
    uint8_t nearest(RGBA c) const {
        uint8_t best = 0;
        long long bestD = std::numeric_limits<long long>::max();
        for (size_t i = 0; i < lut.size(); ++i) {
            const long long dr = int(lut[i].r) - c.r, dg = int(lut[i].g) - c.g;
            const long long db = int(lut[i].b) - c.b, da = int(lut[i].a) - c.a;
            const long long d  = dr*dr + dg*dg + db*db + 64*da*da;
            if (d < bestD) { bestD = d; best = static_cast<uint8_t>(i); }
        }
        return best;
    }
};

void applyColorTable(GDALRasterBand* band, const RgbaIndexer& ix,
                     const std::string& filename)
{
    GDALColorTable ct(GPI_RGB);
    for (size_t i = 0; i < ix.lut.size(); ++i) {
        const RGBA& c = ix.lut[i];
        // The 4th component is set for readers that do carry alpha
        // through a color table copy; GTiff itself ignores it, which is
        // exactly why index 0 is also marked NoData below.
        GDALColorEntry e{ short(c.r), short(c.g), short(c.b), short(c.a) };
        ct.SetColorEntry(static_cast<int>(i), &e);
    }

    band->SetColorInterpretation(GCI_PaletteIndex);
    if (band->SetColorTable(&ct) != CE_None)
        std::cerr << "[ERROR] SetColorTable failed for " << filename
                  << " -- single-band raster will have no colours\n";
    if (band->SetNoDataValue(0.0) != CE_None)
        std::cerr << "[WARN] SetNoDataValue(0) failed for " << filename
                  << " -- uncovered cells may render opaque\n";

    if (ix.overflow)
        std::cerr << "[WARN] " << filename << ": palette exceeded 256 entries ("
                  << ix.overflow << " colour(s) snapped to the nearest slot)\n";
}

// ----------------------------------------------------------------------
// Core writer: Byte GeoTIFF, LZW-compressed, in one of two band layouts
// (see BandLayout in operations.h):
//
//   BandLayout::Rgba              -- 4-band (R,G,B,A), true alpha.
//   BandLayout::SingleBandPalette -- 1-band + GDAL color table, NoData=0
//       standing in for alpha (see the "Single-band paletted output"
//       section above for why).
//
// `src` is a single-band classification buffer (cols*rows bytes,
// row-major, row 0 = north/top). Two coloring modes, independent of the
// band layout above:
//
//   palette == nullptr  -- each cell is `coveredColor` if it equals
//       coveredValue, else `otherColor`. Downsampling rule: a block is
//       covered if ANY source cell in it is. Used by the status/mask
//       writers, where the non-covered values are unordered enum tags
//       and taking a max over them would be meaningless.
//
//   palette != nullptr  -- 256-entry lookup, cell value indexes it
//       directly. Downsampling rule: MAX value in the block wins, so
//       higher palette indices are drawn on top when cells collide.
//       coveredValue/coveredColor/otherColor are ignored in this mode.
//
// `reductionPercent` (0-100, 0 = native) shrinks the raster by that
// percentage along each axis. Neither rule is nearest-neighbour or
// majority, both of which erase single-cell slivers and thin ridgeline
// coverage -- exactly the detail worth keeping when you're zoomed out
// far enough to need a reduced raster in the first place.
//
// `outCrs` controls what gets written to disk: the native EPSG:7755
// raster (default), an EPSG:4326 warp of it, or both. See the
// "Reprojection to EPSG:4326" section above for why this is a
// warp-after-write rather than a per-cell reprojection.
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
    int                    reductionPercent = 0,  // 0-100, 0 = full res
    int                    markerRow = -1,        // -1 = no marker
    int                    markerCol = -1,
    int                    markerRadiusCells = 8,
    RGBA                   markerColor = {255, 0, 0, 255},
    RGBA                   markerRingColor = {255, 255, 255, 255},
    int                    markerRingThicknessCells = 2,
    const RGBA*            palette = nullptr,
    int                    srcStride = -1,     // 256 entries, or null
    OutputCrs              outCrs = OutputCrs::Native7755,
    BandLayout             bands = BandLayout::Rgba)
{
    if (cols <= 0 || rows <= 0) {
        std::cerr << "[ERROR] Invalid raster dimensions "
                  << cols << " x " << rows << " for " << filename << "\n";
        return false;
    }

    const int stride = (srcStride > 0) ? srcStride : cols;
    const bool paletted = (bands == BandLayout::SingleBandPalette);

    reductionPercent = clampReduction(reductionPercent, "writeRGBAGeoTIFF");
    const int outCols = reducedDim(cols, reductionPercent);
    const int outRows = reducedDim(rows, reductionPercent);

    GDALAllRegister();

    GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("GTiff");
    if (!driver) {
        std::cerr << "[ERROR] GDAL GTiff driver not available\n";
        return false;
    }

    char** opts = nullptr;
    opts = CSLSetNameValue(opts, "COMPRESS", "LZW");
    opts = CSLSetNameValue(opts, "TILED",    "YES");
    opts = CSLSetNameValue(opts, "BIGTIFF",  "IF_SAFER");
    if (paletted) {
        opts = CSLSetNameValue(opts, "PHOTOMETRIC", "PALETTE");
    } else {
        opts = CSLSetNameValue(opts, "ALPHA",       "YES");
        opts = CSLSetNameValue(opts, "PHOTOMETRIC", "RGB");
    }

    const CrsPlan plan = planCrsOutputs(filename, outCrs);


    GDALDataset* ds = driver->Create(
        plan.nativePath.c_str(), outCols, outRows,
        /*bands=*/paletted ? 1 : 4, GDT_Byte, opts);
    CSLDestroy(opts);

    if (!ds) {
        std::cerr << "[ERROR] Failed to create GeoTIFF: " << plan.nativePath << "\n";
        return false;
    }

    // Same top-left corner, same total ground extent — only the pixel
    // size grows. Derived from the ACTUAL dimension ratio rather than
    // the requested percentage, since reducedDim() rounds.
    const double pixelW = resolution * static_cast<double>(cols) / outCols;
    const double pixelH = resolution * static_cast<double>(rows) / outRows;

    const double gt[6] = {
        originX, pixelW, 0.0,
        originY, 0.0,    -pixelH
    };
    applyGeoreferencing(ds, gt, plan.nativePath);

    // Marker geometry lives in OUTPUT pixel space: only the center is
    // rescaled, the radius stays put, so the marker keeps the same
    // on-screen size instead of shrinking away with the raster. Computed
    // before band setup because the palette (if any) needs to know
    // whether marker colours will actually appear.
    const bool hasMarker = (markerRow >= 0 && markerCol >= 0 &&
                            markerRow < rows && markerCol < cols);
    const int outMarkerRow = hasMarker ? scaleIndex(markerRow, rows, outRows) : -1;
    const int outMarkerCol = hasMarker ? scaleIndex(markerCol, cols, outCols) : -1;

    const long long outerRSq = (long long)markerRadiusCells * markerRadiusCells;
    const long long innerRSq = (long long)(markerRadiusCells - markerRingThicknessCells)
                              * (markerRadiusCells - markerRingThicknessCells);

    GDALRasterBand *rBand = nullptr, *gBand = nullptr,
                   *bBand = nullptr, *aBand = nullptr, *iBand = nullptr;
    RgbaIndexer ix;

    if (paletted) {
        iBand = ds->GetRasterBand(1);
        // Register every colour that can possibly be emitted, before any
        // pixel is written. Markers first: they're few and must never be
        // the ones that get squeezed out by a large source palette.
        if (hasMarker) { ix.index(markerColor); ix.index(markerRingColor); }
        if (palette) { for (int i = 0; i < 256; ++i) ix.index(palette[i]); }
        else         { ix.index(otherColor); ix.index(coveredColor); }
        applyColorTable(iBand, ix, plan.nativePath);
    } else {
        rBand = ds->GetRasterBand(1);
        gBand = ds->GetRasterBand(2);
        bBand = ds->GetRasterBand(3);
        aBand = ds->GetRasterBand(4);
        rBand->SetColorInterpretation(GCI_RedBand);
        gBand->SetColorInterpretation(GCI_GreenBand);
        bBand->SetColorInterpretation(GCI_BlueBand);
        aBand->SetColorInterpretation(GCI_AlphaBand);
    }

    // Column blocks are identical for every row, so resolve them once
    // instead of redoing the integer math per pixel.
    std::vector<int> colLo(outCols), colHi(outCols);
    for (int x = 0; x < outCols; ++x)
        srcSpan(x, outCols, cols, colLo[x], colHi[x]);

    // Bounded to O(outCols) per iteration -- one scanline at a time, not
    // the full image at once.
    std::vector<uint8_t> rowR, rowG, rowB, rowA, rowI;
    if (paletted) rowI.resize(outCols);
    else { rowR.resize(outCols); rowG.resize(outCols);
           rowB.resize(outCols); rowA.resize(outCols); }

    bool ok = true;
    for (int y = 0; y < outRows && ok; ++y) {
        int sy0, sy1;
        srcSpan(y, outRows, rows, sy0, sy1);

        // Skip the marker distance-check entirely for rows far from it --
        // this loop still runs once per row regardless, so keep the
        // common (non-marker) case cheap.
        const bool rowNearMarker = hasMarker &&
            std::abs(y - outMarkerRow) <= markerRadiusCells;

        for (int x = 0; x < outCols; ++x) {
            RGBA c;

            if (palette) {
                // Highest value in the block wins, so break-point dots
                // survive on top of the segments they sit on.
                uint8_t best = 0;
                for (int sy = sy0; sy < sy1; ++sy) {
                    const uint8_t* srcRow = src + static_cast<size_t>(sy) * stride;
                    for (int sx = colLo[x]; sx < colHi[x]; ++sx)
                        if (srcRow[sx] > best) best = srcRow[sx];
                }
                c = palette[best];
            } else {
                // "Any covered cell in the block wins." Degenerates to a
                // plain lookup when the block is 1x1 (reduction = 0).
                bool covered = false;
                for (int sy = sy0; sy < sy1 && !covered; ++sy) {
                    const uint8_t* srcRow = src + static_cast<size_t>(sy) * stride;
                    for (int sx = colLo[x]; sx < colHi[x]; ++sx) {
                        if (srcRow[sx] == coveredValue) { covered = true; break; }
                    }
                }
                c = covered ? coveredColor : otherColor;
            }

            if (rowNearMarker) {
                long long dx = x - outMarkerCol, dy = y - outMarkerRow;
                long long distSq = dx * dx + dy * dy;
                if (distSq <= outerRSq)
                    c = (distSq >= innerRSq) ? markerRingColor : markerColor;
            }

            if (paletted) {
                rowI[x] = ix.index(c);
            } else {
                rowR[x] = c.r; rowG[x] = c.g; rowB[x] = c.b; rowA[x] = c.a;
            }
        }

        if (paletted) {
            ok = ok && iBand->RasterIO(GF_Write, 0, y, outCols, 1, rowI.data(), outCols, 1, GDT_Byte, 0, 0) == CE_None;
        } else {
            ok = ok && rBand->RasterIO(GF_Write, 0, y, outCols, 1, rowR.data(), outCols, 1, GDT_Byte, 0, 0) == CE_None;
            ok = ok && gBand->RasterIO(GF_Write, 0, y, outCols, 1, rowG.data(), outCols, 1, GDT_Byte, 0, 0) == CE_None;
            ok = ok && bBand->RasterIO(GF_Write, 0, y, outCols, 1, rowB.data(), outCols, 1, GDT_Byte, 0, 0) == CE_None;
            ok = ok && aBand->RasterIO(GF_Write, 0, y, outCols, 1, rowA.data(), outCols, 1, GDT_Byte, 0, 0) == CE_None;
        }
    }

    GDALClose(ds);

    if (!ok) {
        std::cerr << "[ERROR] RasterIO write failed for " << plan.nativePath << "\n";
        return false;
    }

    std::cout << "[INFO]   Wrote GeoTIFF: " << plan.nativePath
              << " (" << outCols << " x " << outRows
              << (paletted ? ", 1-band palette" : ", RGBA");
    if (reductionPercent > 0)
        std::cout << ", reduced " << reductionPercent << "% from "
                  << cols << " x " << rows
                  << ", pixel " << pixelW << "m";
    std::cout << ")\n";

    return finalizeCrsOutputs(plan, paletted);
}

// Alpha=160/255 (~63% opaque) so the basemap underneath remains visible
// through the coverage overlay. Bump toward 255 for a more solid fill,
// or lower for more see-through. (In SingleBandPalette mode this alpha
// value is only meaningful for the color-table entry itself; the
// per-pixel transparency is still all-or-nothing via NoData.)
constexpr RGBA GREEN       { 0, 200, 0, 160 };
constexpr RGBA TRANSPARENT { 0,   0, 0,   0 };

// Reverse-raycast palette, mirroring writeReverseBreakPointsPng's
// colors: green covered segments, red break-point dots on top. Higher
// index wins under downsampling, so the dots never get averaged away.
constexpr uint8_t V_EMPTY      = 0;
constexpr uint8_t V_SEGMENT    = 1;
constexpr uint8_t V_BREAKPOINT = 2;
constexpr RGBA    BREAK_DOT { 220, 30, 30, 255 };

} // namespace

// ----------------------------------------------------------------------
// ViewshedEngine — binary mask (0 = not visible, 1 = visible)
// ----------------------------------------------------------------------
void ViewshedEngine::writeViewshedToGeoTIFF(
    const std::string& filename, int reductionPercent, OutputCrs outCrs,
    BandLayout bands)
{
    int rowMin, rowMax, colMin, colMax;
    if (!computeCoverageBBox(viewshedMask.data(), gridCols, gridRows,
                             /*coveredValue=*/1, rowMin, rowMax, colMin, colMax)) {
        std::cerr << "[viewshed-tif] WARN: nothing visible -- skipping "
                  << filename << "\n";
        return;
    }

    // Keep the observer marker inside the window even if the origin
    // cell itself isn't flagged visible (self-occlusion edge case).
    rowMin = std::min(rowMin, centerX);
    rowMax = std::max(rowMax, centerX);
    colMin = std::min(colMin, centerY);
    colMax = std::max(colMax, centerY);

    // Pad so the marker ring isn't clipped at the window edge, then
    // clamp back into the full grid.
    constexpr int kPad = 10;
    rowMin = std::max(0,            rowMin - kPad);
    colMin = std::max(0,            colMin - kPad);
    rowMax = std::min(gridRows - 1, rowMax + kPad);
    colMax = std::min(gridCols - 1, colMax + kPad);

    const int winCols = colMax - colMin + 1;
    const int winRows = rowMax - rowMin + 1;

    const double originX = static_cast<double>(g.conf.grid.origin_x)
                          + static_cast<double>(colMin) * RESOLUTION;
    const double originY = static_cast<double>(g.conf.grid.origin_y)
                          - static_cast<double>(rowMin) * RESOLUTION;

    writeRGBAGeoTIFF(
        filename,
        viewshedMask.data() + static_cast<size_t>(rowMin) * gridCols + colMin,
        /*coveredValue=*/1, GREEN, TRANSPARENT,
        winCols, winRows,
        originX, originY,
        static_cast<double>(RESOLUTION),
        reductionPercent,
        /*markerRow=*/centerX - rowMin, /*markerCol=*/centerY - colMin,
        /*markerRadiusCells=*/8, {255,0,0,255}, {255,255,255,255}, 2,
        /*palette=*/nullptr,
        /*srcStride=*/gridCols,
        outCrs,
        bands);
}

// ----------------------------------------------------------------------
// RadarEngine
// ----------------------------------------------------------------------

// ----------------------------------------------------------------------
// RadarEngine — full diagnostic (per-cell CoverageStatus)
// ----------------------------------------------------------------------
void RadarEngine::writeCoverageToGeoTIFF(
    const RadarCoverageResult& res, const std::string& filename,
    int reductionPercent, OutputCrs outCrs, BandLayout bands)
{
    static_assert(sizeof(CoverageStatus) == sizeof(uint8_t),
        "CoverageStatus must be byte-sized for the reinterpret_cast below "
        "to be safe -- check the enum's underlying type in operations.h");

    const uint8_t* statusBuf = reinterpret_cast<const uint8_t*>(res.status.data());
    const uint8_t coveredValue = static_cast<uint8_t>(CoverageStatus::COVERED);

    int rowMin, rowMax, colMin, colMax;
    if (!computeCoverageBBox(statusBuf, res.gridCols, res.gridRows,
                             coveredValue, rowMin, rowMax, colMin, colMax)) {
        std::cerr << "[radar-tif] WARN: nothing covered -- skipping "
                  << filename << "\n";
        return;
    }

    // Keep the site marker inside the window even if the site cell
    // itself isn't flagged COVERED (e.g. masked at its own location).
    rowMin = std::min(rowMin, res.srcRow);
    rowMax = std::max(rowMax, res.srcRow);
    colMin = std::min(colMin, res.srcCol);
    colMax = std::max(colMax, res.srcCol);

    constexpr int kPad = 10;  // >= markerRadiusCells so the ring isn't clipped
    rowMin = std::max(0,               rowMin - kPad);
    colMin = std::max(0,               colMin - kPad);
    rowMax = std::min(res.gridRows - 1, rowMax + kPad);
    colMax = std::min(res.gridCols - 1, colMax + kPad);

    const int winCols = colMax - colMin + 1;
    const int winRows = rowMax - rowMin + 1;

    const double originX = static_cast<double>(g.conf.grid.origin_x)
                          + static_cast<double>(colMin) * RESOLUTION;
    const double originY = static_cast<double>(g.conf.grid.origin_y)
                          - static_cast<double>(rowMin) * RESOLUTION;

    writeRGBAGeoTIFF(
        filename,
        statusBuf + static_cast<size_t>(rowMin) * res.gridCols + colMin,
        coveredValue, GREEN, TRANSPARENT,
        winCols, winRows,
        originX, originY,
        static_cast<double>(RESOLUTION),
        reductionPercent,
        /*markerRow=*/res.srcRow - rowMin, /*markerCol=*/res.srcCol - colMin,
        /*markerRadiusCells=*/8, {255,0,0,255}, {255,255,255,255}, 2,
        /*palette=*/nullptr,
        /*srcStride=*/res.gridCols,
        outCrs,
        bands);
}

// ----------------------------------------------------------------------
// RadarEngine — binary fast-mode mask (0/1) from computeFast()
// ----------------------------------------------------------------------
void RadarEngine::writeFastToGeoTIFF(
    const std::vector<uint8_t>& fastMask, const std::string& filename,
    int reductionPercent, OutputCrs outCrs, BandLayout bands)
{
    int rowMin, rowMax, colMin, colMax;
    if (!computeCoverageBBox(fastMask.data(), gridCols, gridRows,
                             /*coveredValue=*/1, rowMin, rowMax, colMin, colMax)) {
        std::cerr << "[radar-fast-tif] WARN: nothing covered -- skipping "
                  << filename << "\n";
        return;
    }

    rowMin = std::min(rowMin, srcRow);
    rowMax = std::max(rowMax, srcRow);
    colMin = std::min(colMin, srcCol);
    colMax = std::max(colMax, srcCol);

    constexpr int kPad = 10;
    rowMin = std::max(0,            rowMin - kPad);
    colMin = std::max(0,            colMin - kPad);
    rowMax = std::min(gridRows - 1, rowMax + kPad);
    colMax = std::min(gridCols - 1, colMax + kPad);

    const int winCols = colMax - colMin + 1;
    const int winRows = rowMax - rowMin + 1;

    const double originX = static_cast<double>(g.conf.grid.origin_x)
                          + static_cast<double>(colMin) * RESOLUTION;
    const double originY = static_cast<double>(g.conf.grid.origin_y)
                          - static_cast<double>(rowMin) * RESOLUTION;

    writeRGBAGeoTIFF(
        filename,
        fastMask.data() + static_cast<size_t>(rowMin) * gridCols + colMin,
        /*coveredValue=*/1, GREEN, TRANSPARENT,
        winCols, winRows,
        originX, originY,
        static_cast<double>(RESOLUTION),
        reductionPercent,
        /*markerRow=*/srcRow - rowMin, /*markerCol=*/srcCol - colMin,
        /*markerRadiusCells=*/8, {255,0,0,255}, {255,255,255,255}, 2,
        /*palette=*/nullptr,
        /*srcStride=*/gridCols,
        outCrs,
        bands);
}

// ======================================================================
// Reverse ray cast (castRayReverse / get_coverage_reverse) exports
// ======================================================================
//
// Both of these emit a WINDOWED raster: only the [rowMin,rowMax) x
// [colMin,colMax) box the rays were cast across, not the full DEM grid.
// The geotransform origin is offset by that box's row/col so the result
// still registers against the full-grid products -- same convention as
// writeOverlapGeoTIFF.
// ----------------------------------------------------------------------

// LocalMask straight out of computeLocalMaskReverse: 0/1 covered, laid
// out over the mask's own window, row 0 = mask.rowMin = north.
// srcRow/srcCol are FULL-GRID coords (the radar site); pass -1 for both
// to omit the site marker.
bool writeReverseMaskGeoTIFF(
    Global&                     g,
    const rc_detail::LocalMask& mask,
    int                         srcRow,
    int                         srcCol,
    const std::string&          filename,
    int                         reductionPercent,
    OutputCrs                   outCrs,
    BandLayout                  bands)
{
    const int cols = mask.localCols;
    const int rows = mask.rowMax - mask.rowMin;

    if (cols <= 0 || rows <= 0) {
        std::cerr << "[reverse-tif] FAIL: invalid mask window "
                  << cols << " x " << rows << "\n";
        return false;
    }
    if (mask.buf.size() != static_cast<size_t>(rows) * cols) {
        std::cerr << "[reverse-tif] FAIL: mask buffer is "
                  << mask.buf.size() << " bytes but the window implies "
                  << static_cast<size_t>(rows) * cols
                  << " -- rowMin/rowMax/localCols disagree with buf\n";
        return false;
    }

    // Window origin = full-grid origin shifted by the mask's own offset.
    const double originX = static_cast<double>(g.conf.grid.origin_x)
                          + static_cast<double>(mask.colMin) * RESOLUTION;
    const double originY = static_cast<double>(g.conf.grid.origin_y)
                          - static_cast<double>(mask.rowMin) * RESOLUTION;

    // Site marker in window-local pixels. Guarded rather than assumed:
    // a clamped R can pull the window in, and a marker outside it would
    // otherwise get snapped onto an edge pixel and quietly lie.
    int markerRow = -1, markerCol = -1;
    if (srcRow >= 0 && srcCol >= 0) {
        const int lr = srcRow - mask.rowMin;
        const int lc = srcCol - mask.colMin;
        if (lr >= 0 && lr < rows && lc >= 0 && lc < cols) {
            markerRow = lr;
            markerCol = lc;
        } else {
            std::cerr << "[reverse-tif] WARN: source (" << srcRow << ","
                      << srcCol << ") is outside the mask window -- "
                         "marker omitted\n";
        }
    }

    return writeRGBAGeoTIFF(
        filename,
        mask.buf.data(),
        /*coveredValue=*/1, GREEN, TRANSPARENT,
        cols, rows,
        originX, originY,
        static_cast<double>(RESOLUTION),
        reductionPercent,
        markerRow, markerCol,
        /*markerRadiusCells=*/8, {255,0,0,255}, {255,255,255,255}, 2,
        /*palette=*/nullptr, /*srcStride=*/-1,
        outCrs,
        bands);
}

// Break points from get_coverage_reverse, rendered the same way
// Bin::writeReverseBreakPointsPng renders them: a green SEGMENT from
// the source out to each break point, with a red dot at the break
// point itself.
//
// The segments are what make the footprint solid. A break point is one
// cell out of ~28 million, so stamping only the points gives you a
// visually empty raster -- the filled blob you see in the PNG comes
// entirely from the source->breakpoint lines.
//
// The Bresenham stepping rule here is character-for-character the one
// in castRayReverse's path walk and in the PNG renderer, so the three
// agree on which cells a given ray touches.
//
// NOTE ON SEMANTICS: this fills every cell from the source up to and
// including the break point, matching what castRayReverse itself writes
// into its mask. Cells that are shadowed but nearer than the farthest
// visible cell are therefore drawn as covered. That's the existing
// definition of this product, not a rendering artifact -- if you want
// strictly-visible cells only, that's computeLocalMaskReverse's mask
// via writeReverseMaskGeoTIFF, which is a different (and slower) thing.
//
// Window args mirror writeReverseBreakPointsPng's, so the call site can
// hand both the same numbers:
//   (winRowMin, winColMin, winCols, winRows, srcRow, srcCol, breakPoints)
bool writeReverseBreakPointsGeoTIFF(
    Global&                       g,
    const std::string&            filename,
    int                           winRowMin,
    int                           winColMin,
    int                           winCols,
    int                           winRows,
    int                           srcRow,
    int                           srcCol,
    const std::vector<GridPoint>& breakPoints,
    int                           reductionPercent,
    int                           dotRadiusCells,
    bool                          fillToBreakPoints,
    OutputCrs                     outCrs,
    BandLayout                    bands)
{
    if (winCols <= 0 || winRows <= 0) {
        std::cerr << "[breakpoint-tif] FAIL: invalid window "
                  << winCols << " x " << winRows << "\n";
        return false;
    }
    if (dotRadiusCells < 0) dotRadiusCells = 0;

    if (breakPoints.empty())
        std::cerr << "[breakpoint-tif] WARN: no break points -- writing an "
                     "empty raster with just the source marker\n";

    // Single-band stamp buffer over the window (1 byte/cell -- the same
    // footprint LocalMask would take for the same R).
    std::vector<uint8_t> buf(static_cast<size_t>(winRows) * winCols, V_EMPTY);

    const int localSrcR = srcRow - winRowMin;
    const int localSrcC = srcCol - winColMin;
    const bool srcInWindow =
        (localSrcR >= 0 && localSrcR < winRows &&
         localSrcC >= 0 && localSrcC < winCols);

    if (!srcInWindow)
        std::cerr << "[breakpoint-tif] WARN: source (" << srcRow << ","
                  << srcCol << ") is outside the window -- segments and "
                     "marker omitted\n";

    // Highest value wins on collision, so a break-point dot is never
    // overwritten by a segment from a neighbouring ray.
    auto stamp = [&](int lr, int lc, uint8_t v) {
        if ((unsigned)lr >= (unsigned)winRows ||
            (unsigned)lc >= (unsigned)winCols) return;
        uint8_t& dst = buf[static_cast<size_t>(lr) * winCols + lc];
        if (v > dst) dst = v;
    };

    // Same stepping rule as castRayReverse's path walk.
    auto drawSegment = [&](int r0, int c0, int r1, int c1) {
        const int adr = std::abs(r1 - r0), adc = std::abs(c1 - c0);
        const int sr  = (r0 < r1) ? 1 : -1;
        const int sc  = (c0 < c1) ? 1 : -1;
        int err = adr - adc;
        int row = r0, col = c0;
        while (true) {
            stamp(row, col, V_SEGMENT);
            if (row == r1 && col == c1) break;
            const int e2 = 2 * err;
            if (e2 > -adc) { err -= adc; row += sr; }
            if (e2 <  adr) { err += adr; col += sc; }
        }
    };

    const long long dotRSq =
        static_cast<long long>(dotRadiusCells) * dotRadiusCells;

    int skipped = 0;
    for (const auto& p : breakPoints) {
        const int lr = p.r - winRowMin;
        const int lc = p.c - winColMin;
        if (lr < 0 || lr >= winRows || lc < 0 || lc >= winCols) { ++skipped; continue; }

        // 1. Covered segment, source -> break point. This is the fill.
        if (fillToBreakPoints && srcInWindow)
            drawSegment(localSrcR, localSrcC, lr, lc);

        // 2. Break-point dot on top.
        const int yLo = std::max(0,           lr - dotRadiusCells);
        const int yHi = std::min(winRows - 1, lr + dotRadiusCells);
        const int xLo = std::max(0,           lc - dotRadiusCells);
        const int xHi = std::min(winCols - 1, lc + dotRadiusCells);

        for (int y = yLo; y <= yHi; ++y) {
            const long long dy  = y - lr;
            const long long dy2 = dy * dy;
            for (int x = xLo; x <= xHi; ++x) {
                const long long dx = x - lc;
                if (dx * dx + dy2 <= dotRSq) stamp(y, x, V_BREAKPOINT);
            }
        }
    }
    if (skipped > 0)
            std::cerr << "[breakpoint-tif] WARN: " << skipped << " of "
                      << breakPoints.size() << " break points fell outside the "
                      << winCols << "x" << winRows << " window and were skipped\n";

        // Tighten to what actually got drawn (segments + dots), rather than
        // writing the caller's full search window verbatim.
        int bRowMin, bRowMax, bColMin, bColMax;
        if (!computeCoverageBBox(buf.data(), winCols, winRows, V_SEGMENT,
                                 bRowMin, bRowMax, bColMin, bColMax)) {
            // V_SEGMENT alone may miss an isolated dot with no fill segment
            // (fillToBreakPoints == false); fall back to scanning for any
            // non-empty cell in that case.
            auto anyDrawn = [&](int& r0, int& r1, int& c0, int& c1) {
                r0 = winRows; r1 = -1; c0 = winCols; c1 = -1;
                for (int y = 0; y < winRows; ++y)
                    for (int x = 0; x < winCols; ++x)
                        if (buf[static_cast<size_t>(y) * winCols + x] != V_EMPTY) {
                            r0 = std::min(r0, y); r1 = std::max(r1, y);
                            c0 = std::min(c0, x); c1 = std::max(c1, x);
                        }
                return r1 >= r0;
            };
            if (!anyDrawn(bRowMin, bRowMax, bColMin, bColMax)) {
                std::cerr << "[breakpoint-tif] WARN: nothing drawn -- skipping "
                          << filename << "\n";
                return false;
            }
        }
        if (srcInWindow) {
            bRowMin = std::min(bRowMin, localSrcR);
            bRowMax = std::max(bRowMax, localSrcR);
            bColMin = std::min(bColMin, localSrcC);
            bColMax = std::max(bColMax, localSrcC);
        }
        const int kPad = std::max(10, dotRadiusCells + 2);
        bRowMin = std::max(0,            bRowMin - kPad);
        bColMin = std::max(0,            bColMin - kPad);
        bRowMax = std::min(winRows - 1, bRowMax + kPad);
        bColMax = std::min(winCols - 1, bColMax + kPad);

        const int cropCols = bColMax - bColMin + 1;
        const int cropRows = bRowMax - bRowMin + 1;

        const double originX = static_cast<double>(g.conf.grid.origin_x)
                              + static_cast<double>(winColMin + bColMin) * RESOLUTION;
        const double originY = static_cast<double>(g.conf.grid.origin_y)
                              - static_cast<double>(winRowMin + bRowMin) * RESOLUTION;

        RGBA palette[256];
        for (auto& e : palette) e = TRANSPARENT;
        palette[V_SEGMENT]    = GREEN;
        palette[V_BREAKPOINT] = BREAK_DOT;

        const bool ok = writeRGBAGeoTIFF(
            filename,
            buf.data() + static_cast<size_t>(bRowMin) * winCols + bColMin,
            /*coveredValue=*/0, TRANSPARENT, TRANSPARENT,
            cropCols, cropRows,
            originX, originY,
            static_cast<double>(RESOLUTION),
            reductionPercent,
            srcInWindow ? localSrcR - bRowMin : -1,
            srcInWindow ? localSrcC - bColMin : -1,
            /*markerRadiusCells=*/8,
            /*markerColor=*/RGBA{255, 0, 255, 255},
            /*markerRingColor=*/RGBA{0, 0, 0, 255},
            /*markerRingThicknessCells=*/2,
            palette,
            /*srcStride=*/winCols,
            outCrs,
            bands);

    if (ok)
        std::cout << "[breakpoint-tif] "
                  << (static_cast<int>(breakPoints.size()) - skipped)
                  << " break point(s), segments "
                  << (fillToBreakPoints ? "filled" : "omitted")
                  << ", dot radius " << dotRadiusCells << " cell(s)\n";
    return ok;
}

// ======================================================================
// AOI overlap counts
// ======================================================================

bool writeOverlapGeoTIFF(
    Global&                        g,
    const FeasibleGrid&            aoiGrid,
    const std::string&             filename,
    const std::vector<GridPoint>&  placements,
    int                            reductionPercent,
    OutputCrs                      outCrs,
    BandLayout                     bands)
{
    if (aoiGrid.height <= 0 || aoiGrid.width <= 0) {
        std::cerr << "[overlap-tif] FAIL: invalid grid dimensions "
                  << aoiGrid.width << " x " << aoiGrid.height << "\n";
        return false;
    }

    const int width  = aoiGrid.width;
    const int height = aoiGrid.height;
    const bool paletted = (bands == BandLayout::SingleBandPalette);

    reductionPercent = clampReduction(reductionPercent, "writeOverlapGeoTIFF");
    const int outWidth  = reducedDim(width,  reductionPercent);
    const int outHeight = reducedDim(height, reductionPercent);

    // Dense local overlap buffer, scoped to the AOI's own bounding box
    // only -- not the full DEM grid (same allocation philosophy as
    // rc_detail::LocalMask). Default 0 == "no coverage" == transparent.
    // Still built at native resolution; the reduction happens at write
    // time so the source counts stay exact.
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
    // (row 0 = north/top, col grows east). Unaffected by downsampling:
    // the top-left corner doesn't move, only the pixel size grows.
    const double originX = static_cast<double>(g.conf.grid.origin_x)
                          + static_cast<double>(aoiGrid.origin_c) * RESOLUTION;
    const double originY = static_cast<double>(g.conf.grid.origin_y)
                          - static_cast<double>(aoiGrid.origin_r) * RESOLUTION;

    const double pixelW = static_cast<double>(RESOLUTION) * width  / outWidth;
    const double pixelH = static_cast<double>(RESOLUTION) * height / outHeight;

    GDALAllRegister();
    GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("GTiff");
    if (!driver) {
        std::cerr << "[ERROR] GDAL GTiff driver not available\n";
        return false;
    }

    char** opts = nullptr;
    opts = CSLSetNameValue(opts, "COMPRESS", "LZW");
    opts = CSLSetNameValue(opts, "TILED",    "YES");
    if (paletted) {
        opts = CSLSetNameValue(opts, "PHOTOMETRIC", "PALETTE");
    } else {
        opts = CSLSetNameValue(opts, "ALPHA",       "YES");
        opts = CSLSetNameValue(opts, "PHOTOMETRIC", "RGB");
    }

    const CrsPlan plan = planCrsOutputs(filename, outCrs);


    GDALDataset* ds = driver->Create(
        plan.nativePath.c_str(), outWidth, outHeight,
        /*bands=*/paletted ? 1 : 4, GDT_Byte, opts);
    CSLDestroy(opts);

    if (!ds) {
        std::cerr << "[ERROR] Failed to create GeoTIFF: " << plan.nativePath << "\n";
        return false;
    }

    const double gt[6] = {
        originX, pixelW, 0.0,
        originY, 0.0,    -pixelH
    };
    applyGeoreferencing(ds, gt, plan.nativePath);

    // Marker centers in this AOI's local pixel space, then rescaled into
    // the reduced raster's pixel space. Computed before band setup so a
    // paletted output can pre-register the marker colours.
    struct LocalMarker { int lr, lc; };
    std::vector<LocalMarker> markers;
    markers.reserve(placements.size());
    for (const auto& p : placements) {
        const int lr = p.r - aoiGrid.origin_r;
        const int lc = p.c - aoiGrid.origin_c;
        markers.push_back({
            scaleIndex(lr, height, outHeight),
            scaleIndex(lc, width,  outWidth)
        });
    }

    constexpr RGBA MARKER_FILL{ 255, 0, 255, 255 };  // magenta, opaque
    constexpr RGBA MARKER_RING{ 0,   0, 0,   255 };  // black outline
    constexpr int  MARKER_RADIUS    = 8;
    constexpr int  MARKER_RING_THICK = 2;
    constexpr long long OUTER_R_SQ = (long long)MARKER_RADIUS * MARKER_RADIUS;
    constexpr long long INNER_R_SQ =
        (long long)(MARKER_RADIUS - MARKER_RING_THICK) * (MARKER_RADIUS - MARKER_RING_THICK);

    GDALRasterBand *rBand = nullptr, *gBand = nullptr,
                   *bBand = nullptr, *aBand = nullptr, *iBand = nullptr;
    RgbaIndexer ix;

    if (paletted) {
        iBand = ds->GetRasterBand(1);
        // Overlap ramp is closed: 0 (transparent) + counts 1..5+ + two
        // marker colours == 8 entries, so this never overflows.
        ix.index(MARKER_FILL); ix.index(MARKER_RING);
        for (int n = 1; n <= 5; ++n) ix.index(overlapColorAlpha(n, 160));
        applyColorTable(iBand, ix, plan.nativePath);
    } else {
        rBand = ds->GetRasterBand(1);
        gBand = ds->GetRasterBand(2);
        bBand = ds->GetRasterBand(3);
        aBand = ds->GetRasterBand(4);
        rBand->SetColorInterpretation(GCI_RedBand);
        gBand->SetColorInterpretation(GCI_GreenBand);
        bBand->SetColorInterpretation(GCI_BlueBand);
        aBand->SetColorInterpretation(GCI_AlphaBand);
    }

    std::vector<int> colLo(outWidth), colHi(outWidth);
    for (int x = 0; x < outWidth; ++x)
        srcSpan(x, outWidth, width, colLo[x], colHi[x]);

    std::vector<uint8_t> rowR, rowG, rowB, rowA, rowI;
    if (paletted) rowI.resize(outWidth);
    else { rowR.resize(outWidth); rowG.resize(outWidth);
           rowB.resize(outWidth); rowA.resize(outWidth); }

    bool ok = true;
    for (int y = 0; y < outHeight && ok; ++y) {
        int sy0, sy1;
        srcSpan(y, outHeight, height, sy0, sy1);

        for (int x = 0; x < outWidth; ++x) {
            // MAX over the block, not mean: averaging would wash a
            // 5-radar hotspot down into the 1-2 range and lose the very
            // thing the color ramp exists to show. Worst case per block
            // is the honest summary here.
            uint8_t best = 0;
            for (int sy = sy0; sy < sy1; ++sy) {
                const uint8_t* srcRow = overlapBuf.data() + static_cast<size_t>(sy) * width;
                for (int sx = colLo[x]; sx < colHi[x]; ++sx)
                    if (srcRow[sx] > best) best = srcRow[sx];
            }
            RGBA c = overlapColorAlpha(best, /*alpha=*/160);
            if (paletted) rowI[x] = ix.index(c);
            else { rowR[x] = c.r; rowG[x] = c.g; rowB[x] = c.b; rowA[x] = c.a; }
        }

        // Stamp markers whose vertical extent touches this row -- bounded
        // to O(placements_near_row * MARKER_RADIUS), never a full-row scan.
        // Radius is in OUTPUT pixels, so markers stay the same visual size
        // no matter how far the raster was reduced.
        for (const auto& mk : markers) {
            const int dy = y - mk.lr;
            if (std::abs(dy) > MARKER_RADIUS) continue;
            const long long dy2 = static_cast<long long>(dy) * dy;

            const int xLo = std::max(0, mk.lc - MARKER_RADIUS);
            const int xHi = std::min(outWidth - 1, mk.lc + MARKER_RADIUS);
            for (int x = xLo; x <= xHi; ++x) {
                const long long dx = x - mk.lc;
                const long long d2 = dx * dx + dy2;
                if (d2 > OUTER_R_SQ) continue;
                const RGBA& c = (d2 >= INNER_R_SQ) ? MARKER_RING : MARKER_FILL;
                if (paletted) rowI[x] = ix.index(c);
                else { rowR[x] = c.r; rowG[x] = c.g; rowB[x] = c.b; rowA[x] = c.a; }
            }
        }

        if (paletted) {
            ok = ok && iBand->RasterIO(GF_Write, 0, y, outWidth, 1, rowI.data(), outWidth, 1, GDT_Byte, 0, 0) == CE_None;
        } else {
            ok = ok && rBand->RasterIO(GF_Write, 0, y, outWidth, 1, rowR.data(), outWidth, 1, GDT_Byte, 0, 0) == CE_None;
            ok = ok && gBand->RasterIO(GF_Write, 0, y, outWidth, 1, rowG.data(), outWidth, 1, GDT_Byte, 0, 0) == CE_None;
            ok = ok && bBand->RasterIO(GF_Write, 0, y, outWidth, 1, rowB.data(), outWidth, 1, GDT_Byte, 0, 0) == CE_None;
            ok = ok && aBand->RasterIO(GF_Write, 0, y, outWidth, 1, rowA.data(), outWidth, 1, GDT_Byte, 0, 0) == CE_None;
        }
    }

    GDALClose(ds);

    if (!ok) {
        std::cerr << "[ERROR] RasterIO write failed for " << plan.nativePath << "\n";
        return false;
    }

    std::cout << "[overlap-tif] Wrote GeoTIFF: " << plan.nativePath
              << " (" << outWidth << " x " << outHeight
              << (paletted ? ", 1-band palette, " : ", RGBA, ")
              << placements.size() << " marker(s)";
    if (reductionPercent > 0)
        std::cout << ", reduced " << reductionPercent << "% from "
                  << width << " x " << height
                  << ", pixel " << pixelW << "m";
    std::cout << ")\n";

    return finalizeCrsOutputs(plan, paletted);
}
