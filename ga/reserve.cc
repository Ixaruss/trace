#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <memory>
#include <fstream>
#include <string>
#include <algorithm>
#include <cstring>
#include <cstdint>
#include "../common.h"
#if defined(_OPENMP)
#include <omp.h>
#endif


// Radar parameters shared across all batch cells
struct BatchParams {
    double antennaHeightM   = 50.0;
    double targetHeightM    = 1000.0;
    double minElevDeg       = 10.0;
    double maxElevDeg       = 30.0;
    double maxRangeKM       = 100.0;
    int    sectorAzimuthDeg = -1;   // -1 = full 360 deg
};

// Binary file layout (packed to avoid any padding surprises)
#pragma pack(push, 1)
struct BinFileHeader {
    char     magic[4];       // "RCOV"
    uint32_t version;        // 1
    uint32_t numResults;
    int32_t  gridRows;
    int32_t  gridCols;
};

struct BinResultHeader {
    int32_t srcRow;
    int32_t srcCol;
    int32_t radiusCells;
    int32_t sectorAzimuthDeg;
};
#pragma pack(pop)

// In-memory representation after loading from disk
struct StoredCoverageResult {
    int srcRow, srcCol;
    int radiusCells;
    int sectorAzimuthDeg;
    int gridRows, gridCols;
    std::vector<uint8_t> status;   // gridRows * gridCols raw CoverageStatus bytes
};

/******************************************************************************
 * CONSTANTS & ENUMS
 ******************************************************************************/
constexpr double EARTH_R     = 6371000.0;
constexpr double INV_2EARTHR = 1.0 / (2.0 * EARTH_R);
constexpr double K_ATMOS     = 4.0 / 3.0;
constexpr double K_EFF_R     = K_ATMOS * EARTH_R;
constexpr double INV_2_KEFFR = 1.0 / (2.0 * K_EFF_R);
constexpr double DEG2RAD     = 3.14159265358979 / 180.0;
constexpr double RAD2DEG     = 180.0 / 3.14159265358979;

enum class CoverageStatus : uint8_t
{
    OUTSIDE_RANGE  = 0,
    COVERED        = 1,
    TERRAIN_MASKED = 2,
    BELOW_MIN_ELEV = 3,
    ABOVE_MAX_ELEV = 4,
    BEYOND_HORIZON = 5,
    OUTSIDE_SECTOR = 6,
};

static_assert(sizeof(CoverageStatus) == sizeof(uint8_t),
              "CoverageStatus must be 1 byte for binary I/O");

enum class LOSStatus : uint8_t
{
    VISIBLE        = 0,
    BEYOND_RANGE   = 1,
    BEYOND_HORIZON = 2,
    BELOW_MIN_ELEV = 3,
    ABOVE_MAX_ELEV = 4,
    TERRAIN_MASKED = 5,
    NO_DATA        = 6,
    OUT_OF_GRID    = 7,
};

/******************************************************************************
 * STRUCTURES
 ******************************************************************************/
struct Point2D { int x, y; };

struct LOSCell
{
    int    row, col;
    double distM;
    float  terrainElev;
    float  correctedElev;
    double terrainSlope;
    double maxSlopeSoFar;
    double targetSlope;
    bool   blocks;
};

struct RadarCoverageResult
{
    std::vector<CoverageStatus> status;
    int srcRow, srcCol;
    int radiusCells;
    int gridCols, gridRows;
    int sectorAzimuthDeg;
    int countCovered, countTerrainMasked;
    int countBelowMinElev, countAboveMaxElev;
    int countBeyondHorizon, countOutsideRange;
    int countOutsideSector;
};

/******************************************************************************
 * RADAR ENGINE CLASS DECLARATION
 ******************************************************************************/
class RadarEngine
{
public:
    Global& g;
    int     gridCols, gridRows;
    int     srcRow,   srcCol;
    int     maxRadiusCells;
    int     rowMin, rowMax, colMin, colMax;

    double  antennaHeight, targetHeight;
    double  minSlopeRatio, maxSlopeRatio;
    double  maxRangeM, horizonM, srcElev;

    std::vector<uint8_t> statusBuf;
    std::vector<uint8_t> fastBuf;
    std::vector<Point2D> boundary;

    std::vector<Point2D> sectorBoundary(int sectorAz) const;
    void prefill(int sectorAz, uint8_t* buf, uint8_t fillInside, uint8_t fillOutside);
    static std::vector<Point2D> circleBoundary(int cr, int cc, int radius);

    static inline double cellAzimuth(int srcRow, int srcCol, int r, int c)
    {
        double drow = r - srcRow;
        double dcol = c - srcCol;
        double az = std::atan2(dcol, -drow) * RAD2DEG;
        if (az < 0.0) az += 360.0;
        return az;
    }

    static inline bool inSectorAz(double az, double centerAz, double halfWidth)
    {
        double diff = std::fabs(az - centerAz);
        if (diff > 180.0) diff = 360.0 - diff;
        return diff <= halfWidth;
    }

    static int snapTo45(int deg)
    {
        deg = ((deg % 360) + 360) % 360;
        return (int)(std::round(deg / 45.0) * 45) % 360;
    }

    // Shared init body (called by both constructors after srcRow/srcCol are set)
    void initEngine(double antennaHeightM, double targetHeightM,
                    double minElevDeg, double maxElevDeg, double maxRangeKM);

public:
    // --- Original constructor: resolves lat/lon → grid coords via BinData ---
    RadarEngine(
        Global& globalInstance,
        double srcLat, double srcLon,
        double antennaHeightM,
        double targetHeightM,
        double minElevDeg, double maxElevDeg,
        double maxRangeKM);

    // --- NEW constructor: accepts grid coords directly (no lat/lon conversion) ---
    RadarEngine(
        Global& globalInstance,
        int srcRow, int srcCol,
        double antennaHeightM,
        double targetHeightM,
        double minElevDeg, double maxElevDeg,
        double maxRangeKM);

    RadarCoverageResult  compute(int sectorAzimuthDeg = -1);
    std::vector<uint8_t> computeFast(int sectorAzimuthDeg = -1);
    void writeToSVG(const RadarCoverageResult& res,
                    const std::string& filename = "radar_coverage.svg");
};

/******************************************************************************
 * CORE RAY CASTING FUNCTIONS (unchanged)
 ******************************************************************************/
static void castRadarRayDiag(
    int srcRow, int srcCol,
    int endRow, int endCol,
    double srcElev, double targetHeight,
    double minSlope, double maxSlope,
    double maxRangeM, double horizonM,
    int maxRadiusCells,
    const short* __restrict__ dem,
    uint8_t*     __restrict__ stat,
    int COLS, int ROWS)
{
    int adr = std::abs(endRow - srcRow), adc = std::abs(endCol - srcCol);
    int sr  = (srcRow < endRow) ? 1 : -1, sc = (srcCol < endCol) ? 1 : -1;
    int err = adr - adc;

    const std::ptrdiff_t strideR = static_cast<std::ptrdiff_t>(sr) * COLS;
    const std::ptrdiff_t strideC = sc;
    std::ptrdiff_t off = static_cast<std::ptrdiff_t>(srcRow) * COLS + srcCol;

    const short* demCell  = dem  + off;
    uint8_t*     statCell = stat + off;

    int    row = srcRow, col = srcCol;
    double dist = 0.0, curvDrop = 0.0;
    double maxTerrainSlope = -1e18;
    int    cellCount = 0;

    while (true)
    {
        if (row == endRow && col == endCol) break;

        int  e2 = 2 * err;
        bool mr = (e2 > -adc), mc = (e2 < adr);
        if (mr) { err -= adc; row += sr; demCell += strideR; statCell += strideR; }
        if (mc) { err += adr; col += sc; demCell += strideC; statCell += strideC; }

        double step = (mr && mc) ? RESOLUTION * 1.41421356237 : RESOLUTION;
        curvDrop += (2.0 * dist * step + step * step) * INV_2_KEFFR;
        dist += step;
        ++cellCount;

        if (cellCount > maxRadiusCells)                                   break;
        if (dist > maxRangeM) { *statCell = (uint8_t)CoverageStatus::OUTSIDE_RANGE; continue; }
        if (dist > horizonM)                                              break;

        short raw = *demCell;
        if (raw < -9000)                                                  break;

        double terrainElev  = static_cast<double>(raw) - curvDrop;
        double terrainSlope = (terrainElev - srcElev) / dist;
        double targetSlope  = terrainSlope + targetHeight / dist;

        uint8_t s;
        if      (maxTerrainSlope > targetSlope) s = (uint8_t)CoverageStatus::TERRAIN_MASKED;
        else if (targetSlope < minSlope)        s = (uint8_t)CoverageStatus::BELOW_MIN_ELEV;
        else if (targetSlope > maxSlope)        s = (uint8_t)CoverageStatus::ABOVE_MAX_ELEV;
        else                                    s = (uint8_t)CoverageStatus::COVERED;

        *statCell = s;
        if (terrainSlope > maxTerrainSlope) maxTerrainSlope = terrainSlope;
    }
}

static void castRadarRayFast(
    int srcRow, int srcCol,
    int endRow, int endCol,
    double srcElev, double targetHeight,
    double minSlope, double maxSlope,
    double horizonM, int maxRadiusCells,
    const short* __restrict__ dem,
    uint8_t*     __restrict__ stat,
    int COLS, int ROWS)
{
    int adr = std::abs(endRow - srcRow), adc = std::abs(endCol - srcCol);
    int sr  = (srcRow < endRow) ? 1 : -1, sc = (srcCol < endCol) ? 1 : -1;
    int err = adr - adc;

    const std::ptrdiff_t strideR = static_cast<std::ptrdiff_t>(sr) * COLS;
    const std::ptrdiff_t strideC = sc;
    std::ptrdiff_t off = static_cast<std::ptrdiff_t>(srcRow) * COLS + srcCol;

    const short* demCell  = dem  + off;
    uint8_t*     statCell = stat + off;

    int    row = srcRow, col = srcCol;
    double dist = 0.0, curvDrop = 0.0;
    double maxTerrainSlope = -1e18;
    int    cellCount = 0;

    while (true)
    {
        if (row == endRow && col == endCol) break;

        int  e2 = 2 * err;
        bool mr = (e2 > -adc), mc = (e2 < adr);
        if (mr) { err -= adc; row += sr; demCell += strideR; statCell += strideR; }
        if (mc) { err += adr; col += sc; demCell += strideC; statCell += strideC; }

        double step = (mr && mc) ? RESOLUTION * 1.41421356237 : RESOLUTION;
        curvDrop += (2.0 * dist * step + step * step) * INV_2_KEFFR;
        dist += step;
        ++cellCount;

        if (cellCount > maxRadiusCells) break;
        if (dist > horizonM)            break;

        short raw = *demCell;
        if (raw < -9000)                break;

        double terrainElev  = static_cast<double>(raw) - curvDrop;
        double terrainSlope = (terrainElev - srcElev) / dist;
        double targetSlope  = terrainSlope + targetHeight / dist;

        if (maxTerrainSlope <= targetSlope && targetSlope >= minSlope && targetSlope <= maxSlope)
            *statCell = 1;

        if (terrainSlope > maxTerrainSlope) maxTerrainSlope = terrainSlope;
    }
}

/******************************************************************************
 * RADAR ENGINE — shared init body
 ******************************************************************************/
void RadarEngine::initEngine(
    double antennaHeightM, double targetHeightM,
    double minElevDeg,     double maxElevDeg,
    double maxRangeKM)
{
    antennaHeight = antennaHeightM;
    targetHeight  = targetHeightM;
    maxRangeM     = maxRangeKM * 1000.0;

    gridCols = g.conf.grid.cols;
    gridRows = g.conf.grid.rows;

    maxRadiusCells = (int)std::ceil(maxRangeM / RESOLUTION);
    int maxPossible = std::min(
        std::min(srcRow, gridRows - srcRow - 1),
        std::min(srcCol, gridCols - srcCol - 1));
    if (maxRadiusCells > maxPossible) {
        std::cout << "[WARN] Radius clamped from " << maxRadiusCells
                  << " to " << maxPossible << "\n";
        maxRadiusCells = maxPossible;
    }

    rowMin = std::max(0,        srcRow - maxRadiusCells - 1);
    rowMax = std::min(gridRows, srcRow + maxRadiusCells + 2);
    colMin = std::max(0,        srcCol - maxRadiusCells - 1);
    colMax = std::min(gridCols, srcCol + maxRadiusCells + 2);

    statusBuf.assign((size_t)gridRows * gridCols, 0);
    fastBuf.assign(  (size_t)gridRows * gridCols, 0);

    horizonM = std::sqrt(2.0 * K_EFF_R * antennaHeight)
             + std::sqrt(2.0 * K_EFF_R * targetHeight);
    if (horizonM > maxRangeM) horizonM = maxRangeM;

    minSlopeRatio = std::tan(minElevDeg * DEG2RAD);
    maxSlopeRatio = std::tan(maxElevDeg * DEG2RAD);

    short srcTerrain = g.bin.elev_matrix[(size_t)srcRow * gridCols + srcCol];
    srcElev = static_cast<double>(srcTerrain) + antennaHeight;

    boundary = circleBoundary(srcRow, srcCol, maxRadiusCells);
}

/******************************************************************************
 * RADAR ENGINE — constructor: lat/lon path (original, unchanged)
 ******************************************************************************/
RadarEngine::RadarEngine(
    Global& globalInstance,
    double srcLat, double srcLon,
    double antennaHeightM,
    double targetHeightM,
    double minElevDeg, double maxElevDeg,
    double maxRangeKM)
    : g(globalInstance)
{
    auto t0 = std::chrono::high_resolution_clock::now();

    if (!g.bin.getGridCoords(srcLat, srcLon, srcRow, srcCol)) {
        std::cerr << "[ERROR] Source outside grid\n";
        srcRow = g.conf.grid.rows / 2;
        srcCol = g.conf.grid.cols / 2;
    }

    initEngine(antennaHeightM, targetHeightM, minElevDeg, maxElevDeg, maxRangeKM);

    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "[RADAR] Init (latlon) : "
              << std::chrono::duration<double, std::milli>(t1 - t0).count() << " ms\n";
    std::cout << "[RADAR]   Src cell    : (" << srcRow << "," << srcCol << ")\n";
    std::cout << "[RADAR]   Src terrain : "
              << g.bin.elev_matrix[(size_t)srcRow * gridCols + srcCol]
              << "m  elev: " << srcElev << "m\n";
    std::cout << "[RADAR]   Horizon     : " << horizonM / 1000.0 << " km\n";
    std::cout << "[RADAR]   Full rays   : " << boundary.size() << "\n";
}

/******************************************************************************
 * RADAR ENGINE — NEW constructor: direct grid coords (no lat/lon conversion)
 ******************************************************************************/
RadarEngine::RadarEngine(
    Global& globalInstance,
    int     srcRow_,  int srcCol_,
    double  antennaHeightM,
    double  targetHeightM,
    double  minElevDeg, double maxElevDeg,
    double  maxRangeKM)
    : g(globalInstance)
{
    auto t0 = std::chrono::high_resolution_clock::now();

    srcRow = srcRow_;
    srcCol = srcCol_;

    // Clamp to grid if out of bounds
    const int gr = g.conf.grid.rows, gc = g.conf.grid.cols;
    if (srcRow < 0 || srcRow >= gr || srcCol < 0 || srcCol >= gc) {
        std::cerr << "[ERROR] Grid cell (" << srcRow << "," << srcCol
                  << ") outside [0," << gr << ")x[0," << gc << ") — clamping to center\n";
        srcRow = gr / 2;
        srcCol = gc / 2;
    }

    initEngine(antennaHeightM, targetHeightM, minElevDeg, maxElevDeg, maxRangeKM);

    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "[RADAR] Init (grid)   : "
              << std::chrono::duration<double, std::milli>(t1 - t0).count() << " ms\n";
    std::cout << "[RADAR]   Src cell    : (" << srcRow << "," << srcCol << ")\n";
    std::cout << "[RADAR]   Src terrain : "
              << g.bin.elev_matrix[(size_t)srcRow * gridCols + srcCol]
              << "m  elev: " << srcElev << "m\n";
    std::cout << "[RADAR]   Horizon     : " << horizonM / 1000.0 << " km\n";
    std::cout << "[RADAR]   Full rays   : " << boundary.size() << "\n";
}

/******************************************************************************
 * RADAR ENGINE — circleBoundary / sectorBoundary / prefill (unchanged)
 ******************************************************************************/
std::vector<Point2D> RadarEngine::circleBoundary(int cr, int cc, int radius)
{
    std::vector<Point2D> pts;
    pts.reserve(8 * radius);
    int x = 0, y = radius, d = 3 - 2 * radius;
    auto push8 = [&](int px, int py) {
        pts.push_back({cr + px, cc + py}); pts.push_back({cr - px, cc + py});
        pts.push_back({cr + px, cc - py}); pts.push_back({cr - px, cc - py});
        pts.push_back({cr + py, cc + px}); pts.push_back({cr - py, cc + px});
        pts.push_back({cr + py, cc - px}); pts.push_back({cr - py, cc - px});
    };
    while (y >= x) {
        push8(x, y); x++;
        if (d > 0) { y--; d += 4 * (x - y) + 10; } else { d += 4 * x + 6; }
    }
    std::sort(pts.begin(), pts.end(), [](const Point2D& a, const Point2D& b) {
        return a.x < b.x || (a.x == b.x && a.y < b.y);
    });
    pts.erase(std::unique(pts.begin(), pts.end(),
        [](const Point2D& a, const Point2D& b) { return a.x == b.x && a.y == b.y; }),
        pts.end());
    return pts;
}

std::vector<Point2D> RadarEngine::sectorBoundary(int sectorAz) const
{
    if (sectorAz < 0) return boundary;
    const double halfWidth = 22.5;
    std::vector<Point2D> out;
    out.reserve(boundary.size() / 8 + 1);
    for (const auto& pt : boundary) {
        double az = cellAzimuth(srcRow, srcCol, pt.x, pt.y);
        if (inSectorAz(az, (double)sectorAz, halfWidth))
            out.push_back(pt);
    }
    return out;
}

void RadarEngine::prefill(int sectorAz, uint8_t* buf, uint8_t fillInside, uint8_t fillOutside)
{
    const long long rSq = (long long)maxRadiusCells * maxRadiusCells;
    const double halfWidth = 22.5;

    #pragma omp parallel for schedule(static, 32)
    for (int r = rowMin; r < rowMax; r++) {
        int dr = r - srcRow;
        uint8_t* rowPtr = buf + (size_t)r * gridCols + colMin;
        for (int c = colMin; c < colMax; c++) {
            int dc = c - srcCol;
            if ((long long)dr * dr + (long long)dc * dc > rSq) {
                rowPtr[c - colMin] = (uint8_t)CoverageStatus::OUTSIDE_RANGE;
                continue;
            }
            if (sectorAz < 0) {
                rowPtr[c - colMin] = fillInside;
            } else {
                double az = cellAzimuth(srcRow, srcCol, r, c);
                rowPtr[c - colMin] = inSectorAz(az, (double)sectorAz, halfWidth)
                                     ? fillInside : fillOutside;
            }
        }
    }
}

/******************************************************************************
 * RADAR ENGINE — compute (unchanged)
 ******************************************************************************/
RadarCoverageResult RadarEngine::compute(int sectorAzimuthDeg)
{
    int sectorAz = (sectorAzimuthDeg < 0) ? -1 : snapTo45(sectorAzimuthDeg);
    auto activeBoundary = sectorBoundary(sectorAz);

    std::cout << "[RADAR] Sector        : ";
    if (sectorAz < 0) std::cout << "360deg (full)\n";
    else std::cout << sectorAz << "deg +/-22.5deg  rays=" << activeBoundary.size() << "\n";

    auto t0 = std::chrono::high_resolution_clock::now();

    prefill(sectorAz, statusBuf.data(),
            (uint8_t)CoverageStatus::BEYOND_HORIZON,
            (uint8_t)CoverageStatus::OUTSIDE_SECTOR);

    statusBuf[(size_t)srcRow * gridCols + srcCol] = (uint8_t)CoverageStatus::COVERED;

    const size_t N      = activeBoundary.size();
    const short* demPtr = g.bin.elev_matrix.data();
    uint8_t*  statPtr   = statusBuf.data();
    const int COLS = gridCols, ROWS = gridRows;
    const double SE = srcElev, TH = targetHeight;
    const double MNS = minSlopeRatio, MXS = maxSlopeRatio;
    const double MXR = maxRangeM, HOR = horizonM;
    const int    MRC = maxRadiusCells;

    #pragma omp parallel for schedule(dynamic, 64) default(none) \
        shared(demPtr, statPtr, activeBoundary) \
        firstprivate(N, COLS, ROWS, SE, TH, MNS, MXS, MXR, HOR, MRC)
    for (size_t i = 0; i < N; ++i) {
        int er = activeBoundary[i].x, ec = activeBoundary[i].y;
        if (er >= 0 && ec >= 0 && er < ROWS && ec < COLS)
            castRadarRayDiag(srcRow, srcCol, er, ec,
                             SE, TH, MNS, MXS, MXR, HOR, MRC,
                             demPtr, statPtr, COLS, ROWS);
    }

    auto t1 = std::chrono::high_resolution_clock::now();

    RadarCoverageResult res;
    res.srcRow = srcRow; res.srcCol = srcCol;
    res.radiusCells = maxRadiusCells;
    res.gridCols = gridCols; res.gridRows = gridRows;
    res.sectorAzimuthDeg = sectorAz;
    res.countCovered = res.countTerrainMasked = 0;
    res.countBelowMinElev = res.countAboveMaxElev = 0;
    res.countBeyondHorizon = res.countOutsideRange = res.countOutsideSector = 0;

    long long rSq = (long long)maxRadiusCells * maxRadiusCells;
    for (int r = rowMin; r < rowMax; r++) {
        for (int c = colMin; c < colMax; c++) {
            int dr = r - srcRow, dc = c - srcCol;
            if ((long long)dr * dr + (long long)dc * dc > rSq) continue;
            switch ((CoverageStatus)statusBuf[(size_t)r * gridCols + c]) {
                case CoverageStatus::COVERED:        res.countCovered++;        break;
                case CoverageStatus::TERRAIN_MASKED: res.countTerrainMasked++;  break;
                case CoverageStatus::BELOW_MIN_ELEV: res.countBelowMinElev++;   break;
                case CoverageStatus::ABOVE_MAX_ELEV: res.countAboveMaxElev++;   break;
                case CoverageStatus::BEYOND_HORIZON: res.countBeyondHorizon++;  break;
                case CoverageStatus::OUTSIDE_SECTOR: res.countOutsideSector++;  break;
                default:                             res.countOutsideRange++;   break;
            }
        }
    }

    res.status.resize((size_t)gridRows * gridCols);
    std::memcpy(res.status.data(), statusBuf.data(), (size_t)gridRows * gridCols);

    int total = res.countCovered + res.countTerrainMasked +
                res.countBelowMinElev + res.countAboveMaxElev + res.countBeyondHorizon;
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    std::cout << "[RADAR] Compute       : " << ms << " ms\n";
    std::cout << "  COVERED        : " << res.countCovered        << "\n";
    std::cout << "  TERRAIN_MASKED : " << res.countTerrainMasked  << "\n";
    std::cout << "  BELOW_MIN_ELEV : " << res.countBelowMinElev   << "\n";
    std::cout << "  ABOVE_MAX_ELEV : " << res.countAboveMaxElev   << "\n";
    std::cout << "  BEYOND_HORIZON : " << res.countBeyondHorizon  << "\n";
    std::cout << "  OUTSIDE_SECTOR : " << res.countOutsideSector  << "\n";
    if (total > 0)
        std::cout << "  Coverage %     : " << (100.0 * res.countCovered / total) << "%\n";
    return res;
}

/******************************************************************************
 * RADAR ENGINE — computeFast (unchanged)
 ******************************************************************************/
std::vector<uint8_t> RadarEngine::computeFast(int sectorAzimuthDeg)
{
    int sectorAz = (sectorAzimuthDeg < 0) ? -1 : snapTo45(sectorAzimuthDeg);
    auto activeBoundary = sectorBoundary(sectorAz);

    auto t0 = std::chrono::high_resolution_clock::now();

    for (int r = rowMin; r < rowMax; r++) {
        std::memset(fastBuf.data() + (size_t)r * gridCols + colMin, 0, colMax - colMin);
    }
    fastBuf[(size_t)srcRow * gridCols + srcCol] = 1;

    const size_t N      = activeBoundary.size();
    const short* demPtr = g.bin.elev_matrix.data();
    uint8_t*  statPtr   = fastBuf.data();
    const int COLS = gridCols, ROWS = gridRows;
    const double SE = srcElev, TH = targetHeight;
    const double MNS = minSlopeRatio, MXS = maxSlopeRatio;
    const double HOR = horizonM;
    const int    MRC = maxRadiusCells;

    #pragma omp parallel for schedule(dynamic, 64) default(none) \
        shared(demPtr, statPtr, activeBoundary) \
        firstprivate(N, COLS, ROWS, SE, TH, MNS, MXS, HOR, MRC)
    for (size_t i = 0; i < N; ++i) {
        int er = activeBoundary[i].x, ec = activeBoundary[i].y;
        if (er >= 0 && ec >= 0 && er < ROWS && ec < COLS)
            castRadarRayFast(srcRow, srcCol, er, ec,
                             SE, TH, MNS, MXS, HOR, MRC,
                             demPtr, statPtr, COLS, ROWS);
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    int vis = 0;
    for (int r = rowMin; r < rowMax; r++)
        for (int c = colMin; c < colMax; c++)
            vis += fastBuf[(size_t)r * gridCols + c];

    std::cout << "[RADAR] Fast ("
              << (sectorAz < 0 ? "360deg" : std::to_string(sectorAz) + "deg")
              << ")  : "
              << std::chrono::duration<double, std::milli>(t1 - t0).count()
              << " ms  covered=" << vis << "\n";
    return fastBuf;
}

/******************************************************************************
 * SVG RENDERING — shared implementation used by writeToSVG and exportResultToSVG
 *
 * Renders the bounding box around the radar source, downsampled to at most
 * MAX_DIM x MAX_DIM pixels. Adjacent cells with the same status are merged
 * horizontally into single <rect> elements to keep file sizes manageable.
 ******************************************************************************/
static void renderCoverageSVGImpl(
    const uint8_t* statusData,
    int gridRows,  int gridCols,
    int srcRow,    int srcCol,
    int radiusCells,
    int sectorAz,
    const std::string& filename)
{
    // Bounding box (same as engine's rowMin/rowMax/colMin/colMax)
    int rMin = std::max(0,        srcRow - radiusCells);
    int rMax = std::min(gridRows, srcRow + radiusCells + 1);
    int cMin = std::max(0,        srcCol - radiusCells);
    int cMax = std::min(gridCols, srcCol + radiusCells + 1);

    int bbH = rMax - rMin;
    int bbW = cMax - cMin;

    // Downsample so the longest edge is at most MAX_DIM pixels
    const int MAX_DIM = 800;
    int step  = std::max(1, (int)std::ceil((double)std::max(bbH, bbW) / MAX_DIM));
    int svgH  = (bbH + step - 1) / step;
    int svgW  = (bbW + step - 1) / step;

    // Colour table indexed by CoverageStatus uint8 value
    static const char* COLOR_TABLE[] = {
        "#0d1117",   // 0  OUTSIDE_RANGE  — background (skipped)
        "#2ecc71",   // 1  COVERED        — green
        "#e74c3c",   // 2  TERRAIN_MASKED — red
        "#f1c40f",   // 3  BELOW_MIN_ELEV — yellow
        "#e67e22",   // 4  ABOVE_MAX_ELEV — orange
        "#3498db",   // 5  BEYOND_HORIZON — blue
        "#7f8c8d",   // 6  OUTSIDE_SECTOR — grey
    };
    const int N_COLORS = (int)(sizeof(COLOR_TABLE) / sizeof(COLOR_TABLE[0]));

    std::ofstream svg(filename);
    if (!svg) {
        std::cerr << "[SVG] Cannot open: " << filename << "\n";
        return;
    }

    // ---- SVG header ----
    svg << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        << "<svg xmlns=\"http://www.w3.org/2000/svg\""
        << " width=\"" << svgW << "\" height=\"" << svgH << "\""
        << " viewBox=\"0 0 " << svgW << " " << svgH << "\">\n"
        << "<rect width=\"100%\" height=\"100%\" fill=\"#0d1117\"/>\n";

    // ---- Coverage pixels (horizontal run-length merged) ----
    for (int sr = 0; sr < svgH; ++sr) {
        int gridR = rMin + sr * step;
        if (gridR >= gridRows) continue;

        int sc = 0;
        while (sc < svgW) {
            int gridC = cMin + sc * step;
            if (gridC >= gridCols) break;

            uint8_t s = statusData[(size_t)gridR * gridCols + gridC];
            if (s == 0) { ++sc; continue; }    // skip background

            // Extend run while same status
            int runEnd = sc + 1;
            while (runEnd < svgW) {
                int gC2 = cMin + runEnd * step;
                if (gC2 >= gridCols)                              break;
                if (statusData[(size_t)gridR * gridCols + gC2] != s) break;
                ++runEnd;
            }

            const char* col = (s < N_COLORS) ? COLOR_TABLE[s] : "#ffffff";
            svg << "<rect x=\""      << sc
                << "\" y=\""         << sr
                << "\" width=\""     << (runEnd - sc)
                << "\" height=\"1\" fill=\"" << col << "\"/>\n";
            sc = runEnd;
        }
    }

    // ---- Source marker (crosshair + circle) ----
    int msx = (srcCol - cMin) / step;
    int msy = (srcRow - rMin) / step;
    svg << "<line x1=\"" << (msx - 7) << "\" y1=\"" << msy
        << "\" x2=\"" << (msx + 7)   << "\" y2=\"" << msy
        << "\" stroke=\"white\" stroke-width=\"1\" opacity=\"0.9\"/>\n"
        << "<line x1=\"" << msx << "\" y1=\"" << (msy - 7)
        << "\" x2=\"" << msx   << "\" y2=\"" << (msy + 7)
        << "\" stroke=\"white\" stroke-width=\"1\" opacity=\"0.9\"/>\n"
        << "<circle cx=\"" << msx << "\" cy=\"" << msy
        << "\" r=\"4\" fill=\"none\" stroke=\"white\" stroke-width=\"1.5\"/>\n";

    // ---- Legend ----
    struct LI { const char* color; const char* label; };
    static const LI LEGEND[] = {
        { "#2ecc71", "Covered"        },
        { "#e74c3c", "Terrain Masked" },
        { "#f1c40f", "Below Min Elev" },
        { "#e67e22", "Above Max Elev" },
        { "#3498db", "Beyond Horizon" },
        { "#7f8c8d", "Outside Sector" },
    };
    const int NL = (int)(sizeof(LEGEND) / sizeof(LEGEND[0]));
    const int LX = 6, LY = 6, LW = 132, LH = NL * 16 + 10;

    svg << "<rect x=\""  << LX << "\" y=\"" << LY
        << "\" width=\"" << LW << "\" height=\"" << LH
        << "\" fill=\"#000000\" fill-opacity=\"0.65\" rx=\"4\"/>\n";

    for (int i = 0; i < NL; ++i) {
        int iy = LY + 6 + i * 16;
        svg << "<rect x=\""   << (LX + 6) << "\" y=\"" << iy
            << "\" width=\"10\" height=\"10\" fill=\"" << LEGEND[i].color << "\"/>\n"
            << "<text x=\""  << (LX + 20) << "\" y=\"" << (iy + 9)
            << "\" font-size=\"10\" fill=\"#e0e0e0\""
            << " font-family=\"monospace\">" << LEGEND[i].label << "</text>\n";
    }

    // ---- Info line at bottom ----
    std::string secStr = (sectorAz < 0) ? "360" : (std::to_string(sectorAz) + " +/-22.5");
    svg << "<text x=\""  << (svgW / 2)
        << "\" y=\""     << (svgH - 5)
        << "\" font-size=\"10\" fill=\"#888888\" font-family=\"monospace\""
        << " text-anchor=\"middle\">"
        << "src=(" << srcRow << "," << srcCol << ")"
        << "  rad=" << radiusCells
        << "  sector=" << secStr << "deg"
        << "  step=" << step
        << "</text>\n"
        << "</svg>\n";

    std::cout << "[SVG]  Written: " << filename
              << "  (" << svgW << " x " << svgH << " px, step=" << step << ")\n";
}

/******************************************************************************
 * RADAR ENGINE — writeToSVG  (was a stub — now delegates to renderCoverageSVGImpl)
 ******************************************************************************/
void RadarEngine::writeToSVG(const RadarCoverageResult& res, const std::string& filename)
{
    // CoverageStatus is : uint8_t so the cast is safe
    renderCoverageSVGImpl(
        reinterpret_cast<const uint8_t*>(res.status.data()),
        res.gridRows, res.gridCols,
        res.srcRow,   res.srcCol,
        res.radiusCells, res.sectorAzimuthDeg,
        filename);
}

/******************************************************************************
 * BATCH COMPUTE AND SAVE
 *
 * Binary layout:
 *   BinFileHeader                             (20 bytes)
 *   For each cell:
 *     BinResultHeader                         (16 bytes)
 *     uint8_t status[gridRows * gridCols]     (gridRows * gridCols bytes)
 ******************************************************************************/
void batchComputeAndSave(
    Global&                          g,
    const std::vector<FeasibleCell>& cells,
    const BatchParams&               params,
    const std::string&               outFile)
{
    std::ofstream ofs(outFile, std::ios::binary);
    if (!ofs) {
        std::cerr << "[BATCH] Cannot open output file: " << outFile << "\n";
        return;
    }

    const int32_t gRows = g.conf.grid.rows;
    const int32_t gCols = g.conf.grid.cols;
    const size_t  bufBytes = (size_t)gRows * gCols;   // 1 byte per cell
    const size_t  nCells   = cells.size();

    // Write file header
    BinFileHeader fh;
    std::memcpy(fh.magic, "RCOV", 4);
    fh.version    = 1;
    fh.numResults = static_cast<uint32_t>(nCells);
    fh.gridRows   = gRows;
    fh.gridCols   = gCols;
    ofs.write(reinterpret_cast<const char*>(&fh), sizeof(fh));

    double totalMB = (sizeof(BinFileHeader) +
                      nCells * (sizeof(BinResultHeader) + bufBytes)) / 1048576.0;
    std::cout << "[BATCH] Start: " << nCells << " cells, "
              << "estimated output ~" << totalMB << " MB -> " << outFile << "\n";

    auto wallT0 = std::chrono::high_resolution_clock::now();

    for (size_t i = 0; i < nCells; ++i) {
        const auto& cell = cells[i];
        std::cout << "\n[BATCH] ---- Cell " << (i + 1) << "/" << nCells
                  << "  r=" << cell.r << " c=" << cell.c << " ----\n";

        // Construct engine with direct grid coords
        RadarEngine engine(
            g,
            cell.r, cell.c,
            params.antennaHeightM,
            params.targetHeightM,
            params.minElevDeg,
            params.maxElevDeg,
            params.maxRangeKM);

        auto res = engine.compute(params.sectorAzimuthDeg);

        // Per-result header
        BinResultHeader rh;
        rh.srcRow           = static_cast<int32_t>(res.srcRow);
        rh.srcCol           = static_cast<int32_t>(res.srcCol);
        rh.radiusCells      = static_cast<int32_t>(res.radiusCells);
        rh.sectorAzimuthDeg = static_cast<int32_t>(res.sectorAzimuthDeg);
        ofs.write(reinterpret_cast<const char*>(&rh), sizeof(rh));

        // Status buffer: CoverageStatus is uint8_t, so bufBytes == nElements
        ofs.write(
            reinterpret_cast<const char*>(res.status.data()),
            static_cast<std::streamsize>(bufBytes));
    }

    auto wallT1 = std::chrono::high_resolution_clock::now();
    std::cout << "\n[BATCH] Finished: " << nCells << " results in "
              << std::chrono::duration<double>(wallT1 - wallT0).count() << " s  ("
              << totalMB << " MB)\n";
}

/******************************************************************************
 * LOAD COVERAGE RESULTS FROM BINARY FILE
 ******************************************************************************/
std::vector<StoredCoverageResult> loadCoverageResults(const std::string& binFile)
{
    std::ifstream ifs(binFile, std::ios::binary);
    if (!ifs) {
        std::cerr << "[LOAD] Cannot open: " << binFile << "\n";
        return {};
    }

    BinFileHeader fh;
    ifs.read(reinterpret_cast<char*>(&fh), sizeof(fh));
    if (!ifs || std::strncmp(fh.magic, "RCOV", 4) != 0) {
        std::cerr << "[LOAD] Bad magic — not an RCOV file\n";
        return {};
    }

    std::cout << "[LOAD] File: " << binFile
              << "  v=" << fh.version
              << "  results=" << fh.numResults
              << "  grid=" << fh.gridRows << "x" << fh.gridCols << "\n";

    const size_t bufBytes = (size_t)fh.gridRows * fh.gridCols;
    std::vector<StoredCoverageResult> results;
    results.reserve(fh.numResults);

    for (uint32_t i = 0; i < fh.numResults; ++i) {
        BinResultHeader rh;
        ifs.read(reinterpret_cast<char*>(&rh), sizeof(rh));
        if (!ifs) {
            std::cerr << "[LOAD] Unexpected EOF at result header " << i << "\n";
            break;
        }

        StoredCoverageResult r;
        r.srcRow           = rh.srcRow;
        r.srcCol           = rh.srcCol;
        r.radiusCells      = rh.radiusCells;
        r.sectorAzimuthDeg = rh.sectorAzimuthDeg;
        r.gridRows         = fh.gridRows;
        r.gridCols         = fh.gridCols;
        r.status.resize(bufBytes);

        ifs.read(reinterpret_cast<char*>(r.status.data()),
                 static_cast<std::streamsize>(bufBytes));
        if (!ifs) {
            std::cerr << "[LOAD] Unexpected EOF reading status buffer " << i << "\n";
            break;
        }

        std::cout << "[LOAD]   [" << (i + 1) << "] src=("
                  << r.srcRow << "," << r.srcCol << ")"
                  << "  radius=" << r.radiusCells
                  << "  sector=" << r.sectorAzimuthDeg << "\n";

        results.push_back(std::move(r));
    }

    std::cout << "[LOAD] Loaded " << results.size() << " results\n";
    return results;
}

/******************************************************************************
 * EXPORT A LOADED RESULT TO SVG  (verification helper)
 ******************************************************************************/
void exportResultToSVG(const StoredCoverageResult& result, const std::string& filename)
{
    renderCoverageSVGImpl(
        result.status.data(),
        result.gridRows, result.gridCols,
        result.srcRow,   result.srcCol,
        result.radiusCells, result.sectorAzimuthDeg,
        filename);
}

/******************************************************************************
 * MAIN
 ******************************************************************************/
int main()
{
    auto tMain0 = std::chrono::high_resolution_clock::now();
    Global g;

    // ------------------------------------------------------------------
    // Original single-engine demo (unchanged)
    // ------------------------------------------------------------------
    std::cout << "--- Initializing Radar Coverage Engine (lat/lon) ---\n";
    RadarEngine engine(
        g,
        33.9871, 74.7742,   // Srinagar airport
        50.0, 1000.0,
        10.0, 30.0, 100.0);

    std::cout << "\n--- Computing Full 360 Coverage ---\n";
    auto res360 = engine.compute();
    engine.writeToSVG(res360, "radar_360.svg");

    std::cout << "\n--- Fast Mode ---\n";
    engine.computeFast();
    engine.computeFast(135);

    // ------------------------------------------------------------------
    // NEW: Batch pipeline — FeasibleCell grid coords, no lat/lon needed
    // ------------------------------------------------------------------
    std::cout << "\n========== BATCH PIPELINE ==========\n";
    auto f = Bin::getFeasibleGrid(g.bin);

    BatchParams params;
    params.antennaHeightM   = 50.0;
    params.targetHeightM    = 1000.0;
    params.minElevDeg       = 10.0;
    params.maxElevDeg       = 30.0;
    params.maxRangeKM       = 100.0;
    params.sectorAzimuthDeg = -1;      // full 360 for all cells

    const std::string binFile = "reserved.bin";

    // --- Step 1: compute and persist ---
    batchComputeAndSave(g, f.cells, params, binFile);

    // --- Step 2: load back ---
    std::cout << "\n--- Loading Results from Disk ---\n";
    auto loaded = loadCoverageResults(binFile);

    // --- Step 3: export first result to SVG for visual verification ---
    if (!loaded.empty()) {
        std::cout << "\n--- Exporting First Result to SVG (verification) ---\n";
        exportResultToSVG(loaded[0], "coverage_verify.svg");

        // Optionally also export all results
        for (size_t i = 1; i < loaded.size(); ++i) {
            std::string svgName = "coverage_result_" + std::to_string(i + 1) + ".svg";
            exportResultToSVG(loaded[i], svgName);
        }
    } else {
        std::cerr << "[MAIN] No results loaded — check binary file\n";
    }

    auto tMain1 = std::chrono::high_resolution_clock::now();
    std::cout << "\n==============================================\n"
              << ">> TOTAL PIPELINE RUN TIME: "
              << std::chrono::duration<double, std::milli>(tMain1 - tMain0).count() / 1000.0
              << " s\n==============================================\n";

    return 0;
}
