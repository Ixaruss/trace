#ifndef OPERATIONS_H_
#define OPERATIONS_H_

#include "common.h"
#include "common/config.h"
#include "preprocess/Config.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

constexpr double EARTH_R = 6371000.0;
constexpr double INV_2EARTHR = 1.0 / (2.0 * EARTH_R);
constexpr double K_ATMOS = 4.0 / 3.0;
constexpr double K_EFF_R = K_ATMOS * EARTH_R;
constexpr double INV_2_KEFFR = 1.0 / (2.0 * K_EFF_R);
constexpr double DEG2RAD = 3.14159265358979 / 180.0;
constexpr double RAD2DEG = 180.0 / 3.14159265358979;

enum class CoverageStatus : uint8_t {
  OUTSIDE_RANGE = 0,
  COVERED = 1,
  TERRAIN_MASKED = 2,
  BELOW_MIN_ELEV = 3,
  ABOVE_MAX_ELEV = 4,
  BEYOND_HORIZON = 5,
  OUTSIDE_SECTOR = 6, // NEW: cell not in the scanned sector
};

enum class LOSStatus : uint8_t {
  VISIBLE = 0,
  BEYOND_RANGE = 1,
  BEYOND_HORIZON = 2,
  BELOW_MIN_ELEV = 3,
  ABOVE_MAX_ELEV = 4,
  TERRAIN_MASKED = 5,
  NO_DATA = 6, // DEM has no data along ray
  OUT_OF_GRID = 7,
};

enum class OutputCrs { Native7755, Wgs84, Both };

// Band layout of the emitted GeoTIFF.
//   Rgba              -- 4-band R,G,B,A Byte raster (default). True alpha,
//                        honoured by everything, ~4x the pixel data.
//   SingleBandPalette -- 1-band Byte raster + color table. ~1/4 the size,
//                        but transparency is expressed as NoData (index 0),
//                        which only GIS-aware readers respect.
enum class BandLayout { Rgba, SingleBandPalette };

static const char *losStatusStr(LOSStatus s) {
  switch (s) {
  case LOSStatus::VISIBLE:
    return "VISIBLE";
  case LOSStatus::BEYOND_RANGE:
    return "BEYOND_RANGE";
  case LOSStatus::BEYOND_HORIZON:
    return "BEYOND_HORIZON";
  case LOSStatus::BELOW_MIN_ELEV:
    return "BELOW_MIN_ELEV";
  case LOSStatus::ABOVE_MAX_ELEV:
    return "ABOVE_MAX_ELEV";
  case LOSStatus::TERRAIN_MASKED:
    return "TERRAIN_MASKED";
  case LOSStatus::NO_DATA:
    return "NO_DATA";
  case LOSStatus::OUT_OF_GRID:
    return "OUT_OF_GRID";
  }
  return "UNKNOWN";
}

/******************************************************************************
 * FAST RESULT — single struct, no allocation
 ******************************************************************************/

/******************************************************************************
 * DETAILED CELL — per-cell breakdown for visualization
 ******************************************************************************/

struct LOSCell {
  int row, col;
  double distM;
  float terrainElev;    // raw terrain elevation (m)
  float correctedElev;  // terrain - curvature drop (m)
  double terrainSlope;  // (correctedElev - srcElev) / dist
  double maxSlopeSoFar; // running max terrain slope before this cell
  double targetSlope;   // (target_abs_elev - srcElev) / targetDist (fixed)
  bool blocks;          // terrainSlope > targetSlope at this cell
};

/******************************************************************************
 * INTERNAL: precompute all ray parameters that don't depend on DEM
 ******************************************************************************/

struct LOSParams {
  // source
  int srcRow, srcCol;
  double srcElev; // absolute source elevation (m)

  // target
  int tgtRow, tgtCol;
  double tgtAbsElev;   // target terrain + targetHeight (m)
  double tgtDist;      // ground distance to target (m)
  double targetSlope;  // (tgtAbsElev - srcElev) / tgtDist — FIXED for ray
  double elevAngleDeg; // elevation angle to target

  // limits
  double maxRangeM;
  double horizonM;
  double minSlope; // tan(minElevDeg)
  double maxSlope; // tan(maxElevDeg)

  // grid
  int COLS, ROWS;
};
struct RadarLOSResult {
  LOSStatus status;
  bool visible;

  double distM;        // ground distance to target (m)
  double elevAngleDeg; // elevation angle to target (degrees)
  double horizonM;     // radar+target LOS horizon distance (m)

  // if TERRAIN_MASKED: where the ray was first blocked
  int blockRow;
  int blockCol;
  double blockDistM;
  float blockElevM; // terrain elevation at block point

  // timings
  double computeUs; // microseconds
};

struct RadarCoverageResult {
  std::vector<CoverageStatus> status;
  int srcRow, srcCol;     // window-local — same space as gridCols/gridRows
                          // below and as status[]. Add winR0/winC0 to get
                          // back to the full grid's coordinates.
  int radiusCells;
  int gridCols, gridRows; // window dims (== full-grid dims under FULL)
  int winR0 = 0, winC0 = 0; // global row/col of this window's (0,0) —
                            // needed by render.cc's writeCoverageToGeoTIFF
                            // to place the output at the right geographic
                            // origin; status[]/srcRow/srcCol themselves
                            // stay window-local, only the georeferencing
                            // math needs this
  double sectorLeftDeg;  // arc start (clockwise from North). < 0 = full 360
  double sectorRightDeg; // arc end   (clockwise from North). < 0 = full 360
  int countCovered, countTerrainMasked;
  int countBelowMinElev, countAboveMaxElev;
  int countBeyondHorizon, countOutsideRange;
  int countOutsideSector;
};

struct RadarLOSDetailedResult {
  RadarLOSResult summary;
  std::vector<LOSCell> cells;
};

struct Point2D {
  int x, y;
};

/**
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
class ViewshedEngine {
private:
  Global &g;
  int gridCols;
  int gridRows;
  int centerX;
  int centerY;
  int winR0 = 0, winC0 = 0;   // global origin of the resident window;
                              // (centerX,centerY), viewshedMask, and every
                              // row/col this class touches are LOCAL to it
  int maxRadiusCells;
  int rowMin, rowMax, colMin, colMax;
  float radarHeightMeters;
  bool useEarthCurvature;
  bool is_full;
  std::vector<uint8_t> viewshedMask;
  std::vector<Point2D> circlePerimeter;

  void precomputeCirclePerimeter();

  // Reset only the active bounding box (rowMin..rowMax, colMin..colMax)
  // of viewshedMask to 0, instead of std::fill-ing the whole grid.
  // Same optimization as RadarEngine::prefill() / the bounding-box
  // reset in RadarEngine::computeFast() — keeps repeated calls cheap
  // even when maxRadiusCells is small relative to the full grid.
  void prefill();

  // Filter circlePerimeter down to the rays whose azimuth (from source)
  // falls on the [leftDeg, rightDeg] arc, swept clockwise, North = 0°.
  // leftDeg < 0 || rightDeg < 0 -> full 360 (returns circlePerimeter).
  // Reuses RadarEngine's azimuth helpers since both engines share the
  // same row/col grid convention.
  std::vector<Point2D> sectorPerimeter(double leftDeg, double rightDeg) const;

public:
  ViewshedEngine(Global &globalInstance, double srcLat, double srcLon,
                 float rangeKM, float txHeightMeters,
                 bool earthCurvature = true)
      : g(globalInstance), radarHeightMeters(txHeightMeters),
        useEarthCurvature(earthCurvature) {
    auto t0 = std::chrono::high_resolution_clock::now();

    // CHANGED: dims come from the RESIDENT WINDOW, not the full grid.
    // g.bin.elev_matrix only holds g.bin.op_window under AOI/UNBOUNDED —
    // indexing it with the full grid's stride (conf.grid.cols) reads past
    // the end of a much smaller vector. Under FULL, op_window IS the full
    // grid, so this is a no-op change for existing callers.
    if (!g.bin.window_valid) {
      std::cerr << "[ERROR] ViewshedEngine: Bin has no resolved window — "
                   "construct Global with an operation type first\n";
    }
    // Grid geometry from data_window — the layout elev_matrix is ACTUALLY
    // in. Tiled: the loaded window. load_full: the whole grid (r0=c0=0).
    // Both the stride (gridCols) and the origin (winR0/winC0) must come from
    // the SAME window. They used to come from different ones in load_full
    // mode — full-grid stride, but op_window's sub-window origin — so
    // render.cc added a window offset to coordinates that were already
    // global, and the GeoTIFF landed displaced by exactly that offset.
    const GridWindow& dw = g.bin.data_window;
    gridCols = dw.cols;
    gridRows = dw.rows;
    winR0    = dw.r0;
    winC0    = dw.c0;
    is_full  = g.bin.load_full;

    int globalCenterX = 0, globalCenterY = 0;
    if (!g.bin.getGridCoords(srcLat, srcLon, globalCenterX, globalCenterY)) {
      std::cerr << "[ERROR] Source outside grid\n";
      centerX = gridRows / 2;
      centerY = gridCols / 2;
    } else if (!dw.contains(globalCenterX, globalCenterY)) {
      std::cerr << "[ERROR] Source (" << srcLat << ", " << srcLon
                << ") is not resident — cell (" << globalCenterX << ", "
                << globalCenterY << ") lies outside the loaded data\n";
      centerX = gridRows / 2;
      centerY = gridCols / 2;
    } else {
      // global -> resident-local. In load_full this is the identity
      // (r0=c0=0), so one path covers both layouts.
      centerX = dw.toLocalRow(globalCenterX);
      centerY = dw.toLocalCol(globalCenterY);
    }
    maxRadiusCells =
        static_cast<int>(std::round(rangeKM * 1000.0 / (double)g.bin.cell_m));

    // FIX: centerX is a ROW and centerY is a COLUMN — the row bound must
    // compare against gridROWS and the column bound against gridCOLS, not
    // the other way around as this previously read. Harmless on a
    // near-square grid, wrong on a rectangular one.
    int maxPossible = std::min(std::min(centerX, gridRows - centerX - 1),
                               std::min(centerY, gridCols - centerY - 1));
    if (maxRadiusCells > maxPossible) {
      std::cout << "[WARN] Radius clamped from " << maxRadiusCells << " to "
                << maxPossible << "\n";
      maxRadiusCells = maxPossible;
    }

    // Bounding box around the source, same convention as
    // RadarEngine::rowMin/rowMax/colMin/colMax (centerX = row,
    // centerY = col here). All window-local now.
    rowMin = std::max(0, centerX - maxRadiusCells - 1);
    rowMax = std::min(gridRows, centerX + maxRadiusCells + 2);
    colMin = std::max(0, centerY - maxRadiusCells - 1);
    colMax = std::min(gridCols, centerY + maxRadiusCells + 2);

    // CHANGED: sized to the window, not the full grid. This is the actual
    // memory win — e.g. ~1GB -> ~11MB for a 50km UNBOUNDED window on the
    // 100m grid.
    viewshedMask.assign((size_t)gridCols * gridRows, 0);
    precomputeCirclePerimeter();

    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "[METRIC] Init                         : "
              << std::chrono::duration<double, std::milli>(t1 - t0).count()
              << " ms\n";
    std::cout << "[INFO]   Center  (" << centerX << ", " << centerY
              << ")  [local; global = (" << (centerX + winR0) << ", "
              << (centerY + winC0) << ")]\n";
    std::cout << "[INFO]   Radius  " << maxRadiusCells
              << " cells = " << (maxRadiusCells * (double)g.bin.cell_m / 1000.0)
              << " km @ " << g.bin.cell_m << " m/cell\n";
    std::cout << "[INFO]   Rays    " << circlePerimeter.size() << "\n";
    std::cout << "[INFO]   Curv    " << (useEarthCurvature ? "ON" : "OFF")
              << "\n";
  }

  // Window-local (row, col) -> global grid (row, col). Needed by any
  // downstream writer (GeoTIFF, overlap merge) that must place cells back
  // in world coordinates rather than window-local ones.
  inline int toGlobalRow(int localRow) const { return localRow + winR0; }
  inline int toGlobalCol(int localCol) const { return localCol + winC0; }

  // leftDeg/rightDeg: arc swept clockwise from leftDeg to rightDeg,
  //                    North = 0°. leftDeg < 0 || rightDeg < 0 -> full 360.
  //                    leftDeg > rightDeg wraps through North.
  void computeViewshed(double leftDeg = -1.0, double rightDeg = -1.0);

  inline const std::vector<uint8_t> &getMask() const { return viewshedMask; }

  void writeViewshedToSVG(const std::string &filename);
  void writeViewshedToPNG(const std::string &filename);
  static std::string base64Encode(const uint8_t *data, size_t length);
  void writeViewshedToGeoTIFF(const std::string &filename,
                              int reductionPercent = 0,
                              OutputCrs outCrs = OutputCrs::Native7755,
                              BandLayout bands = BandLayout::Rgba);
};

class RadarEngine {
public:
  Global &g;
  int gridCols, gridRows;
  int srcRow, srcCol;
  int winR0 = 0, winC0 = 0;   // global origin of the resident window; every
                              // row/col in this class (srcRow/srcCol,
                              // statusBuf, fastBuf, boundary) is LOCAL to it
  int maxRadiusCells;
  int rowMin, rowMax, colMin, colMax;
  double leftDeg = -1.0;
  double rightDeg = -1.0;
  bool use_radar_model;
  int peak_power;
  int antenna_gain;
  double wavelength;
  int crossSection;
  double detachable_threshold;
  bool is_full;
  double antennaHeight, targetHeight;
  double minSlopeRatio, maxSlopeRatio;
  double maxRangeM, horizonM, srcElev;
  bool use_curvature;
  bool use_amsl;
  bool use_sector;
  std::vector<uint8_t> statusBuf;
  std::vector<uint8_t> fastBuf;
  std::vector<Point2D> boundary; // full 360 boundary

  // ------------------------------------------------------------------
  // Filter boundary cells to those on the [leftDeg, rightDeg] arc,
  // swept clockwise from leftDeg to rightDeg (North = 0°).
  // leftDeg < 0 (or rightDeg < 0) → full 360°.
  // ------------------------------------------------------------------
  std::vector<Point2D> sectorBoundary(double leftDeg, double rightDeg) const;
  void prefill(double leftDeg, double rightDeg, uint8_t *buf,
               uint8_t fillInside, uint8_t fillOutside);

  /******************************************************************************
   * AZIMUTH HELPER
   *
   * North-based clockwise azimuth from source to (row,col).
   * North = 0°, East = 90°, South = 180°, West = 270°
   ******************************************************************************/

  static inline double cellAzimuth(int srcRow, int srcCol, int r, int c) {
    double drow = r - srcRow;
    double dcol = c - srcCol;
    // atan2(east, north) = clockwise from north
    double az = std::atan2(dcol, -drow) * RAD2DEG;
    if (az < 0.0)
      az += 360.0;
    return az;
  }

  /******************************************************************************
   * ARC HELPER
   *
   * Returns true if azimuth `az` lies on the arc swept clockwise from
   * `leftDeg` to `rightDeg` (both inclusive). North = 0°, clockwise positive.
   *
   * Handles the 0°/360° wraparound: if leftDeg > rightDeg the arc is taken
   * to cross North (e.g. left=350°, right=20° covers 350°..360°,0°..20°).
   * All inputs are normalised into [0,360) first, so callers may pass any
   * value (including negative angles or values > 360).
   ******************************************************************************/

  static inline double normAz(double deg) {
    double a = std::fmod(deg, 360.0);
    if (a < 0.0)
      a += 360.0;
    return a;
  }

  static inline bool inArc(double az, double leftDeg, double rightDeg) {
    az = normAz(az);
    leftDeg = normAz(leftDeg);
    rightDeg = normAz(rightDeg);

    if (leftDeg <= rightDeg)
      return az >= leftDeg && az <= rightDeg;
    // arc wraps through 0°/360°
    return az >= leftDeg || az <= rightDeg;
  }

public:
  RadarEngine(Global &globalInstance) : g(globalInstance) {
    if (!Rconf::is_available()) {
      std::cerr << "[ERROR] radar conf not found\n";
      return;
    }
    auto rconf = Rconf::load();

    use_radar_model = rconf.use_radar_model;
    peak_power = rconf.peak_power;
    antenna_gain = rconf.antenna_gain;
    wavelength = rconf.wavelength;
    crossSection = rconf.crossSection;
    detachable_threshold = rconf.detachable_threshold;
    use_sector = rconf.use_sector;
    antennaHeight = rconf.antennaHeight;
    targetHeight = rconf.target;
    maxRangeM = rconf.range * 1000, use_curvature = rconf.use_curv;
    use_amsl = rconf.useAMSL;
    leftDeg = rconf.loa;
    rightDeg = rconf.roa;
    auto t0 = std::chrono::high_resolution_clock::now();

    // CHANGED: window dims, not full-grid — see ViewshedEngine for why.
    if (!g.bin.window_valid) {
      std::cerr << "[ERROR] RadarEngine: Bin has no resolved window — "
                   "construct Global with an operation type first\n";
    }
    // Same resident-layout rule as ViewshedEngine: stride AND origin from
    // data_window, so render.cc's georeferencing is right in both layouts.
    const GridWindow& dw = g.bin.data_window;
    gridCols = dw.cols;
    gridRows = dw.rows;
    winR0    = dw.r0;
    winC0    = dw.c0;
    is_full  = g.bin.load_full;

    int globalSrcRow = 0, globalSrcCol = 0;
    if (!g.bin.getGridCoords(rconf.srcLat, rconf.srcLon, globalSrcRow, globalSrcCol)) {
      std::cerr << "[ERROR] Source outside grid\n";
      srcRow = gridRows / 2;
      srcCol = gridCols / 2;
    } else if (!dw.contains(globalSrcRow, globalSrcCol)) {
      std::cerr << "[ERROR] Source (" << rconf.srcLat << ", " << rconf.srcLon
                << ") is not resident — cell (" << globalSrcRow << ", "
                << globalSrcCol << ") lies outside the loaded data\n";
      srcRow = gridRows / 2;
      srcCol = gridCols / 2;
    } else {
      srcRow = dw.toLocalRow(globalSrcRow);
      srcCol = dw.toLocalCol(globalSrcCol);
    }

    maxRadiusCells = (int)std::ceil(maxRangeM / (double)g.bin.cell_m);
    int maxPossible = std::min(std::min(srcRow, gridRows - srcRow - 1),
                               std::min(srcCol, gridCols - srcCol - 1));
    if (maxRadiusCells > maxPossible) {
      std::cout << "[WARN] Radius clamped to " << maxPossible << "\n";
      maxRadiusCells = maxPossible;
    }

    rowMin = std::max(0, srcRow - maxRadiusCells - 1);
    rowMax = std::min(gridRows, srcRow + maxRadiusCells + 2);
    colMin = std::max(0, srcCol - maxRadiusCells - 1);
    colMax = std::min(gridCols, srcCol + maxRadiusCells + 2);

    // CHANGED: sized to the window, not the full grid.
    statusBuf.assign((size_t)gridRows * gridCols, 0);
    fastBuf.assign((size_t)gridRows * gridCols, 0);

    horizonM = std::sqrt(2.0 * K_EFF_R * antennaHeight) +
               std::sqrt(2.0 * K_EFF_R * targetHeight);
    if (horizonM > maxRangeM)
      horizonM = maxRangeM;

    minSlopeRatio = std::tan(rconf.minAngle * DEG2RAD);
    maxSlopeRatio = std::tan(rconf.maxAngle * DEG2RAD);

    // srcRow/srcCol and gridCols are now consistently window-local, so this
    // index is correct against the (possibly windowed) elev_matrix without
    // any other change.
    short srcTerrain = g.bin.elev_matrix[(size_t)srcRow * gridCols + srcCol];
    srcElev = static_cast<double>(srcTerrain) + antennaHeight;

    boundary = circleBoundary(srcRow, srcCol, maxRadiusCells);

    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "[RADAR] Init          : "
              << std::chrono::duration<double, std::milli>(t1 - t0).count()
              << " ms\n";
    std::cout << "[RADAR]   Src terrain : " << srcTerrain
              << "m  elev: " << srcElev << "m\n";
    std::cout << "[RADAR]   Src (local) : (" << srcRow << ", " << srcCol
              << ")  [global = (" << (srcRow + winR0) << ", "
              << (srcCol + winC0) << ")]\n";
    std::cout << "[RADAR]   Horizon     : " << horizonM / 1000.0 << " km\n";
    std::cout << "[RADAR]   Full rays   : " << boundary.size() << "\n";
  }

  // Window-local (row, col) -> global grid (row, col).
  inline int toGlobalRow(int localRow) const { return localRow + winR0; }
  inline int toGlobalCol(int localCol) const { return localCol + winC0; }

  inline void reload_conf() {
    if (!Rconf::is_available()) {
      std::cerr << "[ERROR] radar conf not found\n";
      return;
    }
    auto rconf = Rconf::load();
    use_radar_model = rconf.use_radar_model;
    peak_power = rconf.peak_power;
    antenna_gain = rconf.antenna_gain;
    wavelength = rconf.wavelength;
    crossSection = rconf.crossSection;
    detachable_threshold = rconf.detachable_threshold;
    use_sector = rconf.use_sector;
    antennaHeight = rconf.antennaHeight;
    targetHeight = rconf.target;
    maxRangeM = rconf.range * 1000, use_curvature = rconf.use_curv;
    use_amsl = rconf.useAMSL;
    leftDeg = rconf.loa;
    rightDeg = rconf.roa;

    // CHANGED: same global->local translation as the constructor. The
    // previous version wrote the GLOBAL row/col into srcRow/srcCol while
    // gridCols/gridRows (and therefore every buffer) stayed window-local —
    // the same out-of-bounds bug as the constructor had, just triggered on
    // a config reload instead of construction.
    int globalSrcRow, globalSrcCol;
    if (!g.bin.getGridCoords(rconf.srcLat, rconf.srcLon, globalSrcRow, globalSrcCol)) {
      std::cerr << "[ERROR] Source outside grid\n";
      srcRow = gridRows / 2;
      srcCol = gridCols / 2;
    } else if (!g.bin.data_window.contains(globalSrcRow, globalSrcCol)) {
      std::cerr << "[ERROR] Source (" << rconf.srcLat << ", " << rconf.srcLon
                << ") is not resident — reload Global with a window that "
                   "covers it before reload_conf()\n";
      srcRow = gridRows / 2;
      srcCol = gridCols / 2;
    } else {
      // data_window, matching the constructor — op_window here would give
      // sub-window-local coords against full-grid buffers in load_full mode
      srcRow = g.bin.data_window.toLocalRow(globalSrcRow);
      srcCol = g.bin.data_window.toLocalCol(globalSrcCol);
    }
    minSlopeRatio = std::tan(rconf.minAngle * DEG2RAD);
    maxSlopeRatio = std::tan(rconf.maxAngle * DEG2RAD);
  }
  // leftDeg/rightDeg: arc swept clockwise from leftDeg to rightDeg, North = 0°.
  // Defaults (-1, -1) mean full 360°. leftDeg > rightDeg wraps through North.
  RadarCoverageResult compute();

  std::vector<uint8_t> computeWithAllPoints();

  void writeToSVG(const RadarCoverageResult &res,
                  const std::string &filename = "radar_coverage.svg");

  static std::vector<Point2D> circleBoundary(int cr, int cc, int radius);

  void writeToPNG(const RadarCoverageResult &res,
                  const std::string &filename = "radar_coverage.png");

public:
  static RadarLOSResult radarLOS(Global &g, double srcLat, double srcLon,
                                 double antennaHeightM, double tgtLat,
                                 double tgtLon, double targetHeightM,
                                 double minElevDeg, double maxElevDeg,
                                 double maxRangeKM);

public:
  static RadarLOSDetailedResult
  radarLOSDetailed(Global &g, double srcLat, double srcLon,
                   double antennaHeightM, double tgtLat, double tgtLon,
                   double targetHeightM, double minElevDeg, double maxElevDeg,
                   double maxRangeKM);

  void writeCoverageToGeoTIFF(const RadarCoverageResult &, const std::string &,
                              int, OutputCrs,
                              BandLayout = BandLayout::SingleBandPalette);
  void writeFastToGeoTIFF(const std::vector<uint8_t> &, const std::string &,
                          int, OutputCrs,
                          BandLayout = BandLayout::SingleBandPalette);

  inline static void printLOSResult(const RadarLOSResult &r, double srcLat,
                                    double srcLon, double tgtLat,
                                    double tgtLon) {
    std::cout << "======================================\n";
    std::cout << "RADAR LOS RESULT\n";
    std::cout << "  Src         : " << srcLat << ", " << srcLon << "\n";
    std::cout << "  Tgt         : " << tgtLat << ", " << tgtLon << "\n";
    std::cout << "  Distance    : " << r.distM / 1000.0 << " km\n";
    std::cout << "  Elev angle  : " << r.elevAngleDeg << "°\n";
    std::cout << "  Horizon     : " << r.horizonM / 1000.0 << " km\n";
    std::cout << "  Status      : " << losStatusStr(r.status) << "\n";
    std::cout << "  Visible     : " << (r.visible ? "YES" : "NO") << "\n";

    if (r.status == LOSStatus::TERRAIN_MASKED) {
      std::cout << "  Block row   : " << r.blockRow << "\n";
      std::cout << "  Block col   : " << r.blockCol << "\n";
      std::cout << "  Block dist  : " << r.blockDistM / 1000.0 << " km\n";
      std::cout << "  Block elev  : " << r.blockElevM << " m\n";
    }
    std::cout << "  Compute     : " << r.computeUs << " µs\n";
    std::cout << "======================================\n";
  }
};

bool writeOverlapGeoTIFF(Global &, const FeasibleGrid &, const std::string &,
                         const std::vector<GridPoint> &, int, OutputCrs,
                         BandLayout = BandLayout::SingleBandPalette);

bool writeOverlapCoreTIFF(
    Global&                        g,
    const uint8_t*                 counts,      // row-major, height*width
    int                            origin_r,
    int                            origin_c,
    int                            height,
    int                            width,
    const std::string&             filename,
    const std::vector<GridPoint>&  placements,
    int                            reductionPercent,
    OutputCrs                      outCrs,
    BandLayout                     bands);
#endif
