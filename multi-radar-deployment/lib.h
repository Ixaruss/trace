#include <cmath>
#include <ogr_core.h>
#include <ogr_geometry.h>
#include <ogr_spatialref.h>
#include <memory>
#include <optional>
#include <vector>
#include "../common.h"
#include "../operations.h"

using GeomPtr = std::unique_ptr<OGRGeometry,
  decltype(&OGRGeometryFactory::destroyGeometry)>;

inline std::optional<std::pair<double,double>> polygonPointOnSurface(const std::vector<std::pair<double,double>>& coords)
{
  if (coords.size() < 3) return std::nullopt;

  OGRLinearRing ring;
  for (const auto& c : coords)
    ring.addPoint(c.first, c.second);
  ring.closeRings();

  OGRPolygon poly;
  poly.addRing(&ring);     // copies the ring

  OGRSpatialReference src, dst;
  if (src.importFromEPSG(4326) != OGRERR_NONE ||
  dst.importFromEPSG(7755) != OGRERR_NONE) {
      std::cout << "[centroid] import failed for epsg.. \n";
  }
  src.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
  dst.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

  std::unique_ptr<OGRCoordinateTransformation,
          decltype(&OCTDestroyCoordinateTransformation)>
    fwd(OGRCreateCoordinateTransformation(&src, &dst),
      OCTDestroyCoordinateTransformation),
    rev(OGRCreateCoordinateTransformation(&dst, &src),
      OCTDestroyCoordinateTransformation);

  if (!fwd || !rev) return std::nullopt;
  if (poly.transform(fwd.get()) != OGRERR_NONE){
      std::cout << "[centroid] polygon transform failed.. \n";
      return std::nullopt;
  }

  OGRPoint pt;
  if(poly.Centroid(&pt)!= OGRERR_NONE)
      return std::nullopt;
  else if(pt.transform(rev.get()) != OGRERR_NONE)
      return std::nullopt;
  // NOTE: under TRADITIONAL_GIS_ORDER this pair is (lon, lat) — i.e. the
  // OPPOSITE convention from conf.area_op's (lat, lon) pairs. Any caller
  // combining a value from here with a value from area_op needs to be
  // deliberate about which slot holds which — see computeBearing below,
  // which currently only works BECAUSE its one call site happens to pass
  // a (lat,lon) point and a (lon,lat) centroid into a signature ordered
  // to match exactly that.
  return std::pair {pt.getX(), pt.getY()};
}


inline bool gridToLatLon(Global &g, int r, int c, double &lat, double &lon) {
  static OGRCoordinateTransformation *inv = nullptr;
  if (!inv) {
    OGRSpatialReference src, dst;
    src.importFromEPSG(7755);
    dst.importFromEPSG(4326);
    src.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    dst.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    inv = OGRCreateCoordinateTransformation(&src, &dst);
    if (!inv)
      return false;
  }
  double x = g.conf.grid.origin_x + (double)c * RESOLUTION;
  double y = g.conf.grid.origin_y - (double)r * RESOLUTION;
  if (!inv->Transform(1, &x, &y))
    return false;
  lon = x;
  lat = y;
  return true;
}

// NOTE ON PARAMETER ORDER: point A is (latA, lonA) — lat first — but point
// B is (lonB, latB) — lon first. This is NOT a typo to "clean up" without
// checking the call site first: it currently only produces correct
// results because get_midpoint_from_two_closest_points_to_sector's one
// call site passes a (lat,lon) polygon point as A and a (lon,lat) centroid
// (see polygonPointOnSurface's return convention above) as B, and this
// signature's mismatched order happens to match that mismatch. If you add
// a second call site, get its argument order right against THIS
// signature, not against "what a bearing function should look like" — or
// better, standardize both call sites on one order and fix this signature
// to match, rather than relying on two wrongs cancelling out.
inline double computeBearing(Global &g,double latA,double lonA,double lonB,double latB) {

    OGRPoint pointA(lonA, latA);
    OGRPoint pointB(lonB, latB);

    // Transform points to the metric projection
    if (pointA.transform(g.bin.poCT) != OGRERR_NONE || pointB.transform(g.bin.poCT) != OGRERR_NONE) {
        std::cerr << "Bearing transformation failed." << std::endl;
        return -1.0;
    }

    double deltaX = pointB.getX() - pointA.getX();
    double deltaY = pointB.getY() - pointA.getY();

    // Calculate angle clockwise from North (Y-axis)
    double bearingDeg = Global::rad2deg(std::atan2(deltaX, deltaY));

    // Normalize to 0° - 360°
    if (bearingDeg < 0) {
        bearingDeg += 360.0;
    }

    return bearingDeg;
}

inline double calc_distance(int x1, int x2, int y1, int y2, int cellSize) {
    double dx = x2 -x1;
    double dy = y2 - y1;
    return std::hypot(dx, dy) * cellSize;
}

// FIX: this previously computed (x1+x2)/2 and (y1+y2)/2 — i.e. it averaged
// a point's own row with its own column (x1,x2 are one point's row/col;
// y1,y2 are the OTHER point's row/col), not the two points' matching
// coordinates. Corrected to average row-with-row and col-with-col.
inline std::pair<int,int> get_midpoint(int row1, int col1, int row2, int col2){
    return {(row1 + row2) / 2, (col1 + col2) / 2};
}

inline std::pair<int,int> get_midpoint_from_two_closest_points_to_sector(Global &g, std::vector<std::pair<double,double>> points, std::pair<double,double> centroid,double loa, double roa){
    std::map<double,std::pair<double,double>> buff1;
    std::map<double,std::pair<double,double>> buff2;
    for (auto pt : points) {
        double ax = computeBearing(g, pt.first, pt.second, centroid.first, centroid.second);
        double d1 = abs(ax - loa);
        double d2 = abs(ax - roa);
        buff1.insert({d1,pt});
        buff2.insert({d2,pt});
    }
    auto it1 = buff1.begin();
    auto it2 = buff2.begin();
    int row1, col1, row2, col2;
    g.bin.getGridCoords(it1->second.first, it1->second.second, row1, col1);
    g.bin.getGridCoords(it2->second.first, it2->second.second, row2, col2);
    // FIX: was `if (x1==y1 && x2==y2) return {x1,y1};` — returning
    // {row, row} (duplicating the row, dropping the column entirely)
    // instead of {row, col} for the equal-points case.
    if (row1 == row2 && col1 == col2) return {row1, col1};
    return get_midpoint(row1, col1, row2, col2);
}
