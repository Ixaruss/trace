#include <iostream>
#include <cmath>
#include <ogr_geometry.h>
#include <ogr_spatialref.h>
#include "../common.h"


// =========================================================================
// FUNCTION 1: Compute Distance Between Two Lat/Lon Points
// =========================================================================
double Global::computeDistance(double lonA, double latA, double lonB, double latB, OGRCoordinateTransformation* poCT) {
    OGRPoint pointA(lonA, latA);
    OGRPoint pointB(lonB, latB);

    // Transform points to the metric projection (e.g., EPSG:7755)
    if (pointA.transform(poCT) != OGRERR_NONE || pointB.transform(poCT) != OGRERR_NONE) {
        std::cerr << "Distance transformation failed." << std::endl;
        return -1.0;
    }

    // Returns distance in meters
    return pointA.Distance(&pointB);
}

// =========================================================================
// FUNCTION 2: Compute Bearing Between Two Lat/Lon Points
// =========================================================================
double Global::computeBearing(double lonA, double latA, double lonB, double latB, OGRCoordinateTransformation* poCT) {
    OGRPoint pointA(lonA, latA);
    OGRPoint pointB(lonB, latB);

    // Transform points to the metric projection
    if (pointA.transform(poCT) != OGRERR_NONE || pointB.transform(poCT) != OGRERR_NONE) {
        std::cerr << "Bearing transformation failed." << std::endl;
        return -1.0;
    }

    double deltaX = pointB.getX() - pointA.getX();
    double deltaY = pointB.getY() - pointA.getY();

    // Calculate angle clockwise from North (Y-axis)
    double bearingDeg = rad2deg(std::atan2(deltaX, deltaY));

    // Normalize to 0° - 360°
    if (bearingDeg < 0) {
        bearingDeg += 360.0;
    }

    return bearingDeg;
}

// =========================================================================
// FUNCTION 3: Find Destination Point from Source, Distance, and Bearing
// =========================================================================
bool Global::dropPointAtDistanceAndBearing(double startLon, double startLat,
                                   double distanceMeters, double bearingDegrees,
                                   OGRCoordinateTransformation* poForwardCT,
                                   OGRCoordinateTransformation* poInverseCT,
                                   double& outLon, double& outLat) {
    OGRPoint pointA(startLon, startLat);

    // Transform start point to metric projection
    if (pointA.transform(poForwardCT) != OGRERR_NONE) {
        std::cerr << "Forward projection failed for starting point." << std::endl;
        return false;
    }

    double bearingRad = deg2rad(bearingDegrees);

    // Calculate new Cartesian coordinates
    double newX = pointA.getX() + (distanceMeters * std::sin(bearingRad));
    double newY = pointA.getY() + (distanceMeters * std::cos(bearingRad));

    OGRPoint pointC(newX, newY);

    // Transform back to WGS84 Lat/Lon
    if (pointC.transform(poInverseCT) != OGRERR_NONE) {
        std::cerr << "Inverse transformation failed for destination point." << std::endl;
        return false;
    }

    outLon = pointC.getX();
    outLat = pointC.getY();
    return true;
}
