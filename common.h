#ifndef SHARED_H_
#define SHARED_H_


#include "preprocess/Config.h"
#include "preprocess/bin.h"
#include "preprocess/lib.h"
#include <sys/stat.h>

#if defined(_MSC_VER) || defined(__MINGW32__)
    // Windows (MSVC or MinGW)
    #define fseek64 _fseeki64
    #define ftell64 _ftelli64
    typedef __int64 off64_t;
#elif defined(__unix__) || defined(__APPLE__) || defined(__linux__)
    // POSIX systems
    #define fseek64 fseeko
    #define ftell64 ftello
#else
    // Fallback for older compilers (may fallback to 32-bit)
    #define fseek64 fseek
    #define ftell64 ftell
#endif

struct RadarQuery
{
    double lat, lon;
    int    heightAbove;
    double radiusMeters;
    bool   earthCurvature;
};

class Global {
public:
    Bin bin;
    Config& conf = bin.conf;

    // No auto-load: Global g; gives you a Bin with nothing resolved yet.
    // Each call site loads exactly what it needs, e.g.:
    //   g.load({ELEVATION}, UNBOUNDED, lat, lng, rangeKm);
    //   g.load({CONSTRAINTS}, AOI);
    //   g.load({ALL}, FULL);
    void load(std::set<LoadOptions> opt,
             OperationType type       = FULL,
             double        center_lat = 0.0,
             double        center_lng = 0.0,
             double        range_km   = 0.0) {
        bin.init(opt, type, center_lat, center_lng, range_km);
    }

    // Move this Global's Bin to a different region without reconstructing
    // it. IMPORTANT: any engine (ViewshedEngine/RadarEngine) built against
    // the old window is now invalid — its gridCols/gridRows, translated
    // centre, and mask/status buffers all refer to the region that was
    // just unloaded. Rebuild engines after calling this, don't reuse them.
    void reload(std::set<LoadOptions> opt,
               OperationType type       = FULL,
               double        center_lat = 0.0,
               double        center_lng = 0.0,
               double        range_km   = 0.0) {
        bin.reload(opt, type, center_lat, center_lng, range_km);
    }

    public: static void preprocess(bool images = 0, bool timer = false, bool seperate = false, int individual = -1);

    // The one real two-point LOS check (common/lineofsight.cc) -- correct
    // for windowed/tiled loading via Bin::getElevationRC(). Its former
    // siblings (lineOfVisibilityopt/lineOfVisibility/lineOfVisibilityOptimized/
    // lineOfVisibilityBatch) were duplicate implementations that never got
    // the same fix and have been removed; radar-specific LOS lives in
    // RadarEngine::radarLOS/radarLOSDetailed (operations.h / radarlos.cc)
    // instead, used by the CLI's `rlos` subcommand.
    public: bool lineOfSight(double lat0, double long0,int h0, double lat1, double long1,int h1);

    public: static std::vector<std::pair<int,int>> bresenham(int x0, int y0, int x1, int y1);

    public: static inline double deg2rad(double deg) { return deg * M_PI / 180.0; };

    public: static inline double rad2deg(double rad) { return rad * 180.0 / M_PI; };

    public: static double computeDistance(double lonA, double latA, double lonB, double latB, OGRCoordinateTransformation* poCT);

    public: static double computeBearing(double lonA, double latA, double lonB, double latB, OGRCoordinateTransformation* poCT);

    public: static bool dropPointAtDistanceAndBearing(double startLon, double startLat,
                                                      double distanceMeters, double bearingDegrees,
                                                      OGRCoordinateTransformation* poForwardCT,
                                                      OGRCoordinateTransformation* poInverseCT,
                                                      double& outLon, double& outLat);

};

#endif
