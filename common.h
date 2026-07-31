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

struct LOSResult
{
    int    x;
    int    y;
    double dist;
    double angle;
    double maxAngle;
    bool   visible;
    int elev;
};
struct RadarQuery
{
    double lat, lon;
    int    heightAbove;
    double radiusMeters;
    bool   earthCurvature;
};

struct LOSQuery
{
    double lat0, lon0; int h0;
    double lat1, lon1; int h1;
};

class Global {
public:
    Bin bin;
    Config& conf = bin.conf;
    Global() {
        bin.init({ ELEVATION });
    }
    public: static void preprocess(bool images = 0, bool timer = false, bool seperate = false, int individual = -1);

    public: std::vector<LOSResult> lineOfVisibilityopt(double lat0, double lon0, int h0, double lat1, double lon1, int h1);

    public: std::vector<LOSResult> lineOfVisibilityOptimized(double lat0, double lon0, int h0, double lat1, double lon1, int h1);

    public: bool lineOfSight(double lat0, double long0,int h0, double lat1, double long1,int h1);

    public: static std::vector<std::pair<int,int>> bresenham(int x0, int y0, int x1, int y1);

    public:  std::vector<LOSResult> lineOfVisibility(double lat0, double long0,int h0, double lat1, double long1,int h1);

    public: static std::vector<std::vector<LOSResult>> lineOfVisibilityBatch(const std::vector<LOSQuery>& queries, int numThreads = 0);  // 0 = auto detect

    public: static inline double deg2rad(double deg) { return deg * M_PI / 180.0; };

    public: static inline double rad2deg(double rad) { return rad * 180.0 / M_PI; };

    public: static double computeDistance(double lonA, double latA, double lonB, double latB, OGRCoordinateTransformation* poCT);

    public: static double computeBearing(double lonA, double latA, double lonB, double latB, OGRCoordinateTransformation* poCT);

    public: static bool dropPointAtDistanceAndBearing(double startLon, double startLat,
                                                      double distanceMeters, double bearingDegrees,
                                                      OGRCoordinateTransformation* poForwardCT,
                                                      OGRCoordinateTransformation* poInverseCT,
                                                      double& outLon, double& outLat);
    static inline FeasibleGrid restrict_to_area_of_operation(Global& g, FeasibleGrid grid) {
        if (g.bin.conf.area_op.size() < 3) return grid; // no polygon defined, nothing to restrict to

        std::vector<std::pair<double,double>> proj_poly;
        proj_poly.reserve(g.bin.conf.area_op.size());
        for (auto& [lat, lon] : g.bin.conf.area_op) {
            double px, py;
            if (g.bin.latlon_to_projected(lat, lon, px, py))
                proj_poly.push_back({px, py});
        }

        std::vector<FeasibleCell> kept;
        kept.reserve(grid.cells.size());
        int dropped = 0;
        for (auto& cell : grid.cells) {
            double cx = g.bin.conf.grid.origin_x + (cell.c + 0.5) * RESOLUTION;
            double cy = g.bin.conf.grid.origin_y - (cell.r + 0.5) * RESOLUTION;
            if (Bin::point_in_polygon(cx, cy, proj_poly)) kept.push_back(cell);
            else dropped++;
        }

        std::cout << "[ao-filter] kept " << kept.size() << " / " << grid.cells.size()
                  << " (" << dropped << " were bounding-box-but-not-AO)\n";
        grid.cells = std::move(kept);
        return grid;
    }

};

#endif
