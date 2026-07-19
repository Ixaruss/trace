#include <debugapi.h>
#define ODLIB_EXPORTS
#include "export.h"
using namespace std;

std::vector<std::vector<int>> bin_matrix;
ThreadPool pool;
GaEngine AdaptDLL(pool);

ODLIB_API int ADAPT_TestFunction(){
    OutputDebugStringA("DLL : Testing Successful...");
    return 123456;
}

ODLIB_API int ADAPT_PreprocessMaps(){
    OutputDebugStringA("DLL : Preprocessing Begins...");
    Global::preprocess(true,false,true,-1);
    OutputDebugStringA("DLL : Preprocessing Ends...");
    return 0;
}

extern "C" {
ODLIB_API int ADAPT_LoadElevationBin(){
    OutputDebugStringA("DLL : 2D Matrix Data Loading Ends...");
    AdaptDLL.g.bin.init();
    return 0;
}
}

ODLIB_API int ADAPT_RunPSO() {
    return 0;
}

ODLIB_API int ADAPT_RunAutoGA() {
    AdaptDLL.run_auto();
    return 0;
}

ODLIB_API int ADAPT_RunHomoGA() {
    AdaptDLL.run_homogenous();
    return 0;
}

ODLIB_API int ADAPT_RunHeteroGA() {
    AdaptDLL.run_heterogenous();
    return 0;
}

ODLIB_API int ADAPT_ViewShed(double s_lat, double s_lon, int s_hgt_m, int rng_km, int loa, int roa, bool is_sector) {
    OutputDebugStringA("DLL : Viewshed Computation Begins...");
    ViewshedEngine VShedDLL(AdaptDLL.g, s_lat,s_lon,rng_km,s_hgt_m,true);
    if (is_sector){
        VShedDLL.computeViewshed(loa,roa);
        VShedDLL.writeViewshedToPNG("./output/viewshed/vs_sector.svg");
    }
    else{
        VShedDLL.computeViewshed();
        VShedDLL.writeViewshedToPNG("./output/viewshed/vs_full.svg");
    }
    OutputDebugStringA("DLL : Viewshed Computation Ends...");
    return 0;
}


ODLIB_API int ADAPT_RadarCoverage(double s_lat, double s_lon, int s_hgt_m, int rng_km, int ba_min, int ba_max, int tgt_hgt_m, int loa, int roa, bool is_sector) {
    RadarEngine RadarDLL(AdaptDLL.g, s_lat,s_lon,s_hgt_m, tgt_hgt_m, ba_min, ba_max, rng_km,true);
    if (is_sector){
        auto res = RadarDLL.compute(loa,roa);
        RadarDLL.writeToPNG(res, "/output/viewshed/rcov_sector.svg");
    }
    else{
        auto res = RadarDLL.compute();
        RadarDLL.writeToPNG(res, "/output/viewshed/rcov_full.svg");
    }
    return 0;
}

ODLIB_API int ADAPT_RadarCoverageFDR(double s_lat, double s_lon, int s_hgt_m, int rng_km, int ba_min, int ba_max, int tgt_hgt_m) {
    int x,y;
    AdaptDLL.g.bin.getGridCoords(s_lat, s_lon, x, y);

    RadarEngine::get_coverage_reverse(AdaptDLL.g, x, y, s_hgt_m, tgt_hgt_m, ba_min, ba_max, rng_km*1000, "./output/radar/fdr.png");
    return 0;
}

ODLIB_API int16_t ADAPT_GetElevation(double lat, double lon) {
    return AdaptDLL.g.bin.getElevation(lat, lon);
}

ODLIB_API double ADAPT_ComputeDistance(double s_lat, double s_lon, double d_Lat, double d_lon) {
    return 0;
}

ODLIB_API double ADAPT_ComputeBearing(double s_lat, double s_lon, double d_Lat, double d_lon) {
    return 0;
}

ODLIB_API int ADAPT_GetStatusRoad(double lat, double lon) {
    return AdaptDLL.g.bin.getStatusRoad(lat, lon);
}

ODLIB_API int ADAPT_GetStatusRail(double lat, double lon) {
    return AdaptDLL.g.bin.getStatusRail(lat, lon);
}

ODLIB_API int ADAPT_GetStatusWaterArea(double lat, double lon) {
    return AdaptDLL.g.bin.getStatusWaterArea(lat, lon);
}

ODLIB_API int ADAPT_GetStatusWaterLine(double lat, double lon) {
    return AdaptDLL.g.bin.getStatusWaterLine(lat, lon);
}

ODLIB_API int ADAPT_GetStatusIB(double lat, double lon) {
    return AdaptDLL.g.bin.getStatusIB(lat, lon);
}
