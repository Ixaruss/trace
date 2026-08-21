#include <debugapi.h>
#define ODLIB_EXPORTS
#include "export.h"
#include "../iad/coverage.h"
using namespace std;

std::vector<std::vector<int>> bin_matrix;

// std::unique_ptr<ThreadPool>pool;
// std::unique_ptr<GaEngine> AdaptDLL;

ThreadPool pool;
GaEngine AdaptDLL(pool);

extern "C" {
ODLIB_API int DLL_Init(){
    // pool = std::make_unique<ThreadPool>();
    // AdaptDLL = std::make_unique<GaEngine>(*pool);
    return 0;
}

ODLIB_API int DLL_TestFunction(){
    OutputDebugStringA("DLL : Testing Successful...");
    return 123456;
}

ODLIB_API int DLL_SurvGridCoverage(bool single_band_flag){
    OutputDebugStringA("DLL : SGrid Coverage Computation Starts...");
    Iad::get_iad_coverage_tif(AdaptDLL.g, OutputCrs::Wgs84, single_band_flag ? BandLayout::SingleBandPalette : BandLayout::Rgba);

    OutputDebugStringA("DLL : SGrid Coverage Ends...");
    return 0;
}

ODLIB_API int DLL_PreprocessMaps(bool ImageFlag, bool IndividualFlag){
    OutputDebugStringA("DLL : Preprocessing Starts...");
    Global::preprocess(ImageFlag,false,IndividualFlag,-1);
    OutputDebugStringA("DLL : Preprocessing Ends...");
    return 0;
}

ODLIB_API int DLL_LoadBinMaps(){
    OutputDebugStringA("DLL : Map Bins Loading Starts...");
    AdaptDLL.g.bin.init();
    OutputDebugStringA("DLL : Map Bins Loading Ends...");
    return 0;
}

ODLIB_API int DLL_RunPSO() {
    return 0;
}

ODLIB_API int DLL_RunAutoGA(int res_reduction_factor, bool single_band_flag) {
    OutputDebugStringA("DLL : GA (Auto Case Starts)...");
    AdaptDLL.reloadConfig();
    AdaptDLL.run_auto(res_reduction_factor, false, OutputCrs::Wgs84, single_band_flag ? BandLayout::SingleBandPalette : BandLayout::Rgba);
    OutputDebugStringA("DLL : GA (Auto Case Ends)...");
    return 0;
}

ODLIB_API int DLL_RunHomoGA(int res_reduction_factor, bool single_band_flag) {
    OutputDebugStringA("DLL : GA (Optimal Deployment - Homogeneous Case Starts)...");
    AdaptDLL.reloadConfig();
    AdaptDLL.run_homogenous(res_reduction_factor, false, OutputCrs::Wgs84, single_band_flag ? BandLayout::SingleBandPalette : BandLayout::Rgba);
    OutputDebugStringA("DLL : GA (Optimal Deployment - Homogeneous Case Ends)...");
    return 0;
}

ODLIB_API int DLL_RunHeteroGA(int res_reduction_factor, bool single_band_flag) {
    OutputDebugStringA("DLL : GA (Optimal Deployment - Heterogeneous Case Starts)...");
    AdaptDLL.reloadConfig();
    AdaptDLL.run_heterogenous(res_reduction_factor, false, OutputCrs::Wgs84, single_band_flag ? BandLayout::SingleBandPalette : BandLayout::Rgba);
    OutputDebugStringA("DLL : GA (Optimal Deployment - Heterogeneous Case Starts)...");
    return 0;
}

ODLIB_API int DLL_ViewShed(double s_lat, double s_lon, int s_hgt_m, int rng_km, int loa, int roa, bool is_sector, bool GeoTifFlag, bool SvgFlag, int res_reduction_factor, bool use_curvature_flag, bool single_band_flag) {
    OutputDebugStringA("DLL : Viewshed Computation Starts...");
    ViewshedEngine VShedDLL(AdaptDLL.g, s_lat,s_lon,rng_km,s_hgt_m, use_curvature_flag);

    if (is_sector){
        VShedDLL.computeViewshed(loa,roa);
        VShedDLL.writeViewshedToPNG("./output/viewshed/viewshed.png");

        if (SvgFlag) VShedDLL.writeViewshedToSVG("./output/viewshed/viewshed.svg");
        if (GeoTifFlag) VShedDLL.writeViewshedToGeoTIFF("./output/viewshed/viewshed.tif", res_reduction_factor, OutputCrs::Wgs84, single_band_flag ? BandLayout::SingleBandPalette : BandLayout::Rgba);
    }
    else{
        VShedDLL.computeViewshed();
        VShedDLL.writeViewshedToPNG("./output/viewshed/viewshed.png");

        if (SvgFlag) VShedDLL.writeViewshedToSVG("./output/viewshed/viewshed.svg");
        if (GeoTifFlag) VShedDLL.writeViewshedToGeoTIFF("./output/viewshed/viewshed.tif", res_reduction_factor, OutputCrs::Wgs84, BandLayout::SingleBandPalette);
    }
    OutputDebugStringA("DLL : Viewshed Computation Ends...");
    return 0;
}

ODLIB_API int DLL_RadarCoverage(double s_lat, double s_lon, double s_hgt_m, double rng_km, double ba_min, double ba_max, double tgt_hgt_m, double loa, double roa, bool is_sector, bool GeoTifFlag, bool SvgFlag, int res_reduction_factor, bool single_band_flag, bool earth_curv_flag) {
    OutputDebugStringA("DLL : Radar Coverage Computation Starts...");

    RadarEngine RadarDLL(AdaptDLL.g, s_lat,s_lon,s_hgt_m, tgt_hgt_m, ba_min, ba_max, rng_km, earth_curv_flag);
    if (is_sector){
        auto res = RadarDLL.compute(loa,roa);
        RadarDLL.writeToPNG(res, "./output/radar/radarcov.png");

        if (SvgFlag) RadarDLL.writeToSVG(res, "./output/radar/radarcov.svg");
        if (GeoTifFlag) RadarDLL.writeCoverageToGeoTIFF(res, "./output/radar/radarcov.tif", res_reduction_factor, OutputCrs::Wgs84, single_band_flag ? BandLayout::SingleBandPalette : BandLayout::Rgba);
    }
    else{
        auto res = RadarDLL.compute();
        RadarDLL.writeToPNG(res, "./output/radar/radarcov.png");

        if (SvgFlag) RadarDLL.writeToSVG(res, "./output/radar/radarcov.svg");
        if (GeoTifFlag) RadarDLL.writeCoverageToGeoTIFF(res, "./output/radar/radarcov.tif", res_reduction_factor, OutputCrs::Wgs84, single_band_flag ? BandLayout::SingleBandPalette : BandLayout::Rgba);
    }
    OutputDebugStringA("DLL : Radar Coverage Computation Ends...");
    return 0;
}

ODLIB_API int DLL_RadarCoverageFDR(double s_lat, double s_lon, int s_hgt_m, int rng_km, int ba_min, int ba_max, int tgt_hgt_m, int res_reduction_factor) {
    OutputDebugStringA("DLL : Radar Coverage Computation (First Detection) Starts...");

    int x,y;
    AdaptDLL.g.bin.getGridCoords(s_lat, s_lon, x, y);

    RadarEngine::get_coverage_reverse(AdaptDLL.g, x, y, s_hgt_m, tgt_hgt_m, ba_min, ba_max, rng_km*1000,"./output/radar/radarfdr.tif", "./output/radar/radarfdr.png", res_reduction_factor);
    OutputDebugStringA("DLL : Radar Coverage Computation(First Detection) Ends...");
    return 0;
}

ODLIB_API int16_t DLL_GetElevation(double lat, double lon) {
    return AdaptDLL.g.bin.getElevation(lat, lon);
}

ODLIB_API double DLL_ComputeDistance(double s_lat, double s_lon, double d_lat, double d_lon) {
    return (Global::computeDistance(s_lat, s_lon, d_lat, d_lon, AdaptDLL.g.bin.poCT));
}

ODLIB_API double DLL_ComputeBearing(double s_lat, double s_lon, double d_lat, double d_lon) {
    return (Global::computeBearing(s_lat, s_lon, d_lat, d_lon, AdaptDLL.g.bin.poCT));
}

ODLIB_API int DLL_GetStatusRoad(double lat, double lon) {
    return AdaptDLL.g.bin.getStatusRoad(lat, lon);
}

ODLIB_API int DLL_GetStatusRail(double lat, double lon) {
    return AdaptDLL.g.bin.getStatusRail(lat, lon);
}

ODLIB_API int DLL_GetStatusWaterArea(double lat, double lon) {
    return AdaptDLL.g.bin.getStatusWaterArea(lat, lon);
}

ODLIB_API int DLL_GetStatusWaterLine(double lat, double lon) {
    return AdaptDLL.g.bin.getStatusWaterLine(lat, lon);
}

ODLIB_API int DLL_GetStatusIB(double lat, double lon) {
    return AdaptDLL.g.bin.getStatusIB(lat, lon);
}
}
