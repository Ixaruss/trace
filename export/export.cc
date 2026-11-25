#include <debugapi.h>
#include <exception>
#include <fstream>
#define ODLIB_EXPORTS
#include "export.h"
#include "../common/config.h"
using namespace std;

std::vector<std::vector<int>> bin_matrix;

ThreadPool pool;
GaEngine TraceDLL(pool);

extern "C" {
ODLIB_API int DLL_Init(){
    try{
        OutputDebugStringA("DLL : Init Begins...");
        OutputDebugStringA("DLL : Init Ends...");
        return 0;
    } catch (const std::exception& e){
            std::ofstream("C:\\tmp\\trace\\error.log", std::ios::app)<<"DLL_Init: " << typeid(e).name() << ": " << e.what() << "\n";
            return -1;
    } catch (...) {
        std::ofstream("C:\\tmp\\trace\\error.log", std::ios::app)<<"DLL_Init: unknown exception\n";
        return -2;
    }
}

ODLIB_API int DLL_TestFunction(){
    OutputDebugStringA("DLL : Testing Successful...");
    return 123456;
}

ODLIB_API int DLL_UnloadBin(bool single_band_flag){
    OutputDebugStringA("DLL : Unloading Bin Starts...");
    TraceDLL.g.bin.unload(single_band_flag);
    OutputDebugStringA("DLL : Unloading Bin Ends...");
    return 0;
}

ODLIB_API int DLL_PreprocessMaps(bool ImageFlag, bool IndividualFlag){
    OutputDebugStringA("DLL : Preprocessing Starts...");
    Global::preprocess(ImageFlag, false, IndividualFlag,-1);
    OutputDebugStringA("DLL : Preprocessing Ends...");
    return 0;
}

ODLIB_API int DLL_RunPSO() {
    return 0;
}

ODLIB_API int DLL_RunAutoGA(int res_reduction_factor, bool single_band_flag) {
    OutputDebugStringA("DLL : GA (Auto Case Starts)...");
    TraceDLL.g.bin.unload(true);
    TraceDLL.reloadConfig();
    TraceDLL.run_auto(res_reduction_factor, false, OutputCrs::Wgs84, single_band_flag ? BandLayout::SingleBandPalette : BandLayout::Rgba);
    OutputDebugStringA("DLL : GA (Auto Case Ends)...");
    return 0;
}

ODLIB_API int DLL_RunHomoGA(int res_reduction_factor, bool single_band_flag) {
    OutputDebugStringA("DLL : GA (Optimal Deployment - Homogeneous Case Starts)...");
    TraceDLL.g.bin.unload(true);
    TraceDLL.reloadConfig();
    TraceDLL.run_homogenous(res_reduction_factor, false, OutputCrs::Wgs84, single_band_flag ? BandLayout::SingleBandPalette : BandLayout::Rgba);
    OutputDebugStringA("DLL : GA (Optimal Deployment - Homogeneous Case Ends)...");
    return 0;
}

ODLIB_API int DLL_RunHeteroGA(int res_reduction_factor, bool single_band_flag) {
    OutputDebugStringA("DLL : GA (Optimal Deployment - Heterogeneous Case Starts)...");
    TraceDLL.g.bin.unload(true);
    TraceDLL.reloadConfig();
    TraceDLL.run_heterogenous(res_reduction_factor, false, OutputCrs::Wgs84, single_band_flag ? BandLayout::SingleBandPalette : BandLayout::Rgba);
    OutputDebugStringA("DLL : GA (Optimal Deployment - Heterogeneous Case Starts)...");
    return 0;
}

ODLIB_API int DLL_ViewShed(double s_lat, double s_lon, int s_hgt_m, int rng_km, int loa, int roa, bool is_sector, bool GeoTifFlag, bool SvgFlag, int res_reduction_factor, bool use_curvature_flag, bool single_band_flag) {
    OutputDebugStringA("DLL : Viewshed Computation Starts...");
    TraceDLL.g.bin.unload(true);
    TraceDLL.g.load({ELEVATION}, OperationType::UNBOUNDED, s_lat, s_lon, rng_km);
    ViewshedEngine VShedDLL(TraceDLL.g, s_lat,s_lon,rng_km,s_hgt_m, use_curvature_flag);

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

ODLIB_API int DLL_RadarCoverage(int res_reduction_factor, bool single_band_flag) {
    OutputDebugStringA("DLL : Radar Coverage Computation Starts...");

    if(!Rconf::is_available()){
        std::cerr<<"[ERROR] radar conf not found\n";
        return -1;
    }
    auto rconf = Rconf::load();

    OutputDebugStringA("DLL : x...");
    TraceDLL.g.bin.unload(true);
    TraceDLL.g.load({ELEVATION}, OperationType::UNBOUNDED, rconf.srcLat, rconf.srcLon, rconf.range);
    OutputDebugStringA("DLL : y...");

    RadarEngine RadarDLL(TraceDLL.g);
    auto res = RadarDLL.compute();
    string output_ext = rconf.output.substr(rconf.output.size()-3, 3);
    string path = rconf.output.substr(0,rconf.output.size()-3);

    if(output_ext == "tif" ) RadarDLL.writeCoverageToGeoTIFF(res, "./output/radar/radarcov.tif", res_reduction_factor, OutputCrs::Wgs84, single_band_flag ? BandLayout::SingleBandPalette : BandLayout::Rgba);
    else if(output_ext == "svg") RadarDLL.writeToSVG(res, "./output/radar/radarcov.svg");
    else RadarDLL.writeToPNG(res, "./output/radar/radarcov.png");;
    string r_output = path + output_ext;

    OutputDebugStringA("DLL : Radar Coverage Computation Ends...");
    return 0;
}

ODLIB_API int16_t DLL_GetElevation(double lat, double lon) {
    return TraceDLL.g.bin.getElevation(lat, lon);
}

ODLIB_API double DLL_ComputeDistance(double s_lat, double s_lon, double d_lat, double d_lon) {
    return (Global::computeDistance(s_lat, s_lon, d_lat, d_lon, TraceDLL.g.bin.poCT));
}

ODLIB_API double DLL_ComputeBearing(double s_lat, double s_lon, double d_lat, double d_lon) {
    return (Global::computeBearing(s_lat, s_lon, d_lat, d_lon, TraceDLL.g.bin.poCT));
}

ODLIB_API bool DLL_GetStatusRoad(double lat, double lon) {
    TraceDLL.g.load({CONSTRAINTS},UNBOUNDED,lat,lon,0.1);
    return TraceDLL.g.bin.getStatusRoad(lat, lon);
}

ODLIB_API bool DLL_GetStatusRail(double lat, double lon) {
    TraceDLL.g.load({CONSTRAINTS},UNBOUNDED,lat,lon,0.1);
    return TraceDLL.g.bin.getStatusRail(lat, lon);
}

ODLIB_API bool DLL_GetStatusWaterArea(double lat, double lon) {
    TraceDLL.g.load({CONSTRAINTS},UNBOUNDED,lat,lon,0.1);
    return TraceDLL.g.bin.getStatusWaterArea(lat, lon);
}

ODLIB_API bool DLL_GetStatusWaterLine(double lat, double lon) {
    TraceDLL.g.load({CONSTRAINTS},UNBOUNDED,lat,lon,0.1);
    return TraceDLL.g.bin.getStatusWaterLine(lat, lon);
}

ODLIB_API bool DLL_GetStatusIB(double lat, double lon) {
    TraceDLL.g.load({CONSTRAINTS},UNBOUNDED,lat,lon,0.1);
    return TraceDLL.g.bin.getStatusIB(lat, lon);
}

ODLIB_API bool DLL_GetLineOfSight(double lat0, double long0, int h0, double lat1, double long1, int h1, bool use_curv) {
    return TraceDLL.g.lineOfSight(lat0, long0, h0, lat1, long1, h1, use_curv, false);
}

}
