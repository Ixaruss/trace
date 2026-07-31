#ifndef EXPORT_H_
#define EXPORT_H_
#include "../cli/cli.h"
#include "../common.h"
#include "../operations.h"
#include "../ga/config.h"
#include "../ga/ga.h"
#include <cstdint>

#ifdef _WIN32
    #ifdef ODLIB_EXPORTS
        #define ODLIB_API extern "C" __declspec(dllexport)
    #else
        #define ODLIB_API extern "C" __declspec(dllimport)
    #endif
#else
    #define ODLIB_API extern "C" __attribute__((visibility("default")))
#endif

ODLIB_API int       DLL_TestFunction();

ODLIB_API int       DLL_PreprocessMaps(bool,bool);
ODLIB_API int       DLL_LoadBinMaps();

ODLIB_API int       DLL_RunPSO();
ODLIB_API int       DLL_RunAutoGA(int);
ODLIB_API int       DLL_RunHomoGA(int);
ODLIB_API int       DLL_RunHeteroGA(int);

ODLIB_API int       DLL_ViewShed(double, double, int, int, int, int, bool, bool, bool, int, bool);
ODLIB_API int       DLL_RadarCoverage(double, double, double, double, double, double, double, double, double, bool, bool, bool, int);
ODLIB_API int       DLL_RadarCoverageFDR(double, double, int, int, int, int, int, int);

ODLIB_API int16_t   DLL_GetElevation(double, double);
ODLIB_API double    DLL_ComputeDistance(double, double,double, double);
ODLIB_API double    DLL_ComputeBearing(double, double,double, double);

ODLIB_API int       DLL_GetStatusRoad(double, double);
ODLIB_API int       DLL_GetStatusRail(double, double);
ODLIB_API int       DLL_GetStatusWaterArea(double, double);
ODLIB_API int       DLL_GetStatusWaterLine(double, double);
ODLIB_API int       DLL_GetStatusIB(double, double);

#endif
