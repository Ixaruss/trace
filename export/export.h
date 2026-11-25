#ifndef EXPORT_H_
#define EXPORT_H_
#include "../cli/cli.h"
#include "../common.h"
#include "../operations.h"
#include "../multi-radar-deployment/config.h"
#include "../multi-radar-deployment/ga.h"

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
ODLIB_API int       DLL_Init();
ODLIB_API int       DLL_UnloadBin(bool);

ODLIB_API int       DLL_RunPSO();
ODLIB_API int       DLL_RunAutoGA(int, bool);
ODLIB_API int       DLL_RunHomoGA(int, bool);
ODLIB_API int       DLL_RunHeteroGA(int, bool);

ODLIB_API int       DLL_ViewShed(double, double, int, int, int, int, bool, bool, bool, int, bool, bool);
ODLIB_API int       DLL_RadarCoverage(int, bool);

ODLIB_API int16_t   DLL_GetElevation(double, double);
ODLIB_API double    DLL_ComputeDistance(double, double,double, double);
ODLIB_API double    DLL_ComputeBearing(double, double,double, double);

ODLIB_API bool       DLL_GetStatusRoad(double, double);
ODLIB_API bool       DLL_GetStatusRail(double, double);
ODLIB_API bool       DLL_GetStatusWaterArea(double, double);
ODLIB_API bool       DLL_GetStatusWaterLine(double, double);
ODLIB_API bool       DLL_GetStatusIB(double, double);

#endif
