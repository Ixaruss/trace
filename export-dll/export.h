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

ODLIB_API int       ADAPT_TestFunction();
ODLIB_API int       ADAPT_PreprocessMaps();

ODLIB_API int       ADAPT_RunPSO();
ODLIB_API int       ADAPT_RunAutoGA();
ODLIB_API int       ADAPT_RunHomoGA();
ODLIB_API int       ADAPT_RunHeteroGA();

ODLIB_API int       ADAPT_ViewShed(double, double, int, int, int, int, bool);
ODLIB_API int       ADAPT_RadarCoverage(double, double, int, int, int, int, int, int, int, bool);
ODLIB_API int       ADAPT_RadarCoverageFDR(double, double, int, int, int, int, int);

ODLIB_API int16_t   ADAPT_GetElevation(double, double);
ODLIB_API double    ADAPT_ComputeDistance(double, double,double, double);
ODLIB_API double    ADAPT_ComputeBearing(double, double,double, double);

ODLIB_API int       ADAPT_GetStatusRoad(double, double);
ODLIB_API int       ADAPT_GetStatusRail(double, double);
ODLIB_API int       ADAPT_GetStatusWaterArea(double, double);
ODLIB_API int       ADAPT_GetStatusWaterLine(double, double);
ODLIB_API int       ADAPT_GetStatusIB(double, double);

ODLIB_API int       ADAPT_GetGridInfo(double, double, int *, int *);
ODLIB_API int       ADAPT_GetLatLonInfo(int, int, double *, double *);

ODLIB_API int       ADAPT_DisplayElevationBin();
ODLIB_API int       ADAPT_LoadElevationBin();



#endif // EXPORT_ADAPT_H
