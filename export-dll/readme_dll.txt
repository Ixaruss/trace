# compile code to create adapt.dll
g++ -shared -O3 -ffast-math -flto -ffp-contract=fast -fopenmp  -march=native -funroll-loops -fprefetch-loop-arrays -fdata-sections -ffunction-sections -Wl,--gc-sections -Wl,-O2 cli/main.cc common/lineofsight.cc common/image.cc common/coverage.cc common/angleviewshed.cc common/radarlos.cc common/mult.cc  common/operations.cc preprocess/jsonParser.cc preprocess/reproject.cc export-dll/export.cc preprocess/preprocess.cc preprocess/lib/*.cc ga/config.cc -o adapt.dll  -ljsoncpp -lgdal -lpng16 -lhdf5 -std=c++23

# check dll functions
# open developer command for vs2017
cd C:\msys64\home\STEP-LE3\optimal-deployment
dumpbin /exports adapt.dll

# check dependencies
dumpbin /dependents adapt.dll

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# working export.cc
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#define ODLIB_EXPORTS
#include "export.h"

#include "../preprocess/bin.h"

using namespace std;

int TestHelloWorld(){
    return 123456;
}

int Bin_BuildElevationMatrix() {
    return Bin::build_elevation_matrix();
}
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# working export.cc
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# working export.h
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#ifndef EXPORT_H_
#define EXPORT_H_

#include <cstdint>

// ── symbol visibility ──────────────────────────────────────────────────────
#ifdef _WIN32
  #ifdef ODLIB_EXPORTS
    #define ODLIB_API extern "C" __declspec(dllexport)
  #else
    #define ODLIB_API extern "C" __declspec(dllimport)
  #endif
#else
  #define ODLIB_API extern "C" __attribute__((visibility("default")))
#endif


ODLIB_API int       TestHelloWorld();
ODLIB_API int       Bin_BuildElevationMatrix ();
// ODLIB_API int       MapsPreprocessing();
#endif
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# working export.h
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# working dll import in visual studio code
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
<DllImport("C:\msys64\home\STEP-LE3\optimal-deployment\adapt.dll", CallingConvention:=CallingConvention.Cdecl, SetLastError:=True)>
Public Shared Function Bin_BuildElevationMatrix() As Integer
End Function

<DllImport("C:\msys64\home\STEP-LE3\optimal-deployment\adapt.dll", CallingConvention:=CallingConvention.Cdecl, SetLastError:=True)>
Public Shared Function TestHelloWorld() As Integer
End Function
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# working dll import in visual studio code
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# function call in visual studio code
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Dim x As Integer = TestHelloWorld()
Console.WriteLine("value of x:")
Console.WriteLine(x)

Dim y As Integer = Bin_BuildElevationMatrix()
Console.WriteLine("\nvalue of y:")
Console.WriteLine(y)
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# function call in visual studio code
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@





















3. g++ -shared -O3 -ffast-math -flto -ffp-contract=fast -fopenmp  -march=native -funroll-loops -fprefetch-loop-arrays -fdata-sections -ffunction-sections -Wl,--gc-sections -Wl,-O2 cli/main.cc common/lineofsight.cc common/image.cc common/coverage.cc common/angleviewshed.cc common/radarlos.cc common/mult.cc  common/operations.cc preprocess/jsonParser.cc preprocess/reproject.cc export-dll//export_adapt.cc -o adapt.dll preprocess/preprocess.cc preprocess/lib/*.cc -o adapt.dll  -ljsoncpp -lgdal -lpng16 -lhdf5 -std=c++23
g++ -shared export-dll//export_adapt.cc  -ljsoncpp -lgdal -lpng16 -lhdf5 -std=c++23

g++ -shared -O3 -ffast-math -flto -ffp-contract=fast -fopenmp  -march=native -funroll-loops -fprefetch-loop-arrays -fdata-sections -ffunction-sections -Wl,--gc-sections -Wl,-O2 export-dll//export_adapt.cc preprocess/preprocess.cc preprocess/lib/*.cc -o adapt.dll  -ljsoncpp -lgdal -lpng16 -lhdf5 -std=c++23



g++ -shared export-dll//export.cc -o adapt.dll  -ljsoncpp -lgdal -lpng16 -lhdf5 -std=c++23


g++ -shared -O3 -ffast-math -flto -ffp-contract=fast -fopenmp  -march=native -funroll-loops -fprefetch-loop-arrays -fdata-sections -ffunction-sections -Wl,--gc-sections -Wl,-O2 export-dll/export.cc -o adapt.dll  -ljsoncpp -lgdal -lpng16 -lhdf5 -std=c++23


g++ -shared -O3 -ffast-math -flto -ffp-contract=fast -fopenmp  -march=native -funroll-loops -fprefetch-loop-arrays -fdata-sections -ffunction-sections -Wl,--gc-sections -Wl,-O2 common/lineofsight.cc common/image.cc common/coverage.cc \
		common/angleviewshed.cc common/radarlos.cc common/mult.cc common/operations.cc \
		preprocess/jsonParser.cc preprocess/reproject.cc preprocess/preprocess.cc \
		preprocess/lib/*.cc ga/config.cc export-dll/export.cc -o adapt.dll -ljsoncpp -lgdal -lpng16 -lhdf5 -std=c++23
