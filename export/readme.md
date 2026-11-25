# export

Wraps the CLI's engines (preprocessing, viewshed, radar coverage, GA siting) as a
flat C API (`extern "C"`, `ODLIB_API`/`__declspec(dllexport)` on Windows) so an
external host app (e.g. a .NET or C++ application) can call into `trace` as a
DLL instead of shelling out to the CLI. The public contract is
[export.h](export.h); the implementation is [export.cc](export.cc).

## Building

```sh
make dll
```

See the `dll` target in the repo-root [Makefile](../Makefile) for the exact
compiler invocation. It needs a MinGW-flavored toolchain — `export.cc` includes
`<debugapi.h>`, a Windows-only header — so build it under MSYS2, or see
[docs/BUILDING.md](../docs/BUILDING.md#windows) for the Windows toolchain setup.

## Exposed functions

All declared in [export.h](export.h):

| Function | Does |
|---|---|
| `DLL_TestFunction()` | Returns `123456` — call this first to confirm the DLL loaded and is callable. |
| `DLL_Init()` | One-time setup; call before anything else. |
| `DLL_UnloadBin(bool dropVec)` | Releases the loaded terrain/constraint grid. `dropVec` is passed straight through to `Bin::unload()` — also frees the in-memory grid vector, not just file handles, when true. |
| `DLL_PreprocessMaps(bool withImages, bool individualLayers)` | Builds the `.bin` files from `config/config.json` (same as CLI `trace build`). |
| `DLL_RunAutoGA(int reduction, bool singleBand)` | Multi-radar GA siting, "auto" case from `config/gaconf.json`. |
| `DLL_RunHomoGA(int reduction, bool singleBand)` | Multi-radar GA siting, "homogenous" case. |
| `DLL_RunHeteroGA(int reduction, bool singleBand)` | Multi-radar GA siting, "heterogenous" case. |
| `DLL_RunPSO()` | Declared for the single-radar PSO case; currently a stub that returns `0` and does nothing. |
| `DLL_ViewShed(lat, lon, antennaM, rangeKm, loa, roa, isSector, writeTif, writeSvg, reduction, useCurvature, singleBand)` | Viewshed from one point; always writes a PNG, optionally also SVG/GeoTIFF. |
| `DLL_RadarCoverage(int reduction, bool singleBand)` | Single fixed-position radar coverage from `config/radar.json` (same as CLI `trace radar`). |
| `DLL_GetElevation(lat, lon)` | Terrain elevation at a point. |
| `DLL_ComputeDistance(lat1, lon1, lat2, lon2)` | Great-circle-ish distance between two points, in the loaded grid's projection. |
| `DLL_ComputeBearing(lat1, lon1, lat2, lon2)` | Bearing from point 1 to point 2. |
| `DLL_GetStatusRoad/Rail/WaterArea/WaterLine/IB(lat, lon)` | Whether a point falls inside that constraint layer's buffer. |

One function exists the other way around: `DLL_GetLineOfSight(lat0, lon0, h0,
lat1, lon1, h1, useCurvature)` is **implemented** in `export.cc` but not
declared in `export.h` — it's still exported from the compiled DLL (same
`ODLIB_API` linkage), just undocumented to anyone only reading the header.

## Inspecting the built DLL

From a Visual Studio developer command prompt:

```
dumpbin /exports trace.dll
dumpbin /dependents trace.dll
```
