# compile code to create adapt.dll
g++ -shared -O3 -ffast-math -flto -ffp-contract=fast -fopenmp  -march=native -funroll-loops -fprefetch-loop-arrays -fdata-sections -ffunction-sections -Wl,--gc-sections -Wl,-O2 common/lineofsight.cc common/image.cc common/coverage.cc common/angleviewshed.cc common/radarlos.cc common/mult.cc  common/operations.cc preprocess/jsonParser.cc preprocess/reproject.cc export-dll/export.cc preprocess/preprocess.cc render.cc common/radarfdr.cc preprocess/lib/*.cc ga/config.cc common/lib.cc -o adapt.dll  -ljsoncpp -lgdal -lpng16 -lhdf5 -static-libgcc -static-libstdc++ -Wl,--kill-at -std=c++23

# check dll functions
# open developer command for vs2017
cd C:\msys64\home\STEP-LE3\optimal-deployment
dumpbin /exports adapt.dll

# check dependencies
dumpbin /dependents adapt.dll
