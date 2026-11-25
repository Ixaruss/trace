CXX      := g++
STD      := -std=c++23

# NOTE: original Makefile used $(gdal-config --cflags --libs) which is not
# valid Make syntax (missing the `shell` function) and silently expanded to
# nothing. Fixed below.
GDAL     := $(shell gdal-config --cflags --libs)
LIBS     := -ljsoncpp -lgdal -lpng16 -lhdf5
HDF5     := -lhdf5

BIN      := x
TRACE    := trace

# Shared preprocessing sources used by every target
COMMON_SRC := preprocess/jsonParser.cc preprocess/reproject.cc preprocess/lib/*.cc common/coverage.cc render.cc

# Aggressive optimization flags, shared by the perf-sensitive targets (tool, multi-radar-deployment)
OPT_FLAGS := -O3 -ffast-math -flto -ffp-contract=fast -march=native \
             -funroll-loops -fprefetch-loop-arrays -fdata-sections \
             -ffunction-sections -Wl,--gc-sections -Wl,-O2

.PHONY: all bin test tool deploy multi-radar-deployment clean

all: tool

test:
	$(CXX) preprocess/test.cc $(COMMON_SRC) -o $(BIN) $(GDAL) $(LIBS) $(HDF5) $(STD) && ./$(BIN)

tool:
	$(CXX) $(OPT_FLAGS) -fopenmp \
		cli/main.cc common/lineofsight.cc common/image.cc common/coverage.cc render.cc \
		common/angleviewshed.cc common/radarlos.cc common/config.cc common/operations.cc \
		preprocess/jsonParser.cc preprocess/reproject.cc preprocess/preprocess.cc \
		preprocess/lib/*.cc multi-radar-deployment/config.cc \
		-o $(TRACE) $(LIBS) $(STD)

multi-radar-deployment:
	$(CXX) $(OPT_FLAGS) -pthread \
		multi-radar-deployment/main.cc multi-radar-deployment/config.cc $(COMMON_SRC) \
		-o $(BIN) $(GDAL) $(LIBS) $(HDF5) $(STD) && ./$(BIN)

clean:
	rm -rf output bin && mkdir -p output/{ga,gbadws,heatmap,path,preprocess,pso,radar,simdata,viewshed}

dll:
	g++ -shared -O3 -ffast-math -flto -ffp-contract=fast -fopenmp  -march=native -funroll-loops -fprefetch-loop-arrays -fdata-sections -ffunction-sections -Wl,--gc-sections -Wl,-O2 common/lineofsight.cc common/image.cc common/coverage.cc common/angleviewshed.cc common/radarlos.cc common/operations.cc preprocess/jsonParser.cc preprocess/reproject.cc export/export.cc preprocess/preprocess.cc render.cc preprocess/lib/*.cc multi-radar-deployment/config.cc common/lib.cc  common/config.cc -o trace.dll  -ljsoncpp -lgdal -lpng16 -lhdf5 -static-libgcc -static-libstdc++ -Wl,--kill-at -std=c++23
