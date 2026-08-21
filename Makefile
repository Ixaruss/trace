CXX      := g++
STD      := -std=c++23

# NOTE: original Makefile used $(gdal-config --cflags --libs) which is not
# valid Make syntax (missing the `shell` function) and silently expanded to
# nothing. Fixed below.
GDAL     := $(shell gdal-config --cflags --libs)
LIBS     := -ljsoncpp -lgdal -lpng16 -lhdf5
HDF5     := -lhdf5

BIN      := x
ADAPT    := adapt

# Shared preprocessing sources used by every target
COMMON_SRC := preprocess/jsonParser.cc preprocess/reproject.cc preprocess/lib/*.cc common/radarfdr.cc  common/coverage.cc render.cc

# Aggressive optimization flags, shared by the perf-sensitive targets (tool, ga)
OPT_FLAGS := -O3 -ffast-math -flto -ffp-contract=fast -march=native \
             -funroll-loops -fprefetch-loop-arrays -fdata-sections \
             -ffunction-sections -Wl,--gc-sections -Wl,-O2

.PHONY: all bin test tool deploy ga clean

all: tool

test:
	$(CXX) preprocess/test.cc $(COMMON_SRC) -o $(BIN) $(GDAL) $(LIBS) $(HDF5) $(STD) && ./$(BIN)

tool:
	$(CXX) $(OPT_FLAGS) -fopenmp \
		cli/main.cc common/lineofsight.cc common/image.cc common/coverage.cc render.cc \
		common/angleviewshed.cc common/radarlos.cc common/mult.cc common/operations.cc \
		preprocess/jsonParser.cc preprocess/reproject.cc preprocess/preprocess.cc \
		preprocess/lib/*.cc ga/config.cc iad/*.cc \
		-o $(ADAPT) $(LIBS) $(STD)


ga:
	$(CXX) $(OPT_FLAGS) -pthread \
		ga/main.cc ga/config.cc $(COMMON_SRC) \
		-o $(BIN) $(GDAL) $(LIBS) $(HDF5) $(STD) && ./$(BIN)

init:
	mkdir  .\output\radar\ .\output\viewshed\ .\output\pso\ .\output\ga\ .\output\radarfd\   .\output\path\  .\output\simdata\   .\output\gbadws\  .\output\heatmap\

clean:
	rm -f $(BIN) $(ADAPT) *.o

dll:
	g++ -shared -O3 -ffast-math -flto -ffp-contract=fast -fopenmp  -march=native -funroll-loops -fprefetch-loop-arrays -fdata-sections -ffunction-sections -Wl,--gc-sections -Wl,-O2 common/lineofsight.cc common/image.cc common/coverage.cc common/angleviewshed.cc common/radarlos.cc common/mult.cc  common/operations.cc preprocess/jsonParser.cc preprocess/reproject.cc export-dll/export.cc preprocess/preprocess.cc render.cc common/radarfdr.cc preprocess/lib/*.cc ga/config.cc common/lib.cc iad/*.cc -o adapt.dll  -ljsoncpp -lgdal -lpng16 -lhdf5 -static-libgcc -static-libstdc++ -Wl,--kill-at -std=c++23
