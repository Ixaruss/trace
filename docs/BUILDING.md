# Building & running locally

New here? [GETTING_STARTED.md](GETTING_STARTED.md) has per-distro/OS dependency install steps and a first compile. This doc covers make-target reference, the Makefile variants, and troubleshooting.

## Dependencies

- A C++23 compiler — GCC 13+ or Clang 16+ (this repo is developed against GCC 16).
- [GDAL](https://gdal.org/) (raster/vector I/O, reprojection)
- [jsoncpp](https://github.com/open-source-parsers/jsoncpp) (config files)
- libpng16 (PNG output)
- HDF5 (used by the GA engine)
- OpenMP (ships with GCC; macOS/Clang needs `libomp` separately, see below)
- GNU Make

CLI11 is vendored as a single header ([cli/cli.h](../cli/cli.h)) — nothing to
install for the CLI argument parsing itself.

There are three Makefiles in the repo:

- **[Makefile](../Makefile)** (repo root) — the one that's actually kept in sync
  with the current source tree. Use this on Linux, and on Windows under MSYS2.
- **[make/makefile_mac](../make/makefile_mac)** — macOS-specific: uses
  `clang++`, Homebrew library paths, and `-Xpreprocessor -fopenmp` + `libomp`
  for OpenMP (Apple Clang has no built-in OpenMP support). Use this on macOS.
- **[make/makefile_linux](../make/makefile_linux)** — an older snapshot of the
  root Makefile, kept for reference. Prefer the root `Makefile` on Linux.

All three expect to run from the repo root and write into `output/<something>/...`
— create that tree once before running anything:

```sh
mkdir -p output/{ga,gbadws,heatmap,path,preprocess,pso,radar,simdata,viewshed}
```

(`make clean` also (re)creates this tree, but it `rm -rf`s `output/` and the
built binaries first — fine for a first-time setup, destructive on a repo with
output you want to keep.)

You'll also need config files under `config/` before running most subcommands —
see [config/README.md](../config/README.md).

## Linux

```sh
# Debian/Ubuntu
sudo apt install build-essential libgdal-dev libjsoncpp-dev libpng-dev libhdf5-dev

make tool        # builds ./trace
./trace --help
```

Other targets: `make test` (sanity-check preprocessing), `make dll` (builds a
Windows DLL). `make dll` hardcodes plain `g++`, so it won't link on a native
Linux toolchain — `export/export.cc` uses Windows-only headers
(`<debugapi.h>`). Run it under MSYS2 (see **Windows** below), or substitute a
mingw-w64 cross-compiler by running its command line directly — see
[export/readme.md](../export/readme.md).

## macOS

```sh
brew install gdal jsoncpp libpng hdf5 libomp

make -f make/makefile_mac tool   # builds ./trace
./trace --help
```

`make/makefile_mac` hardcodes Homebrew's **Apple Silicon** prefix
(`/opt/homebrew`). On an Intel Mac, Homebrew installs under `/usr/local` instead
— either symlink `/opt/homebrew` to your actual prefix, or edit the
`HOMEBREW_INCLUDES` / `HOMEBREW_LIBS` / `OMP_FLAGS` / `OMP_LIBS` / `HDF5_LIBS`
variables at the top of that file to match `brew --prefix`.

The `multi-radar-deployment` target in this Makefile builds just the GA engine
standalone (`x_ga`); the `single-radar-deployment` target is a no-op reminder —
PSO is header-only now and runs through `./trace run single` after `make tool`,
not as its own binary.

## Windows

The CLI (and the rest of the Makefile targets) are POSIX-flavored — the
practical options are:

**MSYS2 (recommended for a native build):**
1. Install [MSYS2](https://www.msys2.org/), then in an MSYS2 **MinGW64** shell:
   ```sh
   pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-gdal \
             mingw-w64-x86_64-jsoncpp mingw-w64-x86_64-libpng \
             mingw-w64-x86_64-hdf5 make
   ```
2. From the repo root in that same shell: `make tool` → `./trace.exe --help`.
   (The `tool` target doesn't use `gdal-config` — it links against whatever
   `-lgdal`/`-ljsoncpp`/etc. resolve to on the default include/lib search path,
   which the `mingw64` packages above populate automatically.)

**WSL (simplest):** install a WSL2 Linux distro and follow the **Linux**
instructions above inside it.

**`build.bat`** at the repo root is a separate, older script: it builds a
different artifact (a `maplib.dll` from `preprocess/export.cc`/`export.h`, not
the `trace` CLI or `export/export.cc`) and its source file list predates
some of the moves in this tree. Treat it as a legacy reference, not a verified
build path — [export/readme.md](../export/readme.md) has the current
command line for the `trace.dll` export build (`make dll`), which still needs a
MinGW-flavored toolchain (MSYS2, or cross-compiling from Linux) because it links
against Windows-only headers.

## Running

```sh
./trace build       # preprocess maps into .bin files (needs config/config.json)
./trace viewshed --lat 33.9871 --lon 74.7742 --png
./trace radar                      # needs config/radar.json
./trace run                        # multi-radar GA, needs config/gaconf.json
./trace run single                 # single-radar PSO, needs config/psoconf.json
```

See [cli/README.md](../cli/README.md) for the full command reference.
