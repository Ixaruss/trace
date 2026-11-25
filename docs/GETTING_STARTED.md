# Getting started

Dependency installation and a first compile, across Linux, macOS, and Windows. For make-target reference, troubleshooting, and the Windows DLL build, see [BUILDING.md](BUILDING.md).

## Dependencies

| Dependency | What it's for |
|---|---|
| A C++23 compiler | GCC 13+ or Clang 16+ (this repo is developed against GCC 16) |
| [GDAL](https://gdal.org/) | Raster/vector I/O, reprojection |
| [jsoncpp](https://github.com/open-source-parsers/jsoncpp) | Config file parsing |
| libpng16 | PNG output |
| HDF5 | Used by the GA engine |
| OpenMP | Parallelism — ships with GCC; macOS/Clang needs `libomp` separately |
| GNU Make | Running the build |

CLI11 (the CLI argument parser) is vendored as a single header ([cli/cli.h](../cli/cli.h)) — nothing to install for it.

## Linux

Pick your distribution's package manager below. Package names are current as of recent releases; if one's missing, search your distro's package index for the library name (e.g. "jsoncpp").

**Debian / Ubuntu** (`apt`):
```sh
sudo apt update
sudo apt install build-essential libgdal-dev libjsoncpp-dev libpng-dev libhdf5-dev
```

**Fedora / RHEL / CentOS Stream** (`dnf`):
```sh
sudo dnf install gcc-c++ make gdal-devel jsoncpp-devel libpng-devel hdf5-devel
```

**Arch Linux / Manjaro** (`pacman`):
```sh
sudo pacman -S base-devel gdal jsoncpp libpng hdf5
```

**openSUSE** (`zypper`):
```sh
sudo zypper install gcc-c++ make libgdal-devel libjsoncpp-devel libpng16-devel hdf5-devel
```

**Alpine** (`apk`):
```sh
sudo apk add build-base gdal-dev jsoncpp-dev libpng-dev hdf5-dev
```

Then build:

```sh
make tool
./trace --help
```

## macOS

Install [Homebrew](https://brew.sh/) if you don't already have it, then:

```sh
brew install gdal jsoncpp libpng hdf5 libomp
```

Build with the macOS-specific Makefile (handles Apple Clang's lack of built-in OpenMP and Homebrew's library layout):

```sh
make -f make/makefile_mac tool
./trace --help
```

If you're on an Intel Mac: `make/makefile_mac` hardcodes Homebrew's Apple Silicon prefix (`/opt/homebrew`). Either symlink `/opt/homebrew` to your actual prefix (`/usr/local` on Intel), or edit the `HOMEBREW_INCLUDES`/`HOMEBREW_LIBS`/`OMP_FLAGS`/`OMP_LIBS`/`HDF5_LIBS` variables at the top of that file to match `brew --prefix`.

## Windows

Install [MSYS2](https://www.msys2.org/), then open an **MSYS2 MinGW64** shell and run:

```sh
pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-gdal \
          mingw-w64-x86_64-jsoncpp mingw-w64-x86_64-libpng \
          mingw-w64-x86_64-hdf5 make
```

Then, from the repo root in that same shell:

```sh
make tool
./trace.exe --help
```

(Prefer not to install MSYS2? A WSL2 Ubuntu install and the Debian/Ubuntu steps above work too — simplest if you're already comfortable with WSL.)

## Next steps

- [cli/README.md](../cli/README.md) — full CLI manual
- [config/README.md](../config/README.md) — writing config files
- [ARCHITECTURE.md](ARCHITECTURE.md) — how the engine works
