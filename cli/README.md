# `trace` CLI manual

`trace` (built from [cli/main.cc](main.cc)) is the command-line entry point for the
toolkit: terrain preprocessing, viewshed/LOS analysis, single-radar coverage, and
the multi-radar (GA) / single-radar (PSO) deployment optimizers.

Build it first — see [../docs/BUILDING.md](../docs/BUILDING.md). All examples below
assume you're running `./trace ...` from the repo root, and that the `output/`
directory tree already exists (`make clean` creates it).

Run `./trace --help` or `./trace <subcommand> --help` any time for the authoritative,
in-binary option list — this file explains what each one does and which config file
(if any) it reads. See [../config/README.md](../config/README.md) for the config files
themselves.

```
./trace [OPTIONS] SUBCOMMAND
```

Subcommands: `build`, `viewshed`, `radar`, `rlos`, `los`, `run` (with `run single`
and `run single benchmark`).

---

## `build` — preprocess maps into `.bin` files

```
./trace build [-a | --road | --rail | --water-areas | --water-lines | -i,--ib | --slope] [-v] [-c]
```

Reads **`config/config.json`** (paths to the DEM + shapefiles, the grid definition,
and feasibility distances — required), warps/rasterizes the inputs, and writes the
`.bin` files listed under `config.json`'s `output` block.

| Flag | Effect |
|---|---|
| *(none)* | Build every layer (elevation + all enabled constraints) into the combined map. |
| `-a, --all` | Also write each constraint layer out as its own standalone `.bin`, in addition to the combined map. |
| `--road` / `--rail` / `--water-areas` / `--water-lines` / `-i, --ib` / `--slope` | Build **only** that one constraint layer (mutually exclusive — pick at most one). |
| `-v, --with-images` | Also render a PNG preview of each bin written. |
| `-c, --count-time` | Print how long elevation vs. constraint-map building took. |

```
./trace build -a -v        # full build, per-layer bins, with preview PNGs
./trace build --road       # just the road buffer layer
```

---

## `viewshed` — line-of-sight visibility from one point

```
./trace viewshed --lat LAT --lon LON [options] (--tif | --svg | --png)
```

No config file — everything comes from flags. `--lat`/`--lon` are required, and
exactly one of `--tif` / `--svg` / `--png` must be given.

| Flag | Default | Meaning |
|---|---|---|
| `--lat`, `--lon` | *(required)* | Observer location (WGS84). |
| `-r, --range` | `300` | Range limit, km. |
| `-a, --antenna` | `30` | Antenna/observer height, m AGL. |
| `-s` | off | Restrict to a sector instead of full 360°. |
| `--loa`, `--roa` | `0`, `90` | Sector left/right-of-axis bounds, degrees (only with `-s`). |
| `--use-curvature` | off | Include earth-curvature correction in the LOS calculation. |
| `--single-band` | off | Single-band palette TIFF instead of multi-band RGBA. |
| `--epsg-7755` | off | Reproject TIFF output to EPSG:7755 instead of WGS84. |
| `-o, --output` | `output/viewshed/viewshed_360.<ext>` or `..._sector.<ext>` | Output path; extension is forced to match `--tif/--svg/--png`. |

```
./trace viewshed --lat 33.9871 --lon 74.7742 -r 150 --png
./trace viewshed --lat 33.9871 --lon 74.7742 -s --loa -45 --roa 45 --tif
```

---

## `radar` — single fixed-position radar coverage

```
./trace radar
```

Takes **no flags**. Reads **`config/radar.json`** entirely — position, radar
parameters, curvature/AMSL/radar-model switches, and the output path (its
extension picks PNG/SVG/TIFF). Fails with `[ERROR] radar conf not found` if the
file is missing. To change anything, edit `config/radar.json` and re-run.

```
./trace radar
```

---

## `rlos` — radar line-of-sight between two points (fast path check)

```
./trace rlos --src-lat LAT --src-lon LON --tgt-lat LAT --tgt-lon LON [options]
```

| Flag | Default | Meaning |
|---|---|---|
| `--src-lat/lon`, `--tgt-lat/lon` | *(required)* | Radar and target positions. |
| `-a, --antenna` | `30` | Radar antenna height, m AGL. |
| `-t, --target` | `1000` | Target height, m AGL. |
| `--min-elev`, `--max-elev` | `0.5`, `30` | Radar elevation-angle limits, degrees. |
| `-r, --range` | `300` | Max calculation range, km. |
| `-d, --detailed` | off | Print per-cell telemetry along the ray instead of just the summary. |

```
./trace rlos --src-lat 34.0173 --src-lon 74.5059 --tgt-lat 34.1691 --tgt-lon 74.7127 -d
```

## `los` — plain line-of-sight between two points

Same source/target/antenna/target flags as `rlos` minus `--min-elev`/`--max-elev`/
`--range`/`-d` (no radar-model gating, no detailed mode — just "is there
unobstructed line of sight between these two points"). Prints source/target
elevation, ray cell count, distance, and compute time, then `Target Status :
VISIBLE`/`BLOCKED`.

```
./trace los --src-lat 34.0173 --src-lon 74.5059 --tgt-lat 34.1691 --tgt-lon 74.7127
```

---

## `run` — multi-radar deployment (genetic algorithm)

```
./trace run [-g N] [--reduce PCT] [--png] [-b] [--single-band] [--epsg-7755]
```

Reads **`config/gaconf.json`**. Exactly one of its `auto` / `heterogenous` /
`homogenous` blocks should have `"case": true` — that's the scenario that runs.
Does nothing (prints nothing, exit 0) if `gaconf.json` is missing.

| Flag | Default | Meaning |
|---|---|---|
| `-g, --gen-simdata` | `0` | If non-zero, write simulation data for the top-`N` solutions instead of running the full optimizer. |
| `--reduce` | `50` | Output raster resolution reduction, percent. |
| `--png` | off | Also write a PNG alongside the GeoTIFF. |
| `-b, --run-benchmark` | off | Run the benchmark sweep (every selection/mutation/crossover combination) for whichever case is active, instead of a single optimization run. |
| `--single-band` | off | Single-band palette TIFF instead of multi-band RGBA. |
| `--epsg-7755` | off | Reproject TIFF output to EPSG:7755 instead of WGS84. |

Output goes to the path in `gaconf.json`'s `"output"` field (default
`output/ga/ga.tif`).

```
./trace run --reduce 25 --png
./trace run -b                 # benchmark the active case's operators
./trace run -g 10               # dump simulation data for the top 10 solutions
```

### `run single` — single-radar siting (PSO)

```
./trace run single
./trace run single benchmark
```

Reads **`config/psoconf.json`** (fixed radar spec; PSO hyperparameters, topology/
inertia/boundary strategy, output path). `run single` finds the one best site for
a single radar. `run single benchmark` instead sweeps every
topology × inertia × boundary combination and prints a ranking — use it to pick a
strategy, then set that strategy in `psoconf.json` for normal runs.

```
./trace run single
./trace run single benchmark
```
