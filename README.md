# Trace (Terrain & Radar Analysis & Coverage Engine)

**A C++23 high performance library for scalable, constant-memory terrain analysis — with radar siting and fleet optimization**

`trace` is, first, a library: a terrain engine that turns a digital elevation model and a handful of vector constraint layers (roads, railways, water, restricted boundaries) into a single tiled, queryable terrain store, plus a terrain-aware line-of-sight coverage model and two deployment optimizers built against it with flexible radar models and multiple configurations. It holds memory to the size of the area under evaluation rather than the size of the DEM it's built from — country-scale elevation data and a 150 km AOI query cost roughly the same either way (see [Benchmarks](#benchmarks)).

Everything in this repo that *runs* — the CLI, the export library — is a product built directly on top of that library, not a separate implementation of it:

- **Library core** — preprocessing (DEM + shapefiles → a tiled elevation/constraint store with one 8-bit feasibility bitmask per cell), the coverage model, and the GA/PSO optimizers below. The preprocessing stage alone has no radar-specific code in it — any terrain-constrained siting problem can be built directly on top of just that stage, without touching the rest — see [Architecture §9](docs/ARCHITECTURE.md#9-building-on-just-the-preprocessing-stage).
- **CLI (`trace`)** — a thin product built on the library: preprocessing, viewshed/line-of-sight, single fixed-site radar coverage, and both deployment optimizers below, as one binary. It's also how the library itself is exercised and benchmarked (see [Benchmarks](#benchmarks)). Full reference: [cli/README.md](cli/README.md).
- **Multi-radar deployment** — genetic-algorithm siting of several (optionally heterogeneous) radars at once, jointly optimizing placement and fleet size.
- **Single-radar deployment** — particle-swarm siting of one radar, with a benchmark mode that sweeps PSO strategy profiles and ranks them.
- **Export library** ([export/](export/)) — another product built on the library: the same engines as a flat C API, for embedding in an external host application instead of shelling out to the CLI.
- **Output formats** — GeoTIFF (georeferenced, EPSG:7755 native with optional EPSG:4326 warp) and PNG for every analysis output; SVG additionally for viewshed and single-site radar coverage.

## Architecture

Terrain discretization and tiled storage, the constraint-bitmask feasibility model, the curvature- and horizon-aware coverage model, and both optimizers are described in detail, mapped to their implementation, in [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md).

## Benchmarks

Full tables (hyperparameter leaderboard, terrain cell counts) in [docs/BENCHMARKS.md](docs/BENCHMARKS.md).

**Memory — full terrain store vs. a windowed AOI read**, at a 150 km AOI, across resolutions:

| DEM resolution | Full terrain store | Windowed AOI read | Reduction |
|---|---|---|---|
| 100 m | 7.8 GB | 17.9 MB | ~435× |
| 90 m | 9.6 GB | 22.1 MB | ~434× |
| 30 m | 86.5 GB | 198.7 MB | ~435× |

**Single-radar coverage-evaluation time**, milliseconds, by DEM resolution and maximum radar range:

| DEM resolution | 750 km | 500 km | 250 km | 150 km | 100 km | 50 km |
|---|---|---|---|---|---|---|
| 30 m | 512.1 | 162.2 | 78.8 | 59.3 | 52.4 | 42.2 |
| 90 m | 113.6 | 51.5 | 16.6 | 7.3 | 5.2 | 3.6 |
| 1000 m | 3.5 | 3.0 | 2.8 | 3.0 | 2.9 | 2.8 |

**GA search quality** across 16 selection/crossover/mutation operator profiles, same scenario and generation budget: coverage ranged 84.75%–95.98% of the AOI, with tournament selection, single-point crossover, and nudge mutation the strongest profile found.

**Point-to-point LOS evaluation time** — consistently sub-millisecond (150–200 µs) across tested distances up to 830 km on a 90 m DEM; full figures and the `trace rlos` command to reproduce them are in [docs/BENCHMARKS.md](docs/BENCHMARKS.md#point-to-point-line-of-sight-los-evaluation-time).

## Quick start

```sh
make tool        # Linux; see docs/BUILDING.md for macOS/Windows
cp config/radar.example.jsonc config/radar.json   # edit it, then:
./trace radar
```

## Docs

- [docs/GETTING_STARTED.md](docs/GETTING_STARTED.md) — install dependencies (Linux / macOS / Windows) and compile.
- [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) — architecture and algorithms, mapped to code.
- [docs/BENCHMARKS.md](docs/BENCHMARKS.md) — full benchmark tables and discussion.
- [docs/BUILDING.md](docs/BUILDING.md) — make-target reference, Makefile variants, troubleshooting.
- [cli/README.md](cli/README.md) — full CLI manual, one section per subcommand.
- [config/README.md](config/README.md) — what each config file does and how to write one.
