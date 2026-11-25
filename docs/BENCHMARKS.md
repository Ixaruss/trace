# Benchmarks

Figures below characterize the terrain engine's memory scaling and
single-site coverage-evaluation cost, and the genetic algorithm's search
quality across operator configurations — the three properties the rest of
this system's performance is built on.

## Terrain cell count vs. resolution

Grid size for an Indian-subcontinent DEM extent used in the reference case
studies, at each evaluated resolution:

| Resolution `ρ` | Rows `R` | Cols `C` | Cells (`R × C`) |
|---|---|---|---|
| 1000 m | 4,450 | 6,266 | ≈2.79 × 10⁷ |
| 100 m | 44,494 | 62,659 | ≈2.79 × 10⁹ |
| 90 m | 49,438 | 69,621 | ≈3.44 × 10⁹ |
| 30 m | 148,312 | 208,861 | ≈3.10 × 10¹⁰ |
| 5 m | 889,871 | 1,253,162 | ≈1.12 × 10¹² |

Cell count scales as `1/ρ²`: moving from 100 m to 30 m resolution over the
same extent is an ~11× increase in cells, and 100 m to 5 m is a ~400× increase
— the motivation for the tiled, windowed storage described in
[Architecture §3](ARCHITECTURE.md#3-tiled-windowed-terrain-storage).

## Resident memory — full terrain store vs. windowed AOI read

Measured for a 150 km × 150 km AOI window against the full terrain extent
above:

| DEM resolution | Full terrain store | Windowed AOI read | Reduction |
|---|---|---|---|
| 100 m | 7.8 GB | 17.9 MB | ~435× |
| 90 m | 9.6 GB | 22.1 MB | ~434× |
| 30 m | 86.5 GB | 198.7 MB | ~435× |

The windowed figure is governed by the AOI's extent and sensing-range margin,
not by the resolution-driven growth of the underlying store — which is why
the reduction factor stays close to constant across resolutions even as both
absolute numbers grow by orders of magnitude.

## Single-radar coverage-evaluation time

Wall-clock time, in milliseconds, for one full radial line-of-sight coverage
sweep at a single candidate site, by DEM resolution and maximum configured
radar range:

| DEM resolution | 750 km | 500 km | 250 km | 150 km | 100 km | 50 km |
|---|---|---|---|---|---|---|
| 30 m | 512.143 | 162.187 | 78.832 | 59.323 | 52.420 | 42.150 |
| 90 m | 113.594 | 51.536 | 16.567 | 7.281 | 5.169 | 3.626 |
| 1000 m | 3.468 | 3.006 | 2.791 | 2.962 | 2.897 | 2.782 |

Cost scales with ray length (range) and with the number of rays cast per unit
of boundary arc (finer resolution → a longer circle boundary at the same
range), which is why the 30 m column dominates the table while the 1000 m row
stays roughly flat across range — at that resolution the ray count is small
enough that fixed per-call overhead, not ray-casting work, dominates.

This is the per-candidate cost paid by both optimizers' fitness function,
thousands of times per run — see
[Architecture §5](ARCHITECTURE.md#5-radar-coverage-model) for the
constant-factor memory optimization applied specifically to that hot path.

## Point-to-point line-of-sight (LOS) evaluation time

Distinct from the radial coverage sweep above: `rlos` checks visibility along
one ray between a fixed source and a fixed target, not every ray around a
boundary circle, and terminates as soon as the ray is found blocked rather
than sweeping to the configured range regardless. It's a different operation
with a different cost shape, not a cheaper version of the same one, so the
two tables aren't directly comparable entry-for-entry.

Measured on a real 90 m-resolution DEM tile set (the only resolution
available in the environment this was run in — unlike the multi-resolution
grid in the table above, this is a single data point per distance, not a
systematic sweep, and should be read as "LOS is consistently sub-millisecond
at this resolution," not as a precise resolution/distance surface):

| Source-target distance | Outcome | Compute time |
|---|---|---|
| 55.4 km | VISIBLE | 160.6 µs |
| 110.9 km | VISIBLE | 160.8 µs |
| 166.3 km | VISIBLE | 151.2 µs |
| 277.3 km | VISIBLE | 200.3 µs |
| 554.5 km | VISIBLE | 185.7 µs |
| 831.6 km | TERRAIN_MASKED | 149.7 µs |

Two things to read correctly here. First, `rlos` terminates at the first
blocking cell — the 831.6 km run didn't walk the full ray, it stopped
partway at whatever cell first broke line of sight, which is why its time
isn't the largest in the table despite having the longest distance; a
`VISIBLE` result at a given distance is a true full-ray-length measurement,
a `TERRAIN_MASKED` one is a lower bound. Second, these runs used an
artificially high target altitude (`-t 40000`) and disabled elevation-angle
gating (`--min-elev -90 --max-elev 90`) specifically to stop the default
~153 km horizon limit (at realistic antenna/target heights) and the
elevation-angle checks from short-circuiting the longer-distance rows before
they ever reached the terrain-masking walk — those are benchmark parameters
to force the ray-casting code path to run, not realistic operating values.

To generate a data point yourself, run the CLI's `rlos` subcommand (the
binary this repo builds is `trace`, via `make tool` — see
[Quick start](../README.md#quick-start)):

```sh
./trace rlos \
  --src-lat 34.017346 --src-lon 74.50587999 \
  --tgt-lat 34.169100 --tgt-lon 74.712671 \
  -a 30 -t 40000 --min-elev -90 --max-elev 90 \
  -r 800
```

`-r`/`--range` (default `300` km) is a hard cutoff checked *before* the ray
is cast — set it at or above the actual source-target distance, or the run
reports `BEYOND_RANGE` near-instantly without walking the ray at all.

Each run prints its own measured distance and timing directly — no separate
distance calculation needed:

```
Distance    : <ray length, km>
...
Compute     : <elapsed, µs>
```

To populate the grid: pick `--tgt-lat`/`--tgt-lon` pairs whose printed
`Distance` lands at (or close to) each column's target distance for a given
source, and reload `config/config.json` + rebuild `.bin`s at each row's DEM
resolution (`./trace build`) before running the column sweep at that
resolution. `-a`/`-t`/`--min-elev`/`--max-elev` only gate the result
classification, not the ray-casting cost itself, so they can be held fixed
across the whole grid.

## Genetic algorithm — hyperparameter search quality

16 selection/crossover/mutation operator configurations, evaluated against
the same scenario, population size (150), and generation budget (1000). Full
ranking:

| Rank | Selection | Crossover | Mutation | Coverage |
|---|---|---|---|---|
| 1 | Tournament | Single-Point | Nudge | 95.98% |
| 2 | Roulette | Multi-Point | Nudge | 93.42% |
| 3 | Tournament | Multi-Point | Nudge | 91.21% |
| 4 | SUS | Uniform | Gaussian | 91.14% |
| 5 | Rank-Based | Single-Point | Nudge | 90.99% |
| 6 | SUS | Uniform | Gaussian | 90.95% |
| 7 | SUS | Single-Point | Single | 90.80% |
| 8 | Tournament | Single-Point | Single | 90.69% |
| 9 | Tournament | Multi-Point | Nudge | 88.86% |
| 10 | Rank-Based | Uniform | Gaussian | 88.54% |
| 11 | SUS | Single-Point | Nudge | 88.32% |
| 12 | SUS | Multi-Point | Gaussian | 88.17% |
| 13 | Roulette | Uniform | Gaussian | 86.39% |
| 14 | Tournament | Single-Point | Nudge | 85.14% |
| 15 | Rank-Based | Single-Point | Single | 85.11% |
| 16 | Rank-Based | Multi-Point | Gaussian | 84.75% |

The 11-point spread (84.75%–95.98%) across operator choices, at a fixed
population size and generation budget, indicates the search is sensitive to
operator selection and not just to how long it's allowed to run — which is
why [`trace run single benchmark`](../cli/README.md) exists for the PSO side,
and why the GA side exposes every operator independently rather than hardcoding
one configuration.

Default configuration used for the reference case studies, also a reasonable
starting point for new scenarios:

| Parameter | Value |
|---|---|
| Population size | 150 |
| Max generations | 1000 |
| Selection | Tournament |
| Crossover | Single-point |
| Mutation | Gaussian |
| Crossover probability | 0.90 |
| Mutation probability | 0.10 |
| Target overlap `τ` | 0.25 |
| Count penalty `λ` | 0.20 |
