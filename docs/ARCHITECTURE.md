# Architecture

This document describes the engine's architecture — the terrain pipeline, the
coverage model, and the two deployment optimizers built on top of it — and maps
each piece to where it actually lives in this repository.

## 1. Overview

Radar deployment planning splits into two related questions: *where* should a
fixed number of radars go to maximize coverage, and *how many* radars are
needed to reach a coverage target at minimum cost. Both are posed here as the
same search problem — placement and (optionally) fleet size are both genes in
a chromosome optimized by a genetic algorithm — rather than as two separate
tools.

That search sits on top of a three-stage pipeline, and each stage is a
self-contained module that the next one consumes through a narrow interface:

1. **Terrain preprocessing** turns a DEM and a set of vector constraint layers
   into one aligned, tiled grid: an elevation raster and an 8-bit per-cell
   constraint bitmask.
2. **Feasible-grid extraction** intersects that grid with an Area of Interest
   (AOI) polygon to produce the candidate cell set a search is allowed to
   place a radar on.
3. **Coverage-aware search** (genetic algorithm for multiple radars, particle
   swarm optimization for one) evaluates candidate placements against a
   terrain-masked, curvature-corrected line-of-sight coverage model and
   iterates toward a placement that maximizes coverage (or coverage at
   minimum fleet cost).

Nothing above stage 1 is radar-specific. Section 9 expands on that.

## 2. Terrain preprocessing & constraint encoding

The AOI is discretized into an `R × C` grid at resolution `ρ` meters, with
cell `(r, c)` centered at `(x₀ + (c + 0.5)ρ, y₀ − (r + 0.5)ρ)` for grid origin
`(x₀, y₀)`. Every downstream stage — elevation lookup, constraint
masking, coverage evaluation — shares this exact resolution and alignment, so
layers are never resampled against each other at query time.

Constraint layers (roads, railways, water areas, water lines, and a
restricted/"IB" boundary) are supplied as vector shapefiles. Each is
rasterized onto the grid, then expanded by its own configured standoff
distance `d_k` using a two-pass row/column distance-propagation scan — forward
top-to-bottom/left-to-right, then backward bottom-to-top/right-to-left,
threshold at `d_k`. This runs in `O(RC)` time with `O(C)`
memory (only the current and previous row are resident), which is why it's
cheap enough to run per constraint type rather than needing a one-shot
brute-force nearest-constraint search. Implementation:
[preprocess/lib/matrix.cc](../preprocess/lib/matrix.cc).

A slope-exceedance flag is computed from the elevation grid using a
central-difference approximation of Horn's method.

All of the above collapses into a single 8-bit bitmask per cell — one bit per
constraint type (road / rail / water area / water line / restricted boundary /
slope, with bits reserved for future constraint types) — so a feasibility
query against *any* combination of constraints, at any later stage, is one
bitwise AND rather than a lookup across several matrices or a database layer.
Implementation: [preprocess/bin.h](../preprocess/bin.h) (`BIT_ROAD`,
`BIT_RAIL`, `BIT_WATER_AREA`, `BIT_WATER_LINE`, `BIT_IB`, `BIT_SLOPE`).

## 3. Tiled, windowed terrain storage

Cell count scales as `1/ρ²` — a DEM at 30 m resolution over a fixed extent
holds roughly 11× the cells of the same extent at 100 m. Holding the
elevation grid, slope grid, and constraint bitmask as monolithic `R × C`
arrays is feasible at coarse resolution; at country scale and fine resolution
it is not, because memory cost grows quadratically with `1/ρ` while the AOI
extent stays fixed.

The engine instead writes each grid to disk as a directory of fixed-size
`τ × τ` tiles plus one metadata file recording `(R, C, ρ, x₀, y₀, τ)`. Tile
files are named by their own tile-grid coordinates and extent (e.g.
`<prefix>_tr<tr>_tc<tc>_r<r0>_c<c0>_h<h>_w<w>.bin`), so opening a store is a
directory scan, not an index lookup. Implementation:
[preprocess/tile.h](../preprocess/tile.h).

At runtime, no stage loads a grid in full. A single rectangular window
`W = [r₀, r₀+H) × [c₀, c₀+W)` is read from each store once, sized to the AOI
polygon's bounding rectangle expanded by a margin equal to the largest sensing
range under consideration. Every subsequent lookup — feasibility
extraction, coverage evaluation — is indexed against that window, never
against the full grid, so peak memory for a given run is `O(|W|)`: governed by
the AOI's extent and the fleet's sensing range, not by how much terrain the
underlying store covers.

The effect is a close to two-orders-of-magnitude reduction in resident memory
for AOI-scoped runs over large terrain stores — see [Benchmarks](../README.md#benchmarks)
in the main README for measured figures across four resolutions.

## 4. Feasible grid extraction

Given the AOI polygon and the constraint bitmask, the coverage target set
`A` and the feasible siting set `F ⊆ A` are produced by: computing the
polygon's axis-aligned bounding rectangle in grid coordinates, testing each
cell center for polygon membership with a standard even-odd ray-casting test,
and applying the constraint bitmask filter — reject a cell from `F` if it has
water, restricted-boundary, or excessive-slope bits set, or (optionally) if it
has neither road nor rail proximity. Every cell of `A` is a
coverage target; only cells of `F` may host a radar. Implementation:
[`Bin::getFeasibleGrid`](../preprocess/lib/query.cc).

## 5. Radar coverage model

A candidate radar's nominal range is obtained from the radar range equation,

```
d_max = [ Pt · G² · λ² · σ / ((4π)³ · S_min) ] ^ (1/4)
```

This is used only when a run is configured to derive range from radar-model
parameters rather than take a configured range directly. What actually gets
filled in inside that nominal disk is a separate question, answered by
terrain-aware line-of-sight analysis rather than assumed.

For a candidate site, rays are cast outward from the source to every boundary
cell of the maximum-range circle using Bresenham's line algorithm, and the
union of all rays covers the interior of the disk. Along each
ray, apparent terrain elevation at ground distance `d` is earth-curvature
corrected,

```
h̃(d) = h_terrain(d) − d² / (2 · k_e · R_E)
```

(effective-earth-radius model, `k_e = 4/3`, optional), the terrain
slope seen from the source is tracked as a running maximum along the ray, and
a cell is classified `COVERED`, `TERRAIN_MASKED`, `BELOW_BEAM`/`ABOVE_BEAM`
(outside the vertical scan envelope), `BEYOND_HORIZON`, or `OUTSIDE_RANGE`/
`OUTSIDE_SECTOR`. Implementation:
[common/coverage.cc](../common/coverage.cc) (`circleBoundary`,
`castRadarRayFast`/`castRadarRayDiag`); status codes in
[operations.h](../operations.h).

For the GA and PSO fitness functions, which evaluate this model once per
candidate placement per generation, a second implementation reuses the same
ray-casting logic but allocates only a local `(2R+3)²` bounding-box buffer per
placement instead of a full-grid-sized buffer — a constant-factor memory
reduction per evaluation that matters because the search evaluates thousands
of candidates per run. Implementation:
[multi-radar-deployment/radar.h](../multi-radar-deployment/radar.h)
(`computeLocalMask`).

Point-to-point (single target) line-of-sight — used by the CLI rather than by
the optimizers — is a related but distinct computation: a fixed target slope
is known up front, so the terrain-masking test can terminate on the first
blocking cell instead of sweeping a full radial boundary. Two variants exist
for two different questions. `rlos` answers "can this radar, with its antenna
height/elevation-angle cone/horizon/curvature model, see this target" —
implementation: [common/radarlos.cc](../common/radarlos.cc)
(`RadarEngine::radarLOS`/`radarLOSDetailed`). `los` answers the plainer
"is there unobstructed line of sight between these two points, at these two
heights" with none of the radar-model gating — implementation:
[common/lineofsight.cc](../common/lineofsight.cc) (`Global::lineOfSight`).

## 6. Multi-radar deployment (genetic algorithm)

Two chromosome encodings are supported. **Fixed-length**, used when the radar
count `n` is specified up front: a chromosome is an ordered `n`-tuple of
feasible-site indices, one radar per gene — a maximal covering location
problem. **Variable-length** ("auto"), used for joint position-and-count
search: a chromosome is a variable-length list of (site, radar-type) genes,
free to vary between a configured minimum and maximum across the population
and across generations.

Two fitness modes score a fixed-length chromosome. *Weighted-sum* mode
scores unique coverage plus a bonus for matching a target overlap
fraction, minus a penalty on 3+-radar overlap. *Overlap-target* mode
instead rewards closeness to a configured redundancy target directly. The
variable-length encoding's fitness additionally subtracts a configured
per-radar cost penalty, trading coverage against fleet size.

Coverage overlap for a candidate solution — unique vs. multiply-covered cells
— is computed directly from each radar's masked line-of-sight visibility
rather than from intersecting nominal range circles, so the overlap accounting
stays consistent with the curvature-corrected, horizon-limited coverage model
used to score coverage in the first place.

Selection supports tournament, fitness-proportionate (roulette-wheel),
rank-based, and stochastic-universal-sampling. Crossover supports
single-point, multi-point, and uniform recombination (plus random-subset and
length-normalized recombination for the variable-length encoding). Mutation
supports single-point, spatial-nudge, and Gaussian-creep perturbation (plus
add/remove/relocate-radar operators for variable-length). Each generation
keeps the best chromosome unchanged (elitism) and the run terminates at a
configured generation limit or when the population's fitness variance drops
below a convergence threshold. Implementation:
[multi-radar-deployment/ga.h](../multi-radar-deployment/ga.h).

Heterogeneous fleets (a catalog of radar types with different range,
antenna height, and scan envelope) are supported directly — a chromosome gene
indexes into a shared type table rather than copying radar parameters per
gene.

## 7. Single-radar deployment (particle swarm optimization)

For the single-radar case, particle swarm optimization searches the same
feasible-grid space for the one site that maximizes unique AOI coverage.
Three independent strategy axes are configurable — topology (`global` /
`ring`), inertia schedule (`constant` / `linear_decay` / `jitter`), and
boundary handling (`clamp` / `bounce`) — and a benchmark mode sweeps every
combination of the three and reports a ranking, so a topology × inertia ×
boundary profile can be chosen empirically before being set for normal runs.
Implementation: [single-radar-deployment/pso.h](../single-radar-deployment/pso.h).

## 8. Outputs

Coverage, viewshed, and overlap rasters are written as GeoTIFF (georeferenced,
EPSG:7755 native with optional EPSG:4326 warp), PNG, or SVG, with an optional
single-band palette layout alongside the default 4-band RGBA. Implementation:
[render.cc](../render.cc), [common/image.cc](../common/image.cc).

Every engine above is also reachable without the CLI, through a flat C API
(`extern "C"`) intended for embedding in an external host application.
Implementation: [export/](../export/), documented in
[export/readme.md](../export/readme.md).

## 9. Building on just the preprocessing stage

Stages 2–4 above (tiled terrain storage, constraint-bitmask generation,
feasible-grid extraction) have no radar-specific code in them — they take a
DEM, a set of constraint shapefiles, and an AOI polygon, and produce a feasible
cell set and an elevation lookup, both addressable at `O(1)` per query and
both bounded in memory by the window under evaluation rather than the full
terrain store. Any siting or facility-location problem that needs
"terrain-aware feasible locations inside a boundary, fast, at country scale"
can be built against [preprocess/bin.h](../preprocess/bin.h)'s `Bin` and
`Global` types directly, without touching the GA, PSO, or coverage-model code
at all — a different fitness function over the same feasible-grid primitive
is enough to repurpose the engine for a different placement problem entirely.

## Case studies

This pipeline has been validated against four scenarios spanning both
planning echelons it's designed for: a fixed-count homogeneous fleet, a
variable-count homogeneous fleet across three AOIs of differing terrain
(returning 5–7 radars depending on terrain), and a variable-count
heterogeneous fleet mixing two radar types. Across all four, terrain
masking — not nominal instrumented range — was found to be the binding
constraint on achieved coverage.
