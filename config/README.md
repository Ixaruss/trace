# Config files

`trace` reads its config as plain JSON from **fixed, relative paths** — always
relative to the current working directory, so run `trace` from the repo root (or
make sure a `config/` folder sits next to wherever you invoke it from). There's no
`--config` flag; the path is hardcoded per subcommand.

| File (must be named exactly this) | Read by | Required for |
|---|---|---|
| `config/config.json` | `Config::load()` | `trace build` |
| `config/gaconf.json` | `GaConf::load()` | `trace run` |
| `config/psoconf.json` | `PsoConf::load()` | `trace run single` |
| `config/radar.json` | `Rconf::load()` | `trace radar` |

None of these ship in the repo (only example templates do, see below) — each
subcommand checks `is_available()` first and fails cleanly (or silently no-ops, in
`run`'s case) if its file is missing.

## Getting started: copy the examples

Every file above has a matching `*.example.jsonc` template in this folder, fully
annotated inline (field type, units, defaults, and gotchas). Copy the one(s) you
need and edit in place — keep the `.json` extension (not `.jsonc`) since that's the
exact filename the loaders look for:

```sh
cp config/config.example.jsonc  config/config.json
cp config/gaconf.example.jsonc  config/gaconf.json
cp config/psoconf.example.jsonc config/psoconf.json
cp config/radar.example.jsonc   config/radar.json
```

You only need the files for the subcommands you actually run — e.g. `trace radar`
only needs `radar.json`.

**Comments are fine to keep.** The JSON parser used throughout this project
(`Json::CharReaderBuilder`) allows `//` and `/* */` comments by default, so you
don't have to strip them from the files you copy — only strip them if you're
feeding the same file to some other, stricter JSON tool.

## What's in each file

Full field-by-field docs live as inline comments in the example files — read
those alongside this summary rather than a duplicated reference here, since
they're what actually gets kept up to date when a field changes.

- **`config.example.jsonc`** — map build input (DEM + shapefile paths), the AOI
  grid definition (rows/cols/resolution/extent in EPSG:7755), which constraint
  layers are enforced (`case`) and their buffer distances (`distance`), and where
  the resulting `.bin` files get written (`output`). `grid.rows/cols/min_*/max_*`
  are normally filled in automatically after a build, not hand-written.

- **`gaconf.example.jsonc`** — the multi-radar GA. Pick exactly one scenario by
  setting `"case": true` on exactly one of `auto` (GA chooses radar count/types
  within a min/max), `heterogenous` (fixed mix of radar types), or `homogenous`
  (N identical radars) — leaving more than one `true` means whichever the code
  checks first wins silently, and leaving all `false` makes `trace run` do
  nothing. Also holds GA hyperparameters (population/generations/crossover/
  mutation/early-stop), fitness weights, and model switches (earth curvature,
  radar-equation model, AMSL, threat sector).

- **`psoconf.example.jsonc`** — the single-radar PSO siting (`trace run single`).
  One radar spec, the model switches, output path, and PSO-specific knobs:
  `pso_hyperparameters` (c1/c2/inertia/velocity clamp), `pso_global_settings`
  (swarm size/iterations/seed), and `pso_strategy` (topology/inertia-schedule/
  boundary-handling — run `trace run single benchmark` to see which combination
  scores best before locking one in here).

- **`radar.example.jsonc`** — one radar at a fixed lat/lon (`trace radar`). Same
  radar-spec fields as the other two, plus the source position and an optional
  threat sector (`sector_flag`/`loa`/`roa`).

## Output paths

Each file's `output` field is a path under `output/` (e.g. `output/ga/ga.tif`).
Those directories aren't created automatically at runtime — make sure the
`output/` tree exists first (see [../docs/BUILDING.md](../docs/BUILDING.md)) or
the write will fail.
