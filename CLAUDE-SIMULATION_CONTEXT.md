# SIMULATION_CONTEXT.md — k-Wave, Tissue Models & Sensor Design

> Referenced when working on: `run_single_field_simulation.m`, `run_standalone_simulation.m`, `create_acoustic_medium.m`, `determine_sensor_mask.m`, `apply_element_averaging.m`, `find_optimal_kwave_size.m`, `run_medium_comparison.m`, `test_time_dependence.m`, and the `+ethos/` class layer.

## `+ethos/` Class Layer (condensed simulation/analysis architecture)

An OOP layer that condenses the repeated CONFIG blocks + orchestration copied
across the `run_*` / `study_*` / `test_*` driver scripts. **It wraps the existing
functions; no physics moved into the classes.** Value classes are `parfor`-safe;
`Simulator`/`StudyRunner` are handles. All live in package folder `+ethos/`, so
each is one `.m` file called as `ethos.ClassName`.

### Classes

| Class | Kind | Role | Wraps / calls |
|---|---|---|---|
| `ethos.SimConfig` | value | Config object replacing the copy-pasted CONFIG blocks. Presets `standalone`/`pipeline`/`convergenceSweep`/`gammaBatch`; `withOverrides` (returns a copy); `validate`; `hash`; `toStruct` (plain struct bridge); static `tissueTables`. `dose_source` + `gamma_default=[3 3]`. | `compute_sim_config_hash` |
| `ethos.DoseSource` | enum | `Simulate` \| `Load` \| `Auto`. Drives `Simulator.resolve`. | — |
| `ethos.SimCase` | value | One immutable dose+geometry+medium+sensor bundle. Factories `fromDoseFile`, `fromProcessed` (**folds in `load_processed_data`**), `list` (lightweight, filtered). `materialize` builds medium+sensor on demand; static `buildSensor`. | `load_processed_data`, `list_processed_field_doses`, `load_field_dose_file`, `create_acoustic_medium`, `determine_sensor_mask` |
| `ethos.Simulator` | handle | Owns the load-or-simulate decision. `resolve` (per `dose_source`), `resolveMany` (parfor over value cases), `run` (force sim, optional recon cache), `runBlind` (forward on true CT, invert on reference CT), `loadField`/`loadTotal`. Always returns a `SimResult`. | `run_single_field_simulation`, `run_second_field_simulation`, `load_recon_dose_data` |
| `ethos.SimResult` | value | **Universal currency** — identical whether simulated or loaded. `recon_dose`, `rs_truth`, `ethos_truth`, `density`, `sensor_mask/info`, `metadata`, `rtplan`, `p0`/`reconPressure` (sim-only), `provenance`. Adapters `fromField`/`fromTotal`/`fromSim`. | `load_recon_dose_data` output |
| `ethos.Analysis` | static | `gamma` (default **3%/3 mm**, global, `limit=2·DTA`, `restrict=1`), `ssim` (scored over the **same 10%-of-truth mask** as gamma), `gammaSweep` (n%/n mm curve). | `CalcGamma`, `ssim` |
| `ethos.StudyRunner` | handle | Source-agnostic studies. Metrics `gamma`/`ssim`/`gammaChange`/`gammaSweep`; studies `compare`, `compareBlind` (reference full recon vs blind recon of the adapted anatomy), `sweep`, `paramSweep` (single-field axis sweep keeping every recon; `ScoreVs` truth/first/reference-run), `gammaBatch` (A1/A2 per dose), `beamSummary` (per-beam gamma/SSIM over all segments). `ensureSensor`; plot helpers `plotResult`/`plotSweep`/`plotParamSweep`/`plotBeamSummary`/`plotPassRateCurve`. | `ethos.Analysis`, `ethos.Simulator`, `ethos.SimPlotter` |
| `ethos.SimPlotter` | static | Single home for the figure primitives (was duplicated, file-local, in every `run_*`). `reconVsTruth`, `changePanels`, `gammaMap`, `sensorPlacement`, `doseDifference`, `convergence`, `all`. Ports `plot_truth_recon_diff_axial` / `plot_sensor_dose_planes` / gamma-map panel. | — |

### Key invariants

- **`SimResult` is the only currency the analysis/plot layer touches**, so
  "some studies simulate, some just analyze" is one flag (`dose_source`), not two
  code paths. `Auto` loads a hash-matching recon if present, else simulates.
- **`toStruct()` is the bridge**: every physics function keeps taking a plain
  `config` struct unchanged. Adopt the classes incrementally driver-by-driver.
- **Hash parity**: `SimConfig.toStruct` emits the same field names a pipeline
  CONFIG carries, so a `SimConfig` hashes identically to the CONFIG that wrote
  the on-disk recons (`load_recon_dose_data` also disk-scans as a fallback).
- **Sensor on loaded results**: `load_recon_dose_data` returns no sensor mask;
  `StudyRunner.ensureSensor(result, simCase)` rebuilds it via
  `SimCase.buildSensor` for sensor-placement plots.
- **Not ported**: the redraw-noise ensemble null stays in
  `study_pass_rates_individual.m` (it depends on forward-bundle internals that
  are file-local, not standalone-callable).

### Minimal usage

```matlab
cfg = ethos.SimConfig.standalone().withOverrides( ...
        'working_dir', WD, 'patient_id', '1194203', 'session', 'Session_1');
sr  = ethos.StudyRunner(cfg);
c   = ethos.SimCase.list('1194203','Session_1',cfg,'Beam',15,'Segment',112,'CTLabel','CT_1');
r   = sr.resolve(c(1));            % load-or-simulate -> SimResult
G   = sr.gamma(r);                 % 3%/3 mm vs RayStation truth
S   = sr.ssim(r);                  % SSIM over the gamma mask
ethos.SimPlotter.reconVsTruth(r, G);   ethos.SimPlotter.sensorPlacement(r);
```

The driver scripts have been migrated to this layer where they map cleanly (see
the migration table below). `run_standalone_simulation.m` is the canonical thin
example: preset + a few `withOverrides`, then `resolve` / `gamma` / `ssim` /
`SimPlotter`.

### Driver migration status

Each kept script carries a `[+ethos]` header note explaining the same rationale.

| Script | Status | Class path / why kept |
|---|---|---|
| `run_standalone_simulation.m` | **Migrated** | `Simulator.run` + `SimPlotter` (thin) |
| `run_standalone_comparison.m` | **Migrated** | `StudyRunner.compareBlind` |
| `run_standalone_analysis.m` | **Migrated** | `gamma`/`gammaChange`/`gammaSweep` (load-only) |
| `run_dt_comparison.m` | **Migrated** | `paramSweep('cfl_number', …, 'ScoreVs','first')` (dt ∝ cfl) |
| `run_air_sound_speed_comparison.m` | **Migrated** | `paramSweep` + `SimConfig.withTissue` (air row) |
| `study_pass_rates_allsegments.m` | **Migrated** | `StudyRunner.beamSummary` |
| `study_pass_rates_individual.m` | **Partial** | A1/A2 → `gammaBatch`; kept for the noise-ensemble null (forward-bundle internals) |
| `run_standalone_simulation_upgraded.m` | Kept | Experimental superset of the standalone run |
| `test_dx_comparison.m` | Kept | downscale changes grid → cross-grid gamma (not in `Analysis.gamma`) |
| `study_optimization_sweeps.m` | Kept | cross-grid downscale + blind + reference-run scoring |
| `run_nt_convergence_sweep.m` | Kept | Nt divisor is not a `run_single_field_simulation` param |
| `study_gamma_index_convergence.m` | Kept | needs per-TR-iteration recon |
| `run_gamma_convergence_batch.m` | Kept | thin orchestrator of the above study |
| `test_time_dependence.m` | Kept | bespoke inline k-Wave time-dependence probe |
| `run_medium_comparison.m` | Kept | synthetic phantom media + custom sensors (not the patient pipeline) |
| `test_RS_doses.m` | Kept | pre-sim RS-vs-ETHOS check; no recon involved |
| `run_single_field_simulation.m` / `run_second_field_simulation.m` | Wrapped | physics FUNCTIONS the classes call (not drivers) |

Common reason the "Kept" scripts resist migration: they parametrize k-Wave
internals `run_single_field_simulation` does not expose (Nt divisor, per-iteration
recon, dt beyond CFL, synthetic media) or change the recon grid size (cross-grid
gamma). Reproducing them means porting the inline forward/recon loop into the
class layer — deferred so nothing numerically unvalidated ships in the very
studies whose purpose is numerical sensitivity.

## k-Wave Simulation Parameters

| Parameter | Value |
|-----------|-------|
| PML size | 10 voxels (default) |
| CFL number | 0.3 |
| Time reversal BCs | Dirichlet boundary conditions |
| GPU acceleration | `use_gpu = true` |
| Typical runtime | ~23 min for 7 fields on 256×256×128 grid with GPU |

## Tissue / Grüneisen Methods (`sensor_placement_method`)

| Method | Description |
|--------|-------------|
| `'uniform'` | Single properties everywhere: Γ=1.0, c=1540 m/s, ρ=1000 kg/m³ |
| `'threshold_1'` | 9 tissue types by HU: air, lung, fat, water, blood, muscle, soft tissue, bone, metal |
| `'threshold_2'` | 4 tissue types by HU: water, fat, soft tissue, bone — **most commonly used** |

## Sensor Design

- **Geometry:** Flat rigid 10×10 cm planar array (not curved). Real arrays are rigid; tissue deforms, not the sensor.
- **Placement:** Anterior abdomen, avoiding all beam field jaw projections on the anterior surface.
- **Coupling:** Water fills gaps between sensor and body; no coupling concerns outside body.
- **Beam exclusion:** Uses divergence geometry to project jaw openings from isocenter onto anterior surface.

### Sensor Placement Methods (`sensor_placement_method`)

| Value | Description |
|-------|-------------|
| `'full_plane_anterior'` | Full YZ plane at `sensor_x_index` *(default)* |
| `'full_plane_lateral'` | Full XZ plane at `sensor_y_index` |
| `'spherical'` | Spherical shell geometry |

## Gotchas

- **Slice alignment streaking:** Misaligned CT/dose grids cause horizontal streaking artifacts in body masks. Always verify grid alignment after resampling.
- **GPU memory:** 256³ grids with complex tissue models can exhaust GPU memory. Monitor with `gpuDevice`; fall back to CPU if needed.

## Gamma Pass-Rate Individual Analysis (`study_pass_rates_individual.m`)

A re-runnable, no-k-Wave-needed (except the optional ensemble) batch analysis that, for an **arbitrary
list of dose files** (`CONFIG.dose_filenames`), runs two per-dose comparisons with individually
toggleable outputs. It is the per-dose counterpart to `study_gamma_pass_rates_chart.m` (which handles a
single recon pair), and it pulls its dose/difference panels and 10%-area/sensor views from
`run_standalone_simulation.m`.

### Per-dose batch model

Each `dose_<...>_B<beam>_<seg>.mat` filename is parsed (`parse_dose_selection`) for beam/segment/plan-type
and loaded on **both** CT images via `load_recon_dose_data` (`Mode='set'`), reusing
`load_field_set` / `find_field_by_ct`. The `_CT_k` token in the filename is **ignored** — both CT_1 and
CT_3 are always located. Per dose we cache `recon_CT1`, `recon_CT3`, `rs_CT1`, `rs_CT3`, the CBCT
geometry, RTPLAN gantry/meterset, the display density, and the (possibly grid-expanded) sensor mask.
Results land in a `RESULTS` struct saved to `CONFIG.output_file`.

### A1 vs A2 semantics

- **A1 — recon vs own truth (CT_1 only).** Reference = `recon_CT1`, target = `rs_CT1` (RayStation field
  dose). Measures how faithfully the detector/reconstruction reproduces the planned dose on the baseline
  CT. CT_3's recon is still loaded (used in A2) but A1 does not analyze it.
- **A2 — CT_1 recon vs CT_3 recon.** Reference = `recon_CT1`, target = `recon_CT3`. Measures the
  *change* the system sees between the reference and adapted anatomy — the actual change-detection task.

### Why `rs_dose` provides the 10% cutoff mask

Both A1 and A2 use the **same** low-dose evaluation mask: `rs_CT1 >= 0.10 * max(rs_CT1)`. Anchoring the
gamma evaluation region to the RayStation *truth* (not either recon) keeps the analyzed voxel set fixed
and physically meaningful across both comparisons, so A1 and A2 pass rates are computed over an identical
region and are directly comparable. (The noise-ensemble null is the one exception — see below.)

### Gamma sweep

For each criterion `n` in `CONFIG.gamma_n` (default `1:0.5:5`), gamma is `CalcGamma(ref, tgt, n, n,
'local', 0, 'limit', n*2, 'restrict', 1)` — i.e. **global** gamma (the `'local',0` flag), distance-to-
agreement search capped at `2n`, with `restrict=1` for speed. Pass rate = `100 * mean(gmap(mask) <= 1)`.
A red 90%-pass-rate reference line is drawn on every chart.

### Noise-ensemble null hypothesis (A2 only)

Enabled by `CONFIG.enable.noise_ensemble` (requires `CONFIG.convolution_kernel > 0`, since electronic
noise enters only through the pulse model). For each dose's **CT_1**:

1. **Forward sim once** — `build_forward_bundle(rs_CT1, cbct_CT1, gantry, ...)` runs one k-Wave forward
   simulation and caches the pre-noise sensor response, transfer-function terms, and grid/medium.
2. **Redraw ×N** — `parfor r = 1:num_realizations`: `redraw_noisy_deconv` adds a fresh noise draw (at
   `conv_noise_level`) + Wiener deconvolution, then `reconstruct_recon_dose` re-reconstructs.
3. **Null gamma** — each noisy recalculation is gamma-swept against the **original on-disk `recon_CT1`**
   (`gamma_sweep_pass_rates`, which uses a 10%-of-reference cutoff). This is the null: *same anatomy,
   noise only.* The mean ± std band is overlaid green on the A2 chart.

**Reading signal vs null:** the A2 signal curve (CT_1 vs CT_3) is the deterministic blue series; the
green null band is the pass rate expected from noise alone. Where the **signal sits below the null band**,
the CT_1→CT_3 difference is larger than noise — a *detectable, real* change. Where the signal stays within
the null band, the apparent difference is consistent with noise.

**Grid-expansion caveat:** sensor placement may water-pad/expand the k-Wave grid, so the bundle's recon
grid (`bundle.gridSize_orig`) can differ from the on-disk recon dims. The original `recon_CT1` is embedded
onto the bundle grid with `embed_on_grid(recon_CT1, bundle.gridSize_orig, bundle.embed_offset, 0)` (a
no-op when no expansion occurred) before the null gamma so both volumes share a grid.

**Single-field only:** the ensemble forward sim uses one gantry/sensor geometry, so it is meaningful for
single-field dose specs, not summed totals.

### Parallelization model

- **Gamma sweeps** — one flattened `parfor` over every `(dose × analysis × criterion)` job (CPU). The
  comparison reference/target volumes are broadcast to workers as `single` to limit broadcast memory;
  each job rebuilds the `CalcGamma` structs and scatters its pass rate back into a `[nComp × K]` matrix.
- **Noise ensemble** — `parfor` over realizations, one GPU per worker (`spmd gpuDevice(mod(labindex-1,
  ngpu)+1)`). The forward sim runs once per dose *before* the loop; the per-realization gamma sweep stays
  serial inside the loop (no nested parfor). The CPU gamma pool is started first and reused.

### 10%-area / sensor visualization

`plot_sensor_dose_planes` renders three orthogonal max-projections (coronal/sagittal/axial): CT density
as a grayscale mean-projection background, the dose mask (`rs_CT1 >= 10% max`) as semi-transparent blue,
and the ultrasound array as solid red. The dose-panel figure (`plot_truth_recon_diff_axial`) shows
`volA | volB | (volA − volB)` at the **axial slice containing the truth max**, with the diff panel using a
symmetric blue-white-red diverging colormap (alpha-ramped by magnitude) so over-/under-estimates are
distinguishable. All display volumes are embedded onto the sensor display grid so the red sensor contour
co-registers with the dose.
