# SIMULATION_CONTEXT.md — k-Wave, Tissue Models & Sensor Design

> Referenced when working on: `run_single_field_simulation.m`, `run_standalone_simulation.m`, `create_acoustic_medium.m`, `determine_sensor_mask.m`, `apply_element_averaging.m`, `find_optimal_kwave_size.m`, `run_medium_comparison.m`, `test_time_dependence.m`

## k-Wave Simulation Parameters

| Parameter | Value |
|-----------|-------|
| PML size | 10 voxels (default) |
| CFL number | 0.3 |
| Time reversal BCs | Dirichlet boundary conditions |
| GPU acceleration | `use_gpu = true` |
| Typical runtime | ~23 min for 7 fields on 256×256×128 grid with GPU |

## Acoustic Medium — single source of truth (MANDATORY)

**All medium construction MUST go through `create_acoustic_medium.m`.** It is the only
builder that applies the `threshold_2` **in-body air row** (low-HU voxels inside the body
contour → real air: ρ≈1.2, c=343, Γ=0, so air generates no PA signal and is later masked
from the recon via `recon_dose(density < 100) = 0`). Do **not** hand-roll a medium from the
HU thresholds in a driver.

- **Coupling bath:** wrap it as `build_medium_with_bath` (create_acoustic_medium, then force
  outside-body / couch voxels to the `uniform_*` medium), exactly as `run_standalone_field`
  and `noise_ensemble_error_bars` do.
- **Known offender:** `study_pass_rates_individual.m` has a local `create_medium` that omits
  the air row. It disagrees with the engine and must not be copied into new code. (This was
  the cause of the noise-null pass rate reading low until the null was pointed back at
  `create_acoustic_medium`.)
- When adding a new simulation/analysis script, take the medium from `create_acoustic_medium`
  (+ bath) — never reimplement HU→property assignment.

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
| `'spherical'` | Spherical shell geometry. The grid is water-padded and the data box re-centered so the sphere **circumscribes** the whole data volume (no p0 clipped in the corners). This expansion can grow the grid substantially (~5x voxels for a cube) — use `downscale_factor` if memory-bound. Handled in `run_single_field_simulation.m`. |

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
