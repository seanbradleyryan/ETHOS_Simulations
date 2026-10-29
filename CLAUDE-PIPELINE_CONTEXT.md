# PIPELINE_CONTEXT.md — Pipeline Steps, Orchestrators & Workflow

> Referenced when working on orchestrator scripts or any file that coordinates pipeline steps.

## Pipeline Step Table

| Step | File | Signature | Purpose |
|------|------|-----------|---------|
| 0 | `step0_sort_dicom.m` | `sct_dir = step0_sort_dicom(patient_id, session, config)` | Sort DICOM by reference chains, extract SCT series |
| 0.5 | `step05_fix_mlc_gaps.m` | `[path, n] = step05_fix_mlc_gaps(patient_id, session, config)` | Correct Halcyon dual-layer MLC minimum gaps |
| 0.6 | `step06_explode_segments.m` | `exploded = step06_explode_segments(patient_id, session, config)` | Explode each beam's segments into individual 2-CP RTPLAN files |
| 1 | *Manual* | — | Import exploded RTPLANs into RayStation, recalculate and export field doses |
| 1.5 | `step15_process_doses.m` | `[fields, sct, total, meta] = step15_process_doses(...)` | Resample CT to dose grid, process per-field doses (**Windows work laptop** via `pipeline_compress.m`) |
| 2 | `run_single_field_simulation.m` | `[recon, results] = run_single_field_simulation(...)` | k-Wave forward + time-reversal for one field |
| 2.05 | `step205_run_full_field_simulation.m` | `results = step205_run_full_field_simulation(patient_id, session, config)` | Each CT's total RS dose (`total_dose_CT_1/3.mat`) simulated as ONE field (CT_1 full access, CT_3 blind on CT_1) with the Step 2 sensor mask; four Step 2.5 gamma + SSIM comparisons on the totals. Saves `full_field_dose_CT_*_<hash>.mat` beside the recons. To be compared with the summed per-field recon in analysis |
| 2.5 | `step25_segment_metrics.m` | `results = step25_segment_metrics(patient_id, session, config)` | Per-beam/segment gamma + local SSIM (4 study comparisons); folds raw masked results into each segment's **CT_1** recon file. `pipeline_simulate` now scores pairs inline in Step 2; step25 builds the summary/noise floor (run by Step 3) |
| 3 | `step3_analysis.m` | `results = step3_analysis(patient_id, session, config)` | Total-dose gamma + local SSIM, and the full per-segment study (tables, statistics, figures) from the Step 2.5 summary; everything saved to `AnalysisResults/<pid>/<session>/<label>/step3_analysis/` |

## Orchestrator Scripts

| Script | Platform | Runs | Notes |
|--------|----------|------|-------|
| `pipeline_setup.m` | **Windows** | Steps 0 / 0.5 / 0.6 | Run before RayStation export. Root: `C:/Users/80030361/Documents/ETHOS_Simulations` |
| `pipeline_compress.m` | **Windows** | Step 1.5 + `prepare_uploads` | Run after field dose DICOMs exported from RayStation. Input: `RayStationFiles/[PatientID]/[Session]/` |
| `pipeline_simulate.m` | **Linux cluster** | Steps 2 / 3 | Run after processed `.mat` files transferred from Windows |

## Supporting Files

| File | Purpose |
|------|---------|
| `run_standalone_simulation.m` | Self-contained single-run simulation for testing |
| `determine_sensor_mask.m` | Physics-based flat sensor placement algorithm |
| `create_acoustic_medium.m` | Builds k-Wave medium struct from CT HU values |
| `apply_element_averaging.m` | Post-processing averaging over sensor elements |
| `find_optimal_kwave_size.m` | Selects efficient grid dimensions for k-Wave FFT |
| `convert_dose.m` | Unit conversion utilities for dose arrays |
| `run_medium_comparison.m` | Compares reconstruction quality across Grüneisen methods |
| `test_time_dependence.m` | Time-dependence sensitivity testing |
| `plot_standalone_results.m` | Visualization helper for standalone runs |
| `CalcGamma.m` | Gamma index calculation (external dependency) |
| `load_processed_data.m` | Loads previously processed dose/CT data |
| `load_recon_dose_data.m` | Loads computed recon doses (single/set/total) + RS, ETHOS truth, CBCT, RTPLAN stats — no re-sim. See CLAUDE.md "Analysis Utilities" |
| `load_beam_total_set.m` | Beam mode: one beam's summed RS dose + beam recon on CT_1/CT_3, shaped like `load_recon_dose_data` `'set'` output (Steps 2.5 / 3) |
| `get_simulation_dir.m` | Path of a session's simulation data: `RayStationFiles/<pid>/<session>/processed/simulation_results/<method>/` |
| `config_label.m` | AnalysisResults folder name for a config: `'canonical'` or the fields differing from `get_default_config` |
| `misc/migrate_results_layout.m` | One-off move of pre-overhaul `SimulationResults/` + hash-named analysis folders to the new layout |
| `compute_local_ssim.m` | Local (per-voxel) SSIM map + masked mean over the 10% region. The only SSIM in the pipeline (Steps 2.5 and 3) |
| `studies/study_pass_rates_allsegments.m` | Superseded by Step 3 (all its analyses now run there); kept for reference |

## Full Workflow (Common Operations)

```
1. [WINDOWS]  Place raw DICOM export in EthosExports/[PatientID]/Pancreas/[Session]/
2. [WINDOWS]  Set patient/session in pipeline_setup.m CONFIG → run it (Steps 0/0.5/0.6)
3. [MANUAL]   Import RTPLAN files from Raystation_Input/ into RayStation
4. [MANUAL]   Recalculate dose for each exploded-segment plan in RayStation
5. [MANUAL]   Export field doses as Plan_Field*_Beam*_B*_S*.dcm to
              C:\Users\80030361\Documents\ETHOS_Simulations\RayStationFiles\[PatientID]\[Session]\
6. [WINDOWS]  Set patient/session in pipeline_compress.m CONFIG → run it
              (Step 1.5 + prepare_uploads — processes DICOMs, packages .mat files)
7. [MANUAL]   Transfer processed .mat files from Windows to Linux cluster
8. [CLUSTER]  Set patient/session in pipeline_simulate.m CONFIG → run it (Steps 2/3)
```

## Loading Processed Data

```matlab
load('sct_resampled.mat');       % Contains: sct_resampled struct
load('total_rs_dose.mat');       % Contains: total_rs_dose 3D array (Gy) — sum of all fields
load('total_dose_CT_1.mat');     % Contains: ct_total / ct_total_sparse + ct_total_dims — total dose from CT_1 fields only
load('total_dose_CT_3.mat');     % Same, for CT_3 fields (label matches ct_label in field filenames)
load('total_dose_CT_1_B6.mat');  % Contains: beam_total_sparse + beam_total_dims (or beam_total), num_segments
                                 % — sum of beam 6's segments on CT_1 (one file per CT label x beam)
load('total_recon_dose_<hash8>.mat');  % Contains: total_recon, metadata, config_hash
                                       % <hash8> = compute_sim_config_hash(CONFIG); per-hash file
                                       % see <simulation dir>/config_registry.json
```

## Results Layout

| What | Where | Named by |
|---|---|---|
| Step 1.5 doses, CBCTs, totals | `RayStationFiles/<pid>/<session>/processed/` | — |
| Step 2/2.05/2.5 data (recons, `total_recon_dose`, `sensor_mask`, manifests, `segment_metrics_summary`, `full_field_dose`, noise-ensemble cache, `config_registry.json`, `simulation_log.txt`, `field_status/`) | `get_simulation_dir(...)` = `processed/simulation_results/<method>/` | `<hash>` in filename |
| Step 3 | `AnalysisResults/<pid>/<session>/<label>/step3_analysis/` | config label |
| `verify_pipeline_simulate` | `AnalysisResults/<pid>/<session>/<label>/verify_pipeline_simulate/` | config label |
| `pipeline_simulate` RESULTS + timing (per session) | `AnalysisResults/<pid>/<session>/<label>/pipeline_simulate/` | config label |
| `verify_pipeline_compress_output`, `verify_plot_dose_from_npz` | `AnalysisResults/<pid>/<session>/<script>/` | none (pre-simulation) |

- **Label** (`utils/config_label.m`): `'canonical'` when every config-hash field equals `get_default_config`;
  otherwise `<field>.<value>` per differing field, joined by `__` (e.g. `reconstruction_method.das`,
  `downscale_factor.2__gruneisen_method.threshold_1`). `unset` = field absent in that run. Compares only the
  hash fields, so label ↔ hash is 1:1. Analysis knobs (gamma criteria etc.) are NOT in the label.
- Functions that only hold a hash (step3, verify) get the label from the registry:
  `config_label(hashtoconfig(hash, config))`. pipeline_simulate also stores `.label` per hash in `config_registry.json`.
- Changing a value in `get_default_config` changes which hash is `canonical`; the old canonical hash then gets
  a diff label.
- Pre-overhaul data: run `misc/migrate_results_layout.m` once (dry run first).

> **Per-CT totals** (`total_dose_CT_*.mat`) are written by `step15_process_doses` whenever NPZ-derived
> field doses carry a CT label in their filename (e.g. `..._adapted_CT_1_B6_103.mat`).
> Sparse reconstruction: `reshape(full(ct_total_sparse), ct_total_dims)`.
> Legacy DICOM inputs (no CT label) only produce `total_rs_dose.mat`.

## Reconstruction Mode (`CONFIG.reconstruction_mode`, default `'segment'`)

Set in `get_default_config`; part of the config hash (`'segment'` is dropped from the hash, so
pre-existing segment caches keep their hashes; `'beam'` gets its own hash/files).

- **`'segment'`**: one Step 2 sim + Step 2.5 gamma/SSIM per segment (thousands). Original behavior.
- **`'beam'`**: one per (CT label, beam) on the beam's summed dose `total_dose_<CT>_B<n>.mat`
  (from `pipeline_compress` / Step 1.5); ~15–20 sims.
  - `pipeline_simulate` swaps `field_index` for `list_beam_doses` output. Recon / sim / lock /
    status files are named `total_dose_<CT>_B<n>_recon_<hash>.mat` etc.
  - `load_beam_dose_file` takes metadata (gantry, CT label, meterset) from one segment.
    step15 stores the RTPLAN **beam** meterset on every segment, so that is already the beam MU.
  - **Mode comes from the hash.** Steps 2.5 and 3 read `reconstruction_mode` from the hash's
    registered config (`hashtoconfig`), not from `config.reconstruction_mode` (warns if they differ).
    Reason: a config left at the `'segment'` default hashes to the old segment run, whose files
    exist, so the loaders silently used segment data. `pipeline_simulate` passes `config_hash` to
    Step 2.5 (as it does for 2.05 and 3). Step 2.5 without `config_hash` still uses the config's mode.
  - Step 2.5 runs the four comparisons once per beam on the beam totals (`seg_num` 0) and folds into
    the CT_1 beam recon. All beams are loaded first, then one `parfor` runs over the beams (not over
    a beam's segments, which would be one item). `metrics_plan_type` is not applied. Beam data is loaded by
    `utils/load_beam_total_set.m` (same shape as `load_recon_dose_data` `'set'` output).
  - Step 3 Part A is unchanged. Part B prints a per-beam table (one value per beam,
    All = mean ± SE over beams) and draws the by-beam plots without per-beam error bars. SNR is
    2 values per beam (CT_1 + CT_3), with one panel per beam (`beam_panel_B<n>.png`). It skips the
    segment statistics (t-test, floor, correlations), the correlation plot and `segment_statistics.csv`.
  - `verify_pipeline_simulate` is skipped.
  - Not exactly equal to the summed segment recons: each sim adds one noise draw (sum of N segments
    = N draws), and the LS gain and >1% dose-region mask are nonlinear.

## Step 2.5 — Per-Segment Metrics (folded into CT_1 recon files)

`step25_segment_metrics.m` adapts `study_pass_rates_allsegments.m` into a pipeline step (no plots — raw
data only). `pipeline_simulate` no longer calls it: Step 3 runs it when `segment_metrics_summary_<hash>.mat`
is missing, or run it standalone like `step3_analysis`.

- **Inline scoring in Step 2 (`CONFIG.run_segment_metrics`):** two of the four comparisons are *cross-CT*,
  so they need BOTH the CT_1 and CT_3 recons of a segment. In `pipeline_simulate`, after a field is saved
  and marked `complete`, `fold_segment_metrics` checks the partner field (same filename, CT number
  swapped). If the partner's status file says `complete`, this worker scores the pair and folds it into
  the CT_1 recon, using the same math and schema as Step 2.5. That CPU work overlaps other workers' GPU
  sims. Gamma is forced onto the CPU (`CalcGamma(...,'cpu',1)`). A run that scores ≥1 pair deletes the
  stale `segment_metrics_summary_<hash>.mat` so Step 3 rebuilds it.
- **What step25 still does:** builds the summary rollup, the noise floor and `full_vs_sum`, and scores any
  pair the inline pass missed (pairs cached from older runs, rare filesystem races).
- **Four comparisons per segment** (reference builds the 10% eval mask): `truth1_vs_truth3`,
  `truth1_vs_recon1`, `truth1_vs_recon3`, `truth3_vs_recon3`. For each, BOTH the global gamma index
  (`CONFIG.gamma_dose_pct`/`gamma_dist_mm`) and the local SSIM map are computed.
- **Storage — masked-region-only, folded into the CT_1 recon file:** the per-segment result is appended as
  a `segment_metrics` variable INTO that segment's CT_1 recon `.mat`
  (`<base_CT_1>_recon_<hash>.mat`) — **nothing is written to the CT_3 recon file**. Only voxels inside each
  comparison's 10% mask are kept: `comparison(d).mask_idx` (uint32 linear indices), `.gamma_vals` /
  `.ssim_vals` (single), `.gamma_pass_rate` / `.ssim_mean` (%), alongside `vol_size`, `spacing`, the applied
  LS `recon_ct1_gain`/`recon_ct3_gain`, and gamma params. Re-expand a dense map on demand with
  `M = nan(sm.vol_size); M(c.mask_idx) = c.gamma_vals;`. This keeps the raw index/SSIM arrays without the
  hundreds of GB dense volumes would cost.
- **Resumable:** a segment whose CT_1 recon file already carries `segment_metrics` is skipped (its scalars
  are still read back into the summary), so the step can be re-run piecewise as more fields finish — exactly
  like the per-field sim. `CONFIG.metrics_overwrite = true` forces recompute.
- **Cross-instance safety:** because several `pipeline_simulate` copies can reach Step 2.5 at once, each
  segment's fold is guarded by an atomic `<ct1_recon>.metrics_lock` directory (same `mkdir` trick as
  `claim_field`) so two instances never append to the same recon `.mat` concurrently. A crashed run can
  leave a stale `*.metrics_lock` dir that blocks that one segment — delete it (or re-run with
  `metrics_overwrite`) to reclaim. No stale-overtake timer (unlike Step 2), since a segment is cheap to redo.
- **Summary rollup:** `segment_metrics_summary_<hash>.mat` (per-beam & pooled mean/std of gamma pass % and
  mean local SSIM %, per comparison) is written beside the recons — the data the study's plots are built
  from. Segments needing a CT_1/CT_3 pair that is missing are skipped (matching the study).
- **Noise-only null floor** (`CONFIG.metrics_noise_floor`, default on): the null hypothesis is the gamma
  pass rate of the CT_1 truth vs a NOISE-ONLY reconstruction. `noise_ensemble_error_bars` runs an ensemble
  of noise-only recons and returns the mean ± std pass rate; it is computed **once per session** (the util
  caches by the sim config hash, so every beam/segment shares it) and stored on the summary as
  `results.noise_floor` (`.mean_pass_rate`, `.std_pass_rate`, `.num_samples`, …). `metrics_noise_minutes`
  (30) is the ensemble `TimeBudgetMin`. The reference geometry/summed truth come from
  `load_recon_dose_data(Mode='total')`. Since Step 2.5 draws no figures, the floor is only *saved*; a
  Step 3 draws `noise_floor.mean_pass_rate` as the horizontal null band on the pass-rate-by-beam plot
  and, on the differentials, the point (`truth1_vs_recon1` pass) − (noise mean). The first
  (uncached) ~30 min compute is guarded by a `<noise_ensemble_cache>.lock` so only one instance pays it;
  siblings skip and pick up the cache on a later run.

## Step 3 — Analysis (totals + per-segment study)

`step3_analysis.m` runs after Step 2.5 (`CONFIG.run_step3`, default off in `pipeline_simulate`) or
standalone. It absorbed every analysis from `study_pass_rates_allsegments.m`.

- **Part A, total doses** (loaded via `load_recon_dose_data(Mode='total')`, hash pinned by
  `config.config_hash`):
  - `ethos_vs_rs`: ETHOS RTDOSE truth vs RS **per-CT** total `total_dose_<analysis_ethos_ct_label>.mat`
    (default `CT_1`). Not `total_rs_dose`, which sums both CTs' fields (≈2 plans).
  - `rs_vs_recon`: `total_rs_dose` vs `total_recon` (both sums of every field). The recon is LS-scaled
    when `metrics_normalize`.
  - Gamma and local SSIM use the same settings as Step 2.5: global `CalcGamma`, `limit = 2*dist_mm`,
    reference-only 10% mask, `compute_local_ssim`. `CalcGamma` widths are passed in array order
    `[dy dx dz]`.
  - ETHOS is resampled by **patient position** (RTDOSE IPP/PixelSpacing/GridFrameOffsetVector → dose-grid
    origin/spacing), not by array size.
  - Figures use orthogonal views (transverse/coronal/sagittal) at the reference max-dose voxel, with the
    CT_1 body contour and the real sensor footprint from `sensor_mask_<hash>.mat` (skipped when that mask
    is not on the dose grid).
- **Part B, per-segment study:** reads `segment_metrics_summary_<hash>.mat` and runs Step 2.5 first if it
  is missing. For **both** metrics (gamma pass %, local SSIM %) it produces:
  - per-beam mean ± SE tables, plus a pooled "All" row
  - pass-rate, differential and own-CT fidelity plots, with the noise-only band drawn for gamma only
  - segment statistics: paired one-sided t-test recon1 vs recon3; % above the floor
    (`analysis_noise_floor_pct`, default = Step 2.5 noise-floor mean; gamma only); Pearson R² and
    Spearman ρ of (recon1 − recon3) vs truth1_vs_truth3 and vs CT_1 SNR. All base MATLAB (`betainc`,
    `corrcoef`), no Statistics Toolbox.

  It also produces per-beam `noise_stats.snr` and `analysis_n_random` random segment panels
  (truth | recon | folded gamma map | folded SSIM map). No metric is recomputed here.
- **Outputs** go to `AnalysisResults/<pid>/<session>/<label>/step3_analysis/`:
  - `step3_results.mat`
  - `beam_summary.csv`
  - `segment_statistics.csv`
  - `step3_console_log.txt`: diary, written to the session folder while running and moved here at the end
  - `figures/*.png`
- Step 3 never loads all field doses at once. The random panels load one beam at a time.

## Gotchas

- **Memory limits:** Never load all field doses simultaneously for large grids. Process one field at a time; save individually.
- **Gamma analysis cutoff:** Default 10% low-dose cutoff excludes voxels below 10% of max dose — intentional. Low-dose regions are clinically less relevant and noisy in reconstruction.
