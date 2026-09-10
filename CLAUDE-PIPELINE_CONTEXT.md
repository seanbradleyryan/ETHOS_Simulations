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
| 2.5 | `step25_segment_metrics.m` | `results = step25_segment_metrics(patient_id, session, config)` | Per-beam/segment gamma + local SSIM (4 study comparisons); folds raw masked results into each segment's **CT_1** recon file |
| 3 | `step3_analysis.m` | `results = step3_analysis(patient_id, session, config)` | Gamma analysis (3%/3mm), SSIM, visualization |

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
| `compute_local_ssim.m` | Local (per-voxel) SSIM map + masked mean (study-style, over the 10% region). NOT `step3`'s global `compute_dose_ssim` |

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
load('total_recon_dose_<hash8>.mat');  % Contains: total_recon, metadata, config_hash
                                       % <hash8> = compute_sim_config_hash(CONFIG); per-hash file
                                       % see SimulationResults/<id>/<session>/<method>/config_registry.json
```

> **Per-CT totals** (`total_dose_CT_*.mat`) are written by `step15_process_doses` whenever NPZ-derived
> field doses carry a CT label in their filename (e.g. `..._adapted_CT_1_B6_103.mat`).
> Sparse reconstruction: `reshape(full(ct_total_sparse), ct_total_dims)`.
> Legacy DICOM inputs (no CT label) only produce `total_rs_dose.mat`.

## Step 2.5 — Per-Segment Metrics (folded into CT_1 recon files)

`step25_segment_metrics.m` adapts `study_pass_rates_allsegments.m` into a pipeline step (no plots — raw
data only). It runs after Step 2 (in `pipeline_simulate`, guarded by `CONFIG.run_step25`) and can also be
run standalone like `step3_analysis`.

- **Why a separate step, not inside the sim `parfor`:** two of the four study comparisons are *cross-CT*
  (`recon_CT3 vs rs_CT1`, `rs_CT1 vs rs_CT3`), so they need BOTH the CT_1 and CT_3 reconstructions of the
  same beam/segment. The Step-2 `parfor` produces one field (one CT) at a time, so the pair isn't available
  there. Step 2.5 pairs the finished recons and parallelizes across CPUs (`parfor` over a beam's segments)
  while the GPU work is already done. Gamma is forced onto the CPU (`CalcGamma(...,'cpu',1)`).
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
  plotting script would draw `noise_floor.mean_pass_rate` as the horizontal null line on the pass-rate-by-
  beam plot and, on the differentials, the point (`truth1_vs_recon1` pass) − (noise mean). The first
  (uncached) ~30 min compute is guarded by a `<noise_ensemble_cache>.lock` so only one instance pays it;
  siblings skip and pick up the cache on a later run.

## Gotchas

- **Memory limits:** Never load all field doses simultaneously for large grids. Process one field at a time; save individually.
- **Gamma analysis cutoff:** Default 10% low-dose cutoff excludes voxels below 10% of max dose — intentional. Low-dose regions are clinically less relevant and noisy in reconstruction.
