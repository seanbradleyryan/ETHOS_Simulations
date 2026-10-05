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
| 3 | `step3_analysis.m` | `results = step3_analysis(patient_id, session, config)` | Total-dose gamma + local SSIM, and the full per-segment study (tables, statistics, figures) from the Step 2.5 summary; everything saved to `AnalysisResults/` |

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
  Step 3 draws `noise_floor.mean_pass_rate` as the horizontal null band on the pass-rate-by-beam plot
  and, on the differentials, the point (`truth1_vs_recon1` pass) − (noise mean). The first
  (uncached) ~30 min compute is guarded by a `<noise_ensemble_cache>.lock` so only one instance pays it;
  siblings skip and pick up the cache on a later run.
- **Overlap with Step 2 (`CONFIG.metrics_overlap_step2`, default on):** so one allocation of GPUs + CPUs
  is never half idle, `pipeline_simulate` starts `step25_metrics_watcher.m` as a `batch` job with its own
  CPU pool (`metrics_watcher_workers`) before the Step 2 GPU `parfor`. Pool sizes are auto-detected when
  left `[]`: `num_parallel_workers` = `gpuDeviceCount` (one k-Wave worker per GPU; all CPUs if no GPU) and
  `metrics_watcher_workers` = CPUs − GPU workers − 1, with CPUs from `SLURM_CPUS_PER_TASK` /
  `SLURM_CPUS_ON_NODE`, else `feature('numcores')`. Every `metrics_watcher_poll_sec`
  it runs `step25_segment_metrics` (noise floor + summary off) on beams whose recon files are ALL on disk
  (beam = smallest unit, because `load_recon_dose_data` Mode `set` errors on a missing recon). After Step 2,
  the client runs the GPU noise floor (`metrics_noise_floor_only = true`) while the watcher finishes, creates
  the watcher's stop flag (`metrics_watcher_stop_<hash>_<pid>.flag` in the sim dir), waits for it, then runs
  the normal Step 2.5 call, which reads the folded scalars back, computes anything the watcher missed, and
  writes the summary. Needs `parcluster().NumWorkers >= 1 + metrics_watcher_workers + num_parallel_workers`;
  otherwise the watcher is shrunk or skipped. `field_index` is already sorted by [beam, segment], so beams
  finish progressively during Step 2.

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
- **Outputs** go to `AnalysisResults/<pid>/<session>/<method>/<hash>/`:
  - `step3_results.mat`
  - `beam_summary.csv`
  - `segment_statistics.csv`
  - `step3_console_log.txt`: diary, written to the session folder while running and moved here at the end
  - `figures/*.png`
- Step 3 never loads all field doses at once. The random panels load one beam at a time.

## Gotchas

- **Memory limits:** Never load all field doses simultaneously for large grids. Process one field at a time; save individually.
- **Gamma analysis cutoff:** Default 10% low-dose cutoff excludes voxels below 10% of max dose — intentional. Low-dose regions are clinically less relevant and noisy in reconstruction.
