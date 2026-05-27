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
load('total_recon_dose.mat');    % Contains: total_recon, metadata
```

> **Per-CT totals** (`total_dose_CT_*.mat`) are written by `step15_process_doses` whenever NPZ-derived
> field doses carry a CT label in their filename (e.g. `..._adapted_CT_1_B6_103.mat`).
> Sparse reconstruction: `reshape(full(ct_total_sparse), ct_total_dims)`.
> Legacy DICOM inputs (no CT label) only produce `total_rs_dose.mat`.

## Gotchas

- **Memory limits:** Never load all field doses simultaneously for large grids. Process one field at a time; save individually.
- **Gamma analysis cutoff:** Default 10% low-dose cutoff excludes voxels below 10% of max dose — intentional. Low-dose regions are clinically less relevant and noisy in reconstruction.
