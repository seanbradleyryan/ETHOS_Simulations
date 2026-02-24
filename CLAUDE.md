# CLAUDE.md — ETHOS Photoacoustic Dose Verification Pipeline

## Project Overview

MATLAB-based photoacoustic dose verification system for Varian ETHOS adaptive radiotherapy (pancreatic cancer). Simulates acoustic wave propagation from radiation dose deposition and reconstructs dose via time-reversal using the k-Wave toolbox. The goal is real-time dose monitoring during IMRT delivery.

**Core physics chain:** Radiation dose → initial pressure (p0 = Dose × Grüneisen × density) → k-Wave forward simulation → sensor data → time-reversal reconstruction → recovered dose

## Directory Structure

```
/mnt/weka/home/80030361/ETHOS_Simulations/       # Working root
├── EthosExports/[PatientID]/Pancreas/[Session]/  # Raw DICOM exports
│   └── sct/                                       # Sorted SCT + matched RT files (Step 0 output)
├── RayStationFiles/[PatientID]/[Session]/         # Field doses from RayStation
│   └── processed/                                  # Processed .mat files (Step 1.5 output)
├── SimulationResults/[PatientID]/[Session]/[method]/ # k-Wave outputs (Step 2 output)
├── AnalysisResults/[PatientID]/[Session]/          # Gamma/SSIM results (Step 3 output)
└── PipelineScripts/                                # All pipeline .m files
```

## Pipeline Steps

| Step | File | Function Signature | Purpose |
|------|------|--------------------|---------|
| 0 | `step0_sort_dicom.m` | `sct_dir = step0_sort_dicom(patient_id, session, config)` | Sort DICOM by reference chains, extract SCT series |
| 0.5 | `step05_fix_mlc_gaps.m` | `[path, n] = step05_fix_mlc_gaps(patient_id, session, config)` | Correct Halcyon dual-layer MLC minimum gaps |
| 1 | *Manual* | — | Export field doses from RayStation |
| 1.5 | `step15_process_doses.m` | `[fields, sct, total, meta] = step15_process_doses(...)` | Resample CT to dose grid, process per-field doses |
| 2 | `run_single_field_simulation.m` | `[recon, results] = run_single_field_simulation(...)` | k-Wave forward + time-reversal for one field |
| 3 | `step3_analysis.m` | `results = step3_analysis(patient_id, session, config)` | Gamma analysis (3%/3mm), SSIM, visualization |

**Supporting files:**
- `ethos_master_pipeline_pseudocode.m` — Orchestrator script with CONFIG and all helper functions
- `ethos_kwave_simulation.m` — Original monolithic simulation script (being superseded by modular pipeline)
- `run_standalone_simulation.m` — Self-contained single-run simulation for testing
- `determine_sensor_mask.m` — Physics-based flat sensor placement algorithm
- `CalcGamma.m` — Gamma index calculation (external dependency)
- `load_processed_data.m` / `visualize_processed_data.m` — Utilities

## Code Conventions

### Architecture Patterns
- **CONFIG-driven:** All parameters flow through a `CONFIG` struct (or `config` argument). No hardcoded magic numbers in processing logic — put them in CONFIG sections at the top.
- **Function signature:** Pipeline steps are `function output = stepN_name(patient_id, session, config)`. patient_id and session are always char arrays.
- **Stateless functions:** Designed for `parfor` safety. No persistent variables, no globals, no shared state.
- **Memory-conscious:** Field doses are processed and saved individually as `field_dose_XXX.mat` files, not held in a single large array. Use `-v7.3` for large .mat files.
- **Only the first function in a .m file is externally visible.** All helper functions go below the main function in the same file.

### Input Validation
Every pipeline function starts with input validation:
```matlab
if ~ischar(patient_id) && ~isstring(patient_id)
    error('stepN_name:InvalidInput', 'patient_id must be a string...');
end
patient_id = char(patient_id);  % Normalize to char
```
Validate config fields exist, set defaults for optional fields.

### Error Handling
Use qualified error IDs: `error('function_name:ErrorType', 'message', ...)`. Wrap major processing blocks in try/catch with informative messages.

### Documentation Style
Every function gets a header block with: PURPOSE, INPUTS (with field descriptions), OUTPUTS, ALGORITHM steps, EXAMPLE usage, DEPENDENCIES, and See also references.

### Naming
- Functions: `snake_case` (e.g., `step0_sort_dicom`, `determine_sensor_mask`)
- Variables: `camelCase` for local variables (e.g., `fieldDoseDir`, `patientID`), `snake_case` for struct fields (e.g., `.dose_Gy`, `.gantry_angle`)
- Config fields: `snake_case` (e.g., `config.dose_per_pulse_cGy`, `config.use_gpu`)
- Constants/CONFIG: `UPPER_CASE` struct name (`CONFIG`), lowercase fields

### Logging
Use `fprintf` with step/section prefixes:
```matlab
fprintf('[STEP 1.5] Processing field %d/%d...\n', i, n);
fprintf('  Grid dimensions: [%d x %d x %d]\n', dims);
```

## Key Technical Details

### Coordinate Systems
- DICOM patient coordinates: (x, y, z) in mm. Origin from DICOM ImagePositionPatient.
- MATLAB array indexing: (row, col, slice) = (Y, X, Z). **This mapping matters everywhere.**
- Grid dimensions stored as `[ny, nx, nz]` (rows, cols, slices) — watch for transposition bugs.
- Spacing: `[dx, dy, dz]` in mm. k-Wave needs meters: `dx_m = spacing(1) * 1e-3`.

### Tissue Models (Grüneisen Methods)
- `'uniform'` — Single property values everywhere (Γ=1.0, c=1540 m/s, ρ=1000 kg/m³)
- `'threshold_1'` — 9 tissue types based on HU thresholds (air, lung, fat, water, blood, muscle, soft tissue, bone, metal)
- `'threshold_2'` — 4 tissue types (water, fat, soft tissue, bone) — most commonly used

### Sensor Design
- **Flat rigid 10×10 cm planar array** — not curved. Real ultrasound arrays are rigid; tissue deforms, not the sensor.
- Positioned on anterior abdomen, avoiding all beam field jaw projections on the surface.
- Water fills gaps between sensor and body (no coupling concerns outside body).
- Beam exclusion uses divergence geometry to project jaw openings from isocenter onto anterior surface.

### k-Wave Specifics
- PML size: 10 voxels (default)
- CFL number: 0.3
- Time reversal with Dirichlet boundary conditions
- GPU acceleration via `use_gpu = true`
- Typical runtime: ~23 min for 7 fields on 256×256×128 grid with GPU

### DICOM Reference Chains
Step 0 matches files via: CT SeriesInstanceUID → RTSTRUCT ReferencedFrameOfReferenceSequence → RTPLAN ReferencedStructureSetSequence → RTDOSE ReferencedRTPlanSequence

### MLC Gap Correction (Halcyon)
Halcyon has a dual-layer MLC (proximal + distal). Step 0.5 enforces minimum gap of 0.5 mm with 0.4 mm expansion per side. Valid leaf range: [-140, 140] mm.

## Visualization Preferences

- **Subplots: Maximum 3 rows on screen.** 3–4 columns is fine.
- Orthogonal views typically show sagittal, coronal, and transverse planes at the max dose location or dose centroid.
- Body contour overlays on dose colormaps. Use consistent colorbars.

## Prerequisites

- MATLAB R2022a+
- k-Wave Toolbox (http://www.k-wave.org)
- Image Processing Toolbox
- Parallel Computing Toolbox (optional)
- matRad (at `/mnt/weka/home/80030361/MATLAB/Addons/matRad`)

## Common Operations

```bash
# Typical patient processing
# 1. Place raw DICOM export in EthosExports/[PatientID]/Pancreas/[Session]/
# 2. Run pipeline steps in order via ethos_master_pipeline_pseudocode.m
# 3. After Step 0.5, manually export field doses from RayStation
# 4. Re-run pipeline from Step 1.5 onward
```

```matlab
% Quick standalone test of a simulation
% Set paths in run_standalone_simulation.m CONFIG section, then run it

% Load and inspect processed data
load('sct_resampled.mat');  % Contains: sct_resampled struct
load('total_rs_dose.mat');  % Contains: total_rs_dose 3D array (Gy)
```

## Gotchas & Known Issues

1. **Slice alignment streaking:** Misaligned CT/dose grids cause horizontal streaking artifacts in body masks. Always verify grid alignment after resampling.
2. **Memory limits:** Never load all field doses simultaneously for large grids. Process one field at a time, save individually.
3. **MATLAB function visibility:** Only the first function in a `.m` file is callable from outside. If you need a helper externally, it must be its own file or the first function.
4. **DICOM collection failures:** `dicomCollection()` can be fragile with malformed DICOM. Fall back to manual `dicominfo()` iteration if needed.
5. **GPU memory:** 256³ grids with complex tissue models can exhaust GPU memory. Monitor with `gpuDevice` and fall back to CPU if needed.
6. **Gamma analysis cutoff:** Default 10% low-dose cutoff excludes voxels below 10% of max dose. This is intentional — low-dose regions are clinically less relevant and noisy in reconstruction.
