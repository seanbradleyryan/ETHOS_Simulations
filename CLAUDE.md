# CLAUDE.md — ETHOS Photoacoustic Dose Verification Pipeline

## Claude Code Context

### Context Loading (MANDATORY — check before every first response)

Scan the user's first message for these keywords. If matched, read that file BEFORE responding. Do NOT ask — just read it silently and proceed.

| Keywords in message | File to read |
|---|---|
| simulate, k-Wave, kgrid, sensor, acoustic, forward sim, run_single, run_standalone, pipeline_simulate | `CLAUDE-SIMULATION_CONTEXT.md` |
| RayStation, RS script, scripting, IronPython | `CLAUDE-RayStation.md` |
| DICOM, MLC, segment, RTPLAN, RT struct, sort, explode | `CLAUDE-DICOM_CONTEXT.md` |
| pipeline, step table, orchestrator, workflow, load_recon | `CLAUDE-PIPELINE_CONTEXT.md` |
| SimConfig, SimCase, Simulator, StudyRunner, SimResult, SimPlotter, +ethos, class layer | `CLAUDE-SIMULATION_CONTEXT.md` |

### Context Files

> - `CLAUDE-PIPELINE_CONTEXT.md` — step table, orchestrators, supporting files, full workflow
> - `CLAUDE-SIMULATION_CONTEXT.md` — k-Wave specifics, tissue models, sensor design (for simulation files)
> - `CLAUDE-DICOM_CONTEXT.md` — DICOM reference chains, MLC correction, segment explosion (for DICOM steps)
> - `CLAUDE-RayStation.md` — RayStation Scripting guidelines for step 1.

## MATLAB Editing Rules
- Do NOT run Python scripts or bash commands to verify or analyze `.m` files
- Do NOT attempt to lint or parse MATLAB using non-MATLAB tools
- Make edits directly based on code reading; do not run verification scripts
- If you cannot verify correctness without running MATLAB, state that explicitly instead

## Project Summary

IRAI (Ionizing Radiation Acoustic Imaging) captures sound waves from radiation dose deposition using a 2D ultrasound array, then reconstructs 3D dose via time-reversal. This pipeline simulates a live IRAI system on real ETHOS CBCT patient data by detecting differences in signals between two dose distributions, and comparing that to their different truth distributions and phantom data. 

**Core physics chain:** Radiation dose → initial pressure (`p0 = Dose × Grüneisen × density`) → k-Wave forward simulation → sensor data → time-reversal reconstruction → recovered dose

## Directory Structure

```
C:/Users/80030361/ETHOS_Simulations/
├── EthosExports/[PatientID]/Pancreas/[Session]/   # Raw DICOM exports
│   └── sct/                                        # Sorted SCT + matched RT files (Step 0 output)
├── Raystation_Input/[PatientID]/[Session]/         # Exploded-segment RTPLAN files (Step 0.6 output)
├── RayStationFiles/[PatientID]/[Session]/          # Field dose DICOMs from RayStation
│   └── processed/                                   # Processed .mat files (Step 1.5 output)
├── SimulationResults/[PatientID]/[Session]/[method]/ # k-Wave outputs (Step 2)
├── AnalysisResults/[PatientID]/[Session]/           # Gamma/SSIM results (Step 3)
└── PipelineScripts/                                 # All pipeline .m files
```

## Architecture Patterns

- **Class layer (`+ethos/`):** `SimConfig` (config presets + `toStruct` bridge + `hash`), `SimCase` (dose+geometry+medium+sensor bundle), `Simulator` (load-or-simulate `resolve`), `SimResult` (universal recon+truth container), `Analysis`/`StudyRunner` (gamma default 3%/3mm + SSIM over the same 10%-of-truth mask), `SimPlotter` (all figures). Value classes are `parfor`-safe; the classes **wrap** the functions below, never replace them. See `CLAUDE-SIMULATION_CONTEXT.md`.
- **CONFIG-driven:** All parameters in a `CONFIG` struct. No hardcoded magic numbers in logic.
- **Function signature:** `function output = stepN_name(patient_id, session, config)` — `patient_id` and `session` are always `char`.
- **Stateless functions:** `parfor`-safe. No persistent variables, no globals, no shared state.
- **Memory-conscious:** Field doses are large; process one at a time, never load all simultaneously.
- **MATLAB visibility:** Only the first function in a `.m` file is externally callable. Helpers go below main in the same file.
- **HIPAA:** Code must be executed on a remote device; do not attempt local testing.

## Naming Conventions

| Scope | Style | Example |
|---|---|---|
| Functions | `snake_case` | `step0_sort_dicom`, `determine_sensor_mask` |
| Local variables | `camelCase` | `fieldDoseDir`, `patientID` |
| Struct fields | `snake_case` | `.dose_Gy`, `.gantry_angle` |
| Config fields | `snake_case` | `config.dose_per_pulse_cGy`, `config.use_gpu` |
| CONFIG struct name | `UPPER_CASE` | `CONFIG` |

## Coordinate System

- **DICOM:** `(x, y, z)` in mm; origin from `ImagePositionPatient`. HFS (head-first supine) convention assumed throughout.
- **MATLAB arrays:** `(row, col, slice)` = `(Y, X, Z)`. **This mapping matters everywhere.**
- **Grid dims stored as** `[ny, nx, nz]` — watch for transposition bugs.
- **Spacing:** `[dx, dy, dz]` in mm; k-Wave requires meters: `dx_m = spacing(1) * 1e-3`.
- **Axis directions** (HFS, all in voxel-index space):
    - `X` (cols, dim 2): higher index = patient LEFT, lower = RIGHT.
    - `Y` (rows, dim 1): **lower** index = ANTERIOR, higher = POSTERIOR.
    - `Z` (slices, dim 3): **higher** index = SUPERIOR (cranial/head), lower = INFERIOR (caudal/feet). `dz > 0`.
  When code talks about "sliding inferior" it means *toward lower iz* (feet). For abdominal/pancreas placement the ribs sit on the superior side, so a sensor that drifts to higher iz to clear an exclusion zone will be acoustically blocked.

## Code Conventions

**Input validation** (every pipeline function starts with this):
```matlab
if ~ischar(patient_id) && ~isstring(patient_id)
    error('stepN_name:InvalidInput', 'patient_id must be a string...');
end
patient_id = char(patient_id);  % Normalize to char
```

**Error IDs:** `error('function_name:ErrorType', 'message', ...)`

**Logging:**
```matlab
fprintf('[STEP 1.5] Processing field %d/%d...\n', i, n);
fprintf('  Grid dimensions: [%d x %d x %d]\n', dims);
```

**Documentation header** (every function): PURPOSE, INPUTS (with field descriptions), OUTPUTS, ALGORITHM steps, EXAMPLE usage, DEPENDENCIES, See also.

## Analysis Utilities

### `load_recon_dose_data.m` — load computed recon doses for analysis

Reusable loader that pulls **already-computed** reconstructed doses from the pipeline outputs so the
standalone-simulation / comparison / study analyses can be re-run **without re-running k-Wave**. Each
dose is associated with its RayStation "truth" field dose, the ETHOS RTDOSE truth (resampled onto the
dose grid), the CBCT geometry, and the RTPLAN statistics.

**Signature:** `out = load_recon_dose_data(patient_id, session, config, Name, Value, ...)`

- `config` — needs `.working_dir`; uses `.treatment_site` (default `'Pancreas'`) and
  `.gruneisen_method` (auto-detected when a single method folder exists). May also carry the full
  simulation CONFIG so the config hash can be computed exactly.

**Name/value options:**

| Option | Values | Notes |
|---|---|---|
| `'Mode'` | `'total'` (default) \| `'set'` \| `'single'` | total recon, filtered set, or one field |
| `'Beam'` | scalar | required for `'single'`; filters `'set'` |
| `'Segment'` | scalar | optional filter; narrows `'single'` when a beam has several segments |
| `'PlanType'` | `'any'` (default) \| `'adapted'` \| `'reference'` | filter |
| `'CTLabel'` | `'any'` (default) \| `'CT_1'` \| `'CT_3'` | filter |
| `'Hash'` | 8-char hex | explicit config-hash override |
| `'IncludeRS'` / `'IncludeEthos'` / `'IncludeCBCT'` | logical (default `true`) | skip heavy items for speed |

**Config hash:** auto-discovers the recon hash on disk; otherwise uses the hash computed from CONFIG;
`'Hash'` overrides. Errors (listing the available hashes / method folders) on ambiguity or when nothing
matches, and warns on an embedded-hash mismatch.

**Returns** `out` with common `.patient_id/.session/.gruneisen_method/.config_hash`, `.metadata`
(grid `dimensions`/`spacing`/`origin` + `beam_metadata`), and `.ethos_truth` (when `IncludeEthos`):
- **`'total'`**: `.recon_dose`, `.rs_dose`, `.cbct.CT_1` / `.cbct.CT_3`.
- **`'set'` / `'single'`**: `.fields` struct array (1 element per matched field), each with
  `.recon_dose`, `.rs_dose`, `.cbct`, `.rtplan` (`beam_num`, `seg_num`, `plan_type`, `ct_label`,
  `beam_name`, `gantry_angle`, `meterset`, `isocenter`, `jaw_x`, `jaw_y`), `.source_mat_filename`,
  `.recon_file`.

`recon_dose`, `rs_dose`, and `ethos_truth` all share the dose grid (`metadata.dimensions`), so results
drop straight into gamma/SSIM with no further resampling.

```matlab
config.working_dir      = '/mnt/weka/home/80030361/ETHOS_Simulations';
config.treatment_site   = 'Pancreas';
config.gruneisen_method = 'threshold_2';

T = load_recon_dose_data('1194203','Session_1',config);                                          % summed total
S = load_recon_dose_data('1194203','Session_1',config,'Mode','single','Beam',15,'Segment',112);  % one field
A = load_recon_dose_data('1194203','Session_1',config,'Mode','set','PlanType','adapted','CTLabel','CT_3'); % filtered set
```

Reuses `compute_sim_config_hash`, `list_processed_field_doses`, and `load_field_dose_file`. Reads
recon outputs from `SimulationResults/[PatientID]/[Session]/[method]/` (`<base>_recon_<hash>.mat`,
`total_recon_dose_<hash>.mat`) and RayStation/CBCT inputs from `RayStationFiles/[PatientID]/[Session]/processed/`.

### `study_pass_rates_individual.m` — per-dose batch gamma analysis

Top-level script (editable `CONFIG` block) that runs gamma pass-rate analysis for an **arbitrary list**
of dose files (`CONFIG.dose_filenames`). Each file selects a beam/segment; the script loads that field's
pre-computed reconstruction on **both** CT images (CT_1 and CT_3) plus the RayStation truth via
`load_recon_dose_data` (`Mode='set'`) — no k-Wave needed except for the optional ensemble. Per dose it
runs two independently toggleable analyses:

- **A1** — recon vs its own RayStation truth (CT_1): per-dose detector accuracy.
- **A2** — CT_1 recon vs CT_3 recon: adapted-vs-reference change detection, with an optional
  **noise-ensemble null** (forward sim once → redraw noise → reconstruct ×N → gamma vs the original
  CT_1 recon) overlaid as the no-change noise floor.

Each analysis emits (gated by `CONFIG.enable.{A1,A2,dose_panels,sensor_view,noise_ensemble}`): a gamma
pass-rate vs `n%/n mm` chart, a 3-panel axial figure (volA | volB | signed difference) at the truth
max-dose slice, and a 10%-dose-area + sensor-model view. Gamma sweeps run as one flattened `parfor` over
all `(dose × analysis × criterion)` jobs; the ensemble runs `parfor` over realizations with one GPU per
worker. Depends on `load_recon_dose_data`, `CalcGamma`, and `determine_sensor_mask`. See
`CLAUDE-SIMULATION_CONTEXT.md` for the full method. (Remote/HIPAA execution; do not run locally.)

## Visualization Preferences

- Subplots: **maximum 3 rows** on screen; 3–4 columns fine.
- Orthogonal views at max dose location or dose centroid (sagittal, coronal, transverse).
- Body contour overlays on dose colormaps; consistent colorbars.

## Other

- CT 2 and CT 3 are synonyms now. Keep the names consistent with the current convention in the code, but in outside conversations, these words refer to the same image. 

## Prerequisites

- MATLAB R2022a+
- k-Wave Toolbox (http://www.k-wave.org)
- Image Processing Toolbox
- Parallel Computing Toolbox (optional)
