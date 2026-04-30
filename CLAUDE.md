# CLAUDE.md — ETHOS Photoacoustic Dose Verification Pipeline

## Claude Code Context

### Loading Rules

When a user or file contents indicate a specific domain, load the relevant context based on the below description. 
To load a context file, ask the user: "Should I load `claude-[domain].md` for better context?" 
Or if it's obvious from the code/task, load it proactively and mention it.
Then read and incorporate that file's guidelines into your responses.

### How to Trigger Loading

You can:
1. **Explicitly ask the user** ("Looks like Simulation work—should I load the simulation context?")
2. **Load proactively** when you detect the domain (then say "I've loaded `claude-PIPELINE_CONTEXT.md` for context")
3. **Have the user request it** ("Load the DICOM context")

### Context Files

> - `CLAUDE-PIPELINE_CONTEXT.md` — step table, orchestrators, supporting files, full workflow
> - `CLAUDE-SIMULATION_CONTEXT.md` — k-Wave specifics, tissue models, sensor design (for simulation files)
> - `CLAUDE-DICOM_CONTEXT.md` — DICOM reference chains, MLC correction, segment explosion (for DICOM steps)
> - `CLAUDE-RayStation.md` — RayStation Scripting guidelines for step 1. 

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

- **DICOM:** `(x, y, z)` in mm; origin from `ImagePositionPatient`.
- **MATLAB arrays:** `(row, col, slice)` = `(Y, X, Z)`. **This mapping matters everywhere.**
- **Grid dims stored as** `[ny, nx, nz]` — watch for transposition bugs.
- **Spacing:** `[dx, dy, dz]` in mm; k-Wave requires meters: `dx_m = spacing(1) * 1e-3`.

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

## Visualization Preferences

- Subplots: **maximum 3 rows** on screen; 3–4 columns fine.
- Orthogonal views at max dose location or dose centroid (sagittal, coronal, transverse).
- Body contour overlays on dose colormaps; consistent colorbars.

## Prerequisites

- MATLAB R2022a+
- k-Wave Toolbox (http://www.k-wave.org)
- Image Processing Toolbox
- Parallel Computing Toolbox (optional)
