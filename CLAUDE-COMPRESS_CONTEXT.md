# CLAUDE-COMPRESS_CONTEXT.md — `pipeline_compress` and Daughter Scripts

Documentation for the **dose-processing stage** of the ETHOS photoacoustic pipeline:
`pipeline_compress.m` and everything it calls. This stage runs on the Linux cluster
**after** the RayStation field-dose export and **before** `pipeline_simulate.m`. It turns
raw RayStation NPZ field doses into the per-field / total / CBCT-resampled `.mat` files
that k-Wave consumes.

> **Where the outputs go:** `RayStationFiles/[PatientID]/[Session]/processed/`
> These `.mat` files are the direct input to `pipeline_simulate.m`. There is **no**
> separate upload-packaging step.

---

## 1. File map

| File | Role |
|---|---|
| `pipeline_compress.m` | **Orchestrator script** (not a function). Loops patients × sessions, runs Steps 1.4 and 1.5, prints a summary. |
| `pipeline/step14_npz_to_mat.m` | **Step 1.4.** Converts RayStation `dose_*.npz` → per-field `raw_field_dose` `.mat`. Includes a self-contained NPY/NPZ reader. |
| `pipeline/step15_process_doses.m` | **Step 1.5.** Loads field doses, resamples both CBCTs to the dose grid, builds per-CBCT tissue/body/couch masks, masks + saves each field, accumulates totals. The bulk of the work. |
| `utils/resolve_input_dir.m` | Helper. Picks input directory: `RayStationFiles` first, `EthosExports/.../sct` fallback. |

`pipeline_compress.m` also defines several **local helper functions** (below `main`) used only
by the orchestrator; Steps 1.4 and 1.5 each carry their own private helpers.

---

## 2. Data flow (one patient/session)

```
RayStation export (Python, calc_beam_plan_doses.py)
        │  dose_<id>_<session>_<plan>_<CT_k>_B<beam>_<seg>.npz
        ▼
[Step 1.4] step14_npz_to_mat
        │  unzip NPZ → read .npy members → permute (nz,ny,nx)→(ny,nx,nz)
        │  cm→mm geometry, parse CT label
        │  writes: dose_*.mat  (variable: raw_field_dose)  → RayStationFiles/<id>/<session>/
        ▼
[Step 1.5] step15_process_doses
        │  discover 2 CBCTs (+RTSTRUCTs) → sort by datetime → CT_1 (earlier), CT_3 (later)
        │  resample each CBCT to dose grid, HU→density, build body/couch/tissue masks
        │  per field: load → validate geometry → zero invalid voxels (per CT label) → save
        │  accumulate total_rs_dose + per-CT totals
        ▼
processed/  (per-field dose_*.mat, CBCT{1,3}_resampled.mat, tissue_masks.mat,
             total_rs_dose.mat, total_dose_CT_{1,3}.mat, metadata.mat)
        ▼
pipeline_simulate.m  (k-Wave forward sim)
```

---

## 3. `pipeline_compress.m` (orchestrator)

A **script** (clears the workspace at top), not a function. Configure at the top, run.

### CONFIG fields

| Field | Default | Meaning |
|---|---|---|
| `CONFIG.patients` | `{'1194203'}` | Cell array of patient IDs to process. |
| `CONFIG.sessions` | `{'Session_1'}` | Cell array of session names. |
| `CONFIG.treatment_site` | `'Pancreas'` | Subfolder under `EthosExports/<id>/`. Used for the fallback path. |
| `CONFIG.working_dir` | `fileparts(mfilename('fullpath'))` | Base dir = the folder this script lives in. Host-agnostic (cluster or laptop) as long as `RayStationFiles`/`EthosExports` are checked out beside the scripts. |
| `CONFIG.apply_dose_masking` | `true` | Zero dose outside body / inside couch. Set `false` only for debugging. |
| `CONFIG.use_sparse_storage` | `true` | Save dose/mask arrays as sparse 2D. Reconstruct 3D with `reshape(full(x), dims)`. |
| `CONFIG.run_step14` | `true` | Run Step 1.4 (NPZ→.mat). |
| `CONFIG.run_step15` | `true` | Run Step 1.5 (process/resample). |
| `CONFIG.skip_completed` | `true` | Resume mode — each step inspects `processed/` and skips already-finished sub-tasks. `false` forces a full re-run. |

### Control flow

1. **Path setup** — adds `working_dir`, `pipeline/` (step functions), and `utils/` (helpers)
   to the MATLAB path via `genpath`.
2. **Patient × session loop** — each iteration wrapped in `try/catch`; a failure records
   `status='error'` and continues to the next session (never aborts the whole run).
3. **Input presence check** (`check_raystation_files`) — verifies field-dose files exist in
   `RayStationFiles` first, else the `EthosExports/.../sct` fallback. If none: sets
   `status='awaiting_raystation'` and skips.
4. **Status report** (`report_processed_status`) — prints which `processed/` outputs already
   exist (purely informational; actual skipping happens inside the step functions).
5. **Step 1.4** → `step14_npz_to_mat` (if `run_step14`).
6. **Step 1.5** → `step15_process_doses` (if `run_step15`).
7. **Summary** (`generate_compress_summary`) — per-session status, field count, max dose.

### Local helper functions (in `pipeline_compress.m`)

| Function | Purpose |
|---|---|
| `check_raystation_files(rs_dir, ethos_dir)` | Returns `(exists, num_files, scan_dir)`. Scans dose-file patterns in priority order (`dose_*.npz`, `dose_*.dcm`, `Plan_Field*`, `RD.*`) at the resolved input directory. |
| `make_result_key(patient_id, session)` | Builds a struct-field-safe key, e.g. `P1194203_Session1`. |
| `init_patient_result(patient_id, session)` | Initializes the per-session RESULTS entry. |
| `report_processed_status(processed_dir, input_dir, skip_completed)` | Prints presence/absence of each expected output + NPZ-vs-.mat counts. Informational only. |
| `present(tf)` | `'present'` / `'MISSING'` label. |
| `generate_compress_summary(results)` | Final per-session rollup. |

---

## 4. Step 1.4 — `step14_npz_to_mat.m`

```matlab
n_converted = step14_npz_to_mat(patient_id, session, config)
```

Bridges the RayStation NPZ export and the MATLAB pipeline. Each `dose_*.npz` in the resolved
input dir is unpacked and written as a per-field `.mat` (same stem, `.mat` extension) into
`RayStationFiles/<id>/<session>/`.

- **Inputs:** `patient_id`, `session` (char/string), `config` (needs `.working_dir`; uses
  `.treatment_site`, `.skip_completed`).
- **Output:** `n_converted` — count of NPZ files successfully converted.
- **Resume:** with `config.skip_completed=true`, an NPZ whose target `.mat` already exists is
  skipped.
- **Robustness:** per-file `try/catch` — a bad NPZ warns and is skipped, not fatal.

### NPZ input schema (from `calc_beam_plan_doses.py`)

| Member | Type / shape | Notes |
|---|---|---|
| `dose` | float32, `(nz, ny, nx)` | **C-order** |
| `voxel_size_cm` | float32, `(3,)` | `[vx, vy, vz]` in cm |
| `corner_cm` | float32, `(3,)` | `[x, y, z]` in cm |
| `nx, ny, nz` | int32 scalars | |

### Output `.mat` schema (variable `raw_field_dose`)

| Field | Type | Meaning |
|---|---|---|
| `.dose_Gy` | double `[ny, nx, nz]` | MATLAB order `(row=Y, col=X, slice=Z)` — permuted from NPZ `(nz,ny,nx)` via `permute(dose, [2 3 1])`. |
| `.origin` | double 3×1, mm | DICOM patient coords (`corner_cm × 10`). |
| `.spacing` | double 3×1, mm | `[dx; dy; dz]` (`voxel_size_cm × 10`). |
| `.dimensions` | double row `[ny, nx, nz]` | |
| `.ct_label` | char | Parsed from filename (`CT_1`, `CT_3`) or `''`. |
| `.source_npz` | char | Original NPZ filename. |

Saved with `-v7.3`.

### Private helpers (self-contained NPY reader)

| Helper | Purpose |
|---|---|
| `local_unpack_npz(npz_path, dest)` | An NPZ is a ZIP of `.npy` members → `unzip`. |
| `local_find_member(extracted, name)` | Find an extracted member by basename. |
| `local_extract_ct_label(npz_name)` | Regex the `CT_k` token from the filename. |
| `local_read_npy(path)` | Minimal NPY 1.0/2.0 reader. Parses magic, version, header; supports `<f4 <f8 <i4 <i8 <u4 <u1 |u1 |i1`. Handles C-order (reshape reversed dims + permute) vs Fortran order. |
| `local_parse_descr/fortran/shape(header)` | Regex the NPY header dict fields. |

> **Why a hand-rolled NPY reader:** avoids requiring a Python bridge or extra toolbox on the
> cluster. It reads exactly the dtypes RayStation emits.

---

## 5. Step 1.5 — `step15_process_doses.m`

```matlab
[field_doses, cbct_resampled, total_rs_dose, metadata] = ...
    step15_process_doses(patient_id, session, config)
```

The heart of the compress stage. ~2300 lines; the main function orchestrates 8 numbered
phases, with ~20 private helpers below it.

### Config fields consumed

| Field | Default | Meaning |
|---|---|---|
| `.working_dir` | *(required)* | Base directory. |
| `.treatment_site` | `'Pancreas'` | For the `EthosExports/.../sct` fallback path. |
| `.apply_dose_masking` | `true` | Zero dose outside body / in couch before saving each field. |
| `.batch_size` | `1000` | Fields per batch; batch subtotal folds into the running total then clears, to cap memory. |
| `.skip_completed` | `true` | Resume — detect existing outputs and skip. |
| `.use_sparse_storage` | `true` | Save dose/mask arrays as sparse 2D `[nRows·nCols, nSlices]`. |

### Outputs

- **`field_doses`** — cell array (1 per input file). Lightweight summary struct per field
  (`filepath`, `beam_num`, `seg_num`, `field_num`, `plan_type`, `gantry_angle`, `meterset`,
  `max_dose_Gy`, `source_file`, `isocenter`, `jaw_x`, `jaw_y`, `body_masked`, `couch_masked`).
  **No** dose array — the full dose lives in the per-field `.mat`.
- **`cbct_resampled`** — `struct('CT_1', CBCT1_resampled, 'CT_3', CBCT3_resampled)`. Each carries
  `cubeHU`, `cubeDensity`, `tissueMask`, `roiNames`, `bodyMask`, `couchMask`, dose-grid
  `origin`/`spacing`/`dimensions`, original CBCT dims/spacing, `series_uid`, `series_datetime`,
  `ct_label`.
- **`total_rs_dose`** — 3D sum of all field doses (Gy), zeroed with the **union** body mask and
  **intersection** couch mask across CT_1/CT_3.
- **`metadata`** — combined geometry (`origin`, `spacing`, `dimensions`), IDs, `num_fields`,
  `beam_metadata` (incl. isocenter + jaws for sensor placement), voxel counts, masking flags.

### Files written to `processed/`

| File | Variable(s) | Notes |
|---|---|---|
| `dose_<id>_<session>_<plan>[_CT_k]_B<beam>_<seg>.mat` | `field_dose` | Per field. Masked before save. Sparse when enabled (`is_sparse`, `dose_dims`). |
| `CBCT1_resampled.mat` / `CBCT3_resampled.mat` | `CBCT1_resampled` / `CBCT3_resampled` | Resampled CBCT + masks. |
| `tissue_masks.mat` | `*_ct1` / `*_ct3` | Both CBCTs' tissue/body/couch masks (sparse `_sp` + `_dims` when enabled). |
| `total_rs_dose.mat` | `total_rs_dose` (or `_sparse`+`_dims`) | Cross-CBCT total. |
| `total_dose_CT_1.mat` / `total_dose_CT_3.mat` | `ct_total` (or `_sparse`+`_dims`) | Per-CT-image totals. |
| `metadata.mat` | `metadata` | |

### The 8 phases

1. **Find field-dose files** — input-format priority: `dose_*.mat` (preferred, from Step 1.4) →
   `dose_*.dcm` → legacy patterns (`Plan_Field*`, `Beam*`, `RD.*`). `input_format` (`'mat'` /
   `'dicom'`) selects the per-field loader branch later.
2. **(1.5) Detect existing outputs** — `detectProcessedOutputs` builds `skip_status` + the
   expected per-field output path list. If **all** outputs exist and `skip_completed`, load and
   return early (`loadExistingStep15Outputs`). `need_total_accum` decides whether the per-field
   loop must re-accumulate totals.
3. **Load RTPLAN** — `loadRtplanMetadata` reads gantry angles, metersets, isocenter, jaws per
   beam (needed downstream by `determine_sensor_mask`).
4. **Reference grid** — first field establishes `(origin, spacing, dims)`; both `.mat` and DICOM
   inputs normalize to the same triple. Initializes `total_rs_dose`, per-CT accumulators, metadata.
5. **Discover + resample CBCTs** — `discoverCbctSeries` finds the two CT series and their
   RTSTRUCTs; sorts by `SeriesDate+SeriesTime` → **CT_1 (earlier), CT_3 (later)**. Each CBCT is
   loaded (`loadCbctImagesFromFiles`), resampled to the dose grid (`resampleSctToDoseGrid`,
   linear interp3, −1000 HU fill), and HU→density (`huToDensity`). Cached CBCTs are reloaded
   when present.
6. **RTSTRUCTs → masks** — `loadRtstructAndCreateMasksFromFile` builds per-CBCT tissue/body/couch
   masks (poly2mask per slice, `fillMaskZGaps` to bridge contour/grid slice misalignment).
   Precomputes per-CBCT invalid-dose masks: `invalid = ~(body & ~couch)`.
7. **Process each field (batched)** — per field: skip if cached; else load → `validateGeometry`
   (resample if dims differ) → `extractBeamInfo` (beam/seg/field/plan from filename) → parse
   `CT_k` label → **select that CBCT's invalid mask** (rejects any label that isn't CT_1/CT_3) →
   `getBeamMetadata` / `getBeamGeometry` → zero invalid voxels → save (sparse) → accumulate into
   batch subtotal + per-CT accumulator. Batch subtotal folds into `total_rs_dose` then clears.
8. **Finalize totals + save** — mask `total_rs_dose` with union-body / intersection-couch; mask
   each per-CT total with its own mask; build `CBCT{1,3}_resampled` structs; save all outputs.
   When totals already exist (`~need_total_accum`), load them instead of recomputing.

### Coordinate & physics notes (carried in the code)

- Array order is `(row=Y, col=X, slice=Z)` throughout; grid dims stored `[ny, nx, nz]`.
- DICOM Z-resolution comes from **`GridFrameOffsetVector`**, never `PixelSpacing`
  (`extractDoseSpacing`).
- HU→density is piecewise-linear (`huToDensity`): air 1.2, lung, soft (`≈1000+HU`), bone, clamped
  to `[1, 7800]` kg/m³.
- **Invalid-dose rule:** a voxel is zeroed if it is outside body **or** inside couch (per the
  field's own CT image). For the cross-CBCT total: body in **either** CBCT and not couch in **both**.

### Private helpers (grouped)

| Group | Helpers |
|---|---|
| **Resume / caching** | `detectProcessedOutputs`, `buildFieldOutputPath`, `loadExistingStep15Outputs`, `loadCachedTissueMasks`, `presentStr` |
| **Geometry** | `extractDoseSpacing`, `validateGeometry`, `resampleDoseToGrid`, `resampleSctToDoseGrid` |
| **Filename parsing** | `extractBeamInfo` (5 patterns, priority order) |
| **RTPLAN metadata** | `loadRtplanMetadata`, `getBeamMetadata`, `getBeamGeometry` |
| **Image loading** | `loadSctImages`, `loadCbctImagesFromFiles`, `discoverCbctSeries`, `discoverCbctSeriesInDir` |
| **Density / masks** | `huToDensity`, `loadRtstructAndCreateMasks`, `loadRtstructAndCreateMasksFromFile`, `fillMaskZGaps` |

> **Only the first function in a `.m` file is externally callable** — all of the above are
> private to Step 1.5. `loadSctImages` and `loadRtstructAndCreateMasks` are retained legacy
> directory-globbing variants; the live path uses the explicit-file variants
> (`loadCbctImagesFromFiles`, `loadRtstructAndCreateMasksFromFile`).

---

## 6. `utils/resolve_input_dir.m`

```matlab
scan_dir = resolve_input_dir(primary_dir, fallback_dir, patterns)
```

Chooses **where to read** inputs: `primary_dir` (`RayStationFiles/<id>/<session>`) if it holds any
file matching a pattern; else `fallback_dir` (`EthosExports/.../<session>/sct`) if IT does; else
`primary_dir` unchanged (so the caller reports its normal "not found"). **Reading only** — outputs
always go to `RayStationFiles/.../processed`. `patterns` is a cellstr (or char) of `dir()` globs,
e.g. `{'dose_*.npz'}`. Local helper `dir_has_any` does the existence check.

---

## 7. Common operations

**Run the whole stage** (edit CONFIG at the top of `pipeline_compress.m`, then run the script).

**Reconstruct a sparse-stored dose:**
```matlab
s = load('dose_1194203_Session_1_adapted_CT_1_B13_0.mat', 'field_dose');
if isfield(s.field_dose,'is_sparse') && s.field_dose.is_sparse
    dose3D = reshape(full(s.field_dose.dose_Gy), s.field_dose.dose_dims);
end
```

**Reconstruct the total dose:**
```matlab
t = load('total_rs_dose.mat');
if isfield(t,'total_rs_dose_sparse')
    total3D = reshape(full(t.total_rs_dose_sparse), t.total_rs_dose_dims);
else
    total3D = t.total_rs_dose;
end
```

**Force a full re-run:** set `CONFIG.skip_completed = false`.
**Debug without masking:** set `CONFIG.apply_dose_masking = false` (dose kept everywhere).

---

## 8. Gotchas

- **`pipeline_compress.m` is a script, not a function** — it `clear`s the workspace. Set CONFIG in
  the file, don't pass args.
- **Every field dose must carry a `CT_1`/`CT_3` label.** Step 1.5 errors on any other label —
  masks are per-CT now, with no default. Legacy DICOM (no label) accumulates into `total_rs_dose`
  but has `ct_label=''` and no per-CT total.
- **CT_1 = earlier CBCT, CT_3 = later CBCT** (by `SeriesDate+SeriesTime`). "CT 2" and "CT 3" are
  synonyms in conversation; the code uses `CT_3`.
- **Per-MATLAB-file visibility:** helpers are private to their file. Don't call
  `discoverCbctSeries` etc. from outside Step 1.5.
- **Do not lint/parse these `.m` files with Python or bash** (project rule). Edit from reading;
  correctness that needs MATLAB should be stated, not script-verified.
