# Context

## Background

You are an expert at RayStation and RayStation scripting. Imported to RayStation are the two CBCT files from `pipeline_setup.m`, the exploded plans, the RTSTRUCT files for the CBCTs, and the frame-of-reference registration between the two images.

## Goal

Compute and export a per-beam dose DICOM for each beam in each plan using a RayStation Python script (`calc_beam_plan_doses.py`).

---

# Script: `calc_beam_plan_doses.py`

## Configuration

| Parameter | Value |
|---|---|
| `BASE_EXPORT_ROOT` | `F:/ETHOS_Simulations/RayStationFiles` |
| `EXAMINATION_LABELS` | `["CT 1", "CT 3"]` (planning exam is normally `"CT 2"`) |
| `DOSE_ALGORITHM` | `"PhotonMonteCarlo"` (fallback: `"CCDose"` for Collapsed Cone plans) |
| `DOSE_VOXEL_SIZE` | `{'x': 0.25, 'y': 0.25, 'z': 0.25}` cm |

## Plan Naming Convention

Plans must follow: **`{Session_N} {adapted|reference} {origbeam}`**

- Example: `"Session_1 adapted B13"` → `session="Session_1"`, `plan_type="adapted"`, `origbeam="B13"`
- Patient ID is read from `patient.PatientID`, **not** embedded in the plan name.
- Plans not matching this pattern are skipped with a printed warning.

## Output File Naming

Each segment's dose (one RS beam = one segment) is saved as a compressed NumPy file:

```
dose_{patient_id}_{session}_{plan_type}_{CT_1|CT_3}_{origbeam}_{segment:02d}.npz
```

Written to: `RayStationFiles/{patient_id}/{session}/`. The DICOM-export path is commented out in the script.

## Workflow (per plan)

1. Parse all plans in `case.TreatmentPlans`; group valid ones by session.
2. For each session, read `beam_plan_export_progress.log` (resume: skips plan/CT pairs already DONE).
3. For each beam plan:
   - Set dose grid, update grid structures, `ComputeDose(..., ForceRecompute=True)` on the planning exam (normally **CT 2**).
   - `ComputeDoseOnAdditionalSets` on every exam in `EXAMINATION_LABELS` (`"CT 1"`, `"CT 3"`) that is not the planning exam.
   - `patient.Save()`
   - For each CT: `find_beam_doses` gets this beam set's per-beam doses, then one NPZ per segment is written and the plan/CT is logged DONE.
4. Once per session: export CT 1 / CT 3 images and their RTSTRUCTs as DICOM.

## Finding the right dose (important)

Additional-set doses live **case-wide** in
`case.TreatmentDelivery.FractionEvaluations[*].DoseOnExaminations[*].DoseEvaluations[*]`, with
**one DoseEvaluation per beam set** on each exam. Always match on
`dose_eval.ForBeamSet.BeamSetIdentifier() == beam_set.BeamSetIdentifier()` (older RS: `ForBeamset`).
Never take `DoseEvaluations[0]`: it is the first plan ever computed on that exam, so every plan would
export beam 1's dose. That bug hit Session_2 and was fixed on 2026-10-01. On the planning exam, use
`beam_set.FractionDose.BeamDoses` instead.

## Resume / Progress Log

- Log file: `{export_folder}/beam_plan_export_progress.log`
- Format: `DONE,{plan_name}|{CT label},{timestamp}` followed by `FILE,{filename},{timestamp}` lines; `DONE,{session}|CT_RTSTRUCT,...` for the image export.
- On re-run, completed plan/CT pairs are skipped. **Delete the log** (and the stale NPZs) to force a re-export.

## Key RayStation API Calls

```python
patient  = get_current("Patient")
case     = get_current("Case")
beam_set = plan.BeamSets[0]
beam_set.SetDefaultDoseGrid(VoxelSize={'x':0.25,'y':0.25,'z':0.25})
beam_set.FractionDose.UpdateDoseGridStructures()
beam_set.ComputeDose(ComputeBeamDoses=True, DoseAlgorithm="PhotonMonteCarlo", ForceRecompute=True)
patient.Save()
case.ScriptableDicomExport(
    ExportFolderPath=export_folder,
    BeamSets=[beam_set.BeamSetIdentifier()],
    PhysicalBeamDosesForBeamSets=[beam_set.BeamSetIdentifier()],
    IgnorePreConditionWarnings=True
)
```

## Linac / Physics

- Linac: **Halcyon SN1293**, 6FFF photons, X/Y jaws, dual-leaf MLC.
- Dose type exported: **physical beam doses** (one RD*.dcm per beam, per plan).
