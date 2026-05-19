from connect import *
from collections import defaultdict
from datetime import datetime
import glob
import os
import re
import numpy as np


# ============================================================
# Configuration
# ============================================================
BASE_EXPORT_ROOT = "F:/RayStationFiles"
BASE_EXPORT_ROOT = r"//Client/C$/Users/80030361/Documents/RayStationFiles"

DOSE_ALGORITHM  = "PhotonMonteCarlo"   # change to "CCDose" if plan used Collapsed Cone
DOSE_VOXEL_SIZE = {'x': 0.25, 'y': 0.25, 'z': 0.25}

# Examinations on which CT-specific copies of each original plan are created
# and on which dose is actually computed/exported.
EXAMINATION_LABELS = ["CT 1", "CT 3"]

# Plan name format:  "{session} {adapted|reference} {origbeam}"
#   e.g.  "Session_1 adapted B13" -> session="Session_1", plan_type="adapted", origbeam="B13"
PLAN_NAME_RE = re.compile(
    r'^(Session_\d+)\s+(adapted|reference)\s+(.+)$',
    re.IGNORECASE
)


# ============================================================
# Helper: parse a plan name into (session, plan_type, origbeam)
# ============================================================
def parse_plan_name(plan_name):
    m = PLAN_NAME_RE.match(plan_name)
    if not m:
        return None
    return m.group(1), m.group(2).lower(), m.group(3)


# ============================================================
# Helper: build (and create) the export directory
# ============================================================
def build_export_dir(patient_id, session):
    try:
        if not patient_id or not session:
            raise ValueError(
                f"Cannot build export path: id={patient_id!r}, session={session!r}"
            )
        path = os.path.join(BASE_EXPORT_ROOT, patient_id, session)
        os.makedirs(path, exist_ok=True)
        return path
    except Exception as exc:
        print(f"  WARNING: Could not build session export directory: {exc}")
        print(f"  WARNING: Falling back to base root: {BASE_EXPORT_ROOT}")
        os.makedirs(BASE_EXPORT_ROOT, exist_ok=True)
        return BASE_EXPORT_ROOT


# --- DICOM export (disabled) ---
# # ============================================================
# # Helper: snapshot / diff RD files for rename-after-export
# # ============================================================
# def snapshot_rd_files(folder):
#     return set(
#         os.path.basename(p)
#         for p in glob.glob(os.path.join(folder, "RD*.dcm"))
#     )
#
#
# def find_new_rd_files(folder, pre_snapshot):
#     current   = glob.glob(os.path.join(folder, "RD*.dcm"))
#     new_files = [p for p in current
#                  if os.path.basename(p) not in pre_snapshot]
#     new_files.sort(key=os.path.getmtime)
#     return new_files
#
#
# # ============================================================
# # Helper: rename exported RD files to the project naming scheme
# # ============================================================
# def rename_beam_exports(new_files, export_folder,
#                         patient_id, session, plan_type, ct_label,
#                         origbeam, beam_set):
#     """
#     Rename each newly-exported RD*.dcm to:
#         dose_{id}_{session}_{plan_type}_{ct_label}_{origbeam}_{segment}.dcm
#
#     {ct_label} has spaces replaced by underscores (e.g. "CT 1" -> "CT_1").
#     {segment} is the zero-based beam index within the plan's beam set,
#     zero-padded to two digits.
#     """
#     beams = list(beam_set.Beams)
#
#     if len(new_files) != len(beams):
#         print(f"  WARNING: Expected {len(beams)} beam dose file(s), "
#               f"found {len(new_files)}. Segment numbering may be mismatched.")
#
#     safe_id       = re.sub(r'[\\/:*?"<>| ]', '_', patient_id)
#     safe_session  = re.sub(r'[\\/:*?"<>| ]', '_', session)
#     safe_type     = re.sub(r'[\\/:*?"<>| ]', '_', plan_type)
#     safe_ct       = re.sub(r'[\\/:*?"<>| ]', '_', ct_label)
#     safe_beam     = re.sub(r'[\\/:*?"<>| ]', '_', origbeam)
#
#     final_paths = []
#     for segment_idx, src_path in enumerate(new_files):
#         desired_name = (
#             f"dose_{safe_id}_{safe_session}_{safe_type}_{safe_ct}_"
#             f"{safe_beam}_{segment_idx:02d}.dcm"
#         )
#         desired_path = os.path.join(export_folder, desired_name)
#
#         if os.path.exists(desired_path):
#             os.remove(desired_path)
#
#         os.rename(src_path, desired_path)
#
#         rs_beam = beams[segment_idx].Name if segment_idx < len(beams) else "unknown"
#         print(f"    segment {segment_idx:02d}  ({rs_beam})  ->  {desired_name}")
#         final_paths.append(desired_path)
#
#     return final_paths
# ---


# ============================================================
# Resume helpers (log key is "{plan_name}|{ct_label}")
# ============================================================
def get_completed_plans(log_path):
    completed = set()
    if os.path.exists(log_path):
        with open(log_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split(',', 2)
                if len(parts) >= 2 and parts[0] == 'DONE':
                    completed.add(parts[1])
    return completed


def log_plan_completion(log_path, plan_name, final_paths):
    with open(log_path, 'a') as f:
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f.write(f"DONE,{plan_name},{ts}\n")
        for fp in final_paths:
            f.write(f"  FILE,{os.path.basename(fp)},{ts}\n")


def init_log(log_path, patient_id, session):
    with open(log_path, 'w') as f:
        f.write(f"# Beam plan dose export log - started {datetime.now()}\n")
        f.write(f"# Patient: {patient_id}  |  Session: {session}\n")
        f.write(f"# dose_{{id}}_{{session}}_{{plan_type}}_{{ct_label}}_{{origbeam}}_{{segment}}.npz\n")




# ============================================================
# Main
# ============================================================
try:
    patient = get_current("Patient")
    case    = get_current("Case")

    patient_id = patient.PatientID

    # --------------------------------------------------------
    # Parse every plan name; group ORIGINAL (template) plans by session.
    # CT-specific plans (ct_label is not None) are ignored here — they are
    # generated outputs, not templates.
    # --------------------------------------------------------
    # sessions: sess -> [(plan_type, origbeam, plan), ...]
    sessions = defaultdict(list)
    skipped  = []

    for plan in case.TreatmentPlans:
        parsed = parse_plan_name(plan.Name)
        if parsed is None:
            skipped.append(plan.Name)
            continue
        sess, plan_type, origbeam = parsed
        sessions[sess].append((plan_type, origbeam, plan))

    if skipped:
        print("Plans skipped (do not match expected naming convention):")
        for s in skipped:
            print(f"  - '{s}'")

    if not sessions:
        raise RuntimeError(
            "No valid template plans found. Ensure plans follow "
            "'{Session_N} {adapted|reference} {origbeam}',  e.g. 'Session_1 adapted B13'."
        )

    for key in sessions:
        sessions[key].sort(key=lambda x: (x[0], x[1]))

    total_plans = sum(len(v) for v in sessions.values())
    print(f"\nFound {len(sessions)} session(s) across {total_plans} template plan(s):")
    print(f"  Patient ID: {patient_id}")
    print(f"  Examinations: {EXAMINATION_LABELS}\n")
    for sess, beam_plans in sorted(sessions.items()):
        print(f"  [{sess}]  ({len(beam_plans)} template plan(s))")
        for plan_type, origbeam, plan in beam_plans:
            bs      = plan.BeamSets[0] if plan.BeamSets else None
            n_beams = len(bs.Beams) if bs else 0
            print(f"    {plan_type:<10s}  origbeam={origbeam:<10s}  plan='{plan.Name}'  beams={n_beams}")

    # ============================================================
    # OUTER LOOP – sessions
    # ============================================================
    for sess, beam_plans in sorted(sessions.items()):

        print(f"\n{'='*60}")
        print(f"SESSION: {patient_id} / {sess}  ({len(beam_plans)} template plan(s))")
        print(f"{'='*60}")

        export_folder   = build_export_dir(patient_id, sess)
        print(f"  Export folder  : {export_folder}")

        progress_log    = os.path.join(export_folder, "beam_plan_export_progress.log")
        completed_plans = get_completed_plans(progress_log)

        if completed_plans:
            print(f"  Resuming – {len(completed_plans)} plan(s) already done: "
                  + ", ".join(sorted(completed_plans)))
        else:
            init_log(progress_log, patient_id, sess)
            print("  No prior exports for this session. Starting fresh.")

        # ========================================================
        # INNER LOOP – template plans
        # ========================================================
        for plan_type, origbeam, plan in beam_plans:
            plan_name = plan.Name

            if not plan.BeamSets:
                print(f"  WARNING: '{plan_name}' has no beam sets. Skipping.")
                continue

            beam_set    = plan.BeamSets[0]
            beam_set_id = beam_set.BeamSetIdentifier()

            print(f"\n  --- {plan_type}  origbeam='{origbeam}'  plan='{plan_name}' ---")
            print(f"    BeamSet   : {beam_set_id}")
            print(f"    Beams ({len(beam_set.Beams)}):")
            for bi, b in enumerate(beam_set.Beams):
                print(f"      [{bi:02d}] {b.Name}  MU={b.BeamMU:.4f}")

            # ---- Determine primary examination ----------------------
            # API name varies across RayStation versions; try the common
            # accessors and fall back to None on failure.
            primary_exam = None
            try:
                primary_exam = beam_set.GetPlanningExamination().Name
            except Exception:
                try:
                    primary_exam = plan.GetStructureSet().OnExamination.Name
                except Exception:
                    pass
            print(f"    Primary examination: {primary_exam!r}")

            # ---- Set default dose grid (once) -----------------------
            vx = DOSE_VOXEL_SIZE
            print(f"    Setting dose grid: {vx['x']} x {vx['y']} x {vx['z']} cm ...")
            beam_set.SetDefaultDoseGrid(VoxelSize=DOSE_VOXEL_SIZE)
            beam_set.FractionDose.UpdateDoseGridStructures()

            # ---- Compute primary dose (per-beam) --------------------
            print(f"    Computing primary dose ({DOSE_ALGORITHM}) ...")
            beam_set.ComputeDose(
                ComputeBeamDoses=True,
                DoseAlgorithm=DOSE_ALGORITHM,
                ForceRecompute=True
            )
            print(f"    Primary dose computation complete.")

            # ---- Compute dose on additional examinations ------------
            additional_exam_labels = [e for e in EXAMINATION_LABELS if e != primary_exam]
            if additional_exam_labels:
                print(f"    Computing dose on additional sets: {additional_exam_labels} ...")
                additional_exam_objects = [
                    e for e in case.Examinations
                    if e.Name in additional_exam_labels
                ]
                beam_set.ComputeDoseOnAdditionalSets(
                    Examinations=additional_exam_objects,
                    ComputeBeamDoses=True,
                    AllowGridExpansion=True,
                )
                print(f"    Additional-set dose computation complete.")

            patient.Save()
            print(f"    Patient saved.")

            # ---- Export per-beam dose for each CT in EXAMINATION_LABELS
            for ct_label in EXAMINATION_LABELS:
                log_key = f"{plan_name}|{ct_label}"

                if log_key in completed_plans:
                    print(f"    SKIPPING '{log_key}' (already exported)")
                    continue

                # --- DICOM export (disabled) ---
                # pre_snapshot = snapshot_rd_files(export_folder)
                #
                # print(f"    Exporting per-beam doses for '{ct_label}' ...")
                # # TODO: exact kwarg for per-beam dose on an additional set
                # # varies by RayStation version. If ct_label == primary_exam,
                # # PhysicalBeamDosesForBeamSets works. Otherwise, the additional-
                # # set form is needed; common candidates:
                # #   PhysicalBeamDosesForBeamSetDoseOnAdditionalSet=
                # #       [{'BeamSetId': beam_set_id, 'ExaminationName': ct_label}]
                # #   AdditionalBeamSetDoseOnAdditionalSets=[...]
                # # Adjust below based on the live RayStation traceback.
                # if ct_label == primary_exam or primary_exam is None:
                #     case.ScriptableDicomExport(
                #         ExportFolderPath=export_folder,
                #         BeamSets=[beam_set_id],
                #         PhysicalBeamDosesForBeamSets=[beam_set_id],
                #         IgnorePreConditionWarnings=True
                #     )
                # else:
                #     case.ScriptableDicomExport(
                #         ExportFolderPath=export_folder,
                #         BeamSets=[beam_set_id],
                #         PhysicalBeamDosesForBeamSetDoseOnAdditionalSet=[
                #             {'BeamSetId': beam_set_id, 'ExaminationName': ct_label}
                #         ],
                #         IgnorePreConditionWarnings=True
                #     )
                #
                # new_files = find_new_rd_files(export_folder, pre_snapshot)
                # print(f"    {len(new_files)} new RD file(s) detected.")
                #
                # if not new_files:
                #     print(f"    WARNING: No new RD files found for '{ct_label}'.")
                #     continue
                #
                # final_paths = rename_beam_exports(
                #     new_files, export_folder,
                #     patient_id, sess, plan_type, ct_label, origbeam, beam_set
                # )
                # ---

                # Direct extraction: find FractionEvaluation whose DoseOnExamination
                # matches ct_label, then pull per-beam doses directly.
                print(f"    Extracting per-beam doses for '{ct_label}' ...")
                fe_idx  = None
                doe_idx = None
                for fi, fe in enumerate(case.TreatmentDelivery.FractionEvaluations):
                    for di, doe in enumerate(fe.DoseOnExaminations):
                        if doe.OnExamination.Name == ct_label:
                            fe_idx  = fi
                            doe_idx = di
                            break
                    if fe_idx is not None:
                        break

                if fe_idx is None:
                    print(f"    WARNING: No FractionEvaluation found for '{ct_label}'. Skipping.")
                    continue

                dose_on_exam = (case.TreatmentDelivery
                                    .FractionEvaluations[fe_idx]
                                    .DoseOnExaminations[doe_idx])
                beam_doses = dose_on_exam.DoseEvaluations[0].BeamDoses

                safe_id      = re.sub(r'[\\/:*?"<>| ]', '_', patient_id)
                safe_session = re.sub(r'[\\/:*?"<>| ]', '_', sess)
                safe_type    = re.sub(r'[\\/:*?"<>| ]', '_', plan_type)
                safe_ct      = re.sub(r'[\\/:*?"<>| ]', '_', ct_label)
                safe_beam    = re.sub(r'[\\/:*?"<>| ]', '_', origbeam)

                final_paths = []
                for j in range(len(beam_doses)):
                    grid = beam_doses[j].InDoseGrid
                    nx   = int(grid.NrVoxels.x)
                    ny   = int(grid.NrVoxels.y)
                    nz   = int(grid.NrVoxels.z)
                    vx_array     = np.array([float(grid.VoxelSize.x),
                                             float(grid.VoxelSize.y),
                                             float(grid.VoxelSize.z)], dtype=np.float32)
                    corner_array = np.array([float(grid.Corner.x),
                                             float(grid.Corner.y),
                                             float(grid.Corner.z)], dtype=np.float32)

                    flat       = beam_doses[j].DoseValues.DoseData
                    dose_array = np.fromiter(flat, dtype=np.float32,
                                             count=nx * ny * nz).reshape(nz, ny, nx)

                    desired_name = (
                        f"dose_{safe_id}_{safe_session}_{safe_type}_{safe_ct}_"
                        f"{safe_beam}_{j:02d}.npz"
                    )
                    save_path = os.path.join(export_folder, desired_name)

                    np.savez_compressed(
                        save_path,
                        dose=dose_array,
                        voxel_size_cm=vx_array,
                        corner_cm=corner_array,
                        nx=np.int32(nx),
                        ny=np.int32(ny),
                        nz=np.int32(nz),
                    )
                    print(f"    Saved: {save_path}")
                    final_paths.append(save_path)

                log_plan_completion(progress_log, log_key, final_paths)
                print(f"    Logged completion for '{log_key}'.")

    # --------------------------------------------------------
    print(f"\n{'='*60}")
    print("All sessions processed successfully!")
    print(f"{'='*60}")

except Exception as e:
    import traceback
    print(f"\nERROR: {str(e)}")
    print(traceback.format_exc())
