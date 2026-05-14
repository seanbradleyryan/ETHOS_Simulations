from connect import *
from collections import defaultdict
from datetime import datetime
import glob
import os
import re


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

# Plan name formats handled by this script:
#   Original (template) plans:   "{session} {adapted|reference} {origbeam}"
#       e.g.  "Session_1 adapted B13"
#   New CT-specific plans:       "{session} {adapted|reference} {ct_label} {origbeam}"
#       e.g.  "Session_1 adapted CT 1 B13"
#
# Group 1 (session)   – Session_<digits>
# Group 2 (plan_type) – "adapted" or "reference"
# Group 3 (ct_label)  – optional, "CT <digits>" (None for original plans)
# Group 4 (origbeam)  – beam label, e.g. "B13"
PLAN_NAME_RE = re.compile(
    r'^(Session_\d+)\s+(adapted|reference)(?:\s+(CT\s+\d+))?\s+(.+)$',
    re.IGNORECASE
)


# ============================================================
# Helper: parse a plan name into its components
# ============================================================
def parse_plan_name(plan_name):
    """
    Parse a plan name. Returns (session, plan_type, ct_label, origbeam),
    where ct_label is None for original template plans and a string
    (e.g. "CT 1") for CT-specific plans. Returns None on no match.
    """
    m = PLAN_NAME_RE.match(plan_name)
    if not m:
        return None
    session   = m.group(1)
    plan_type = m.group(2).lower()
    ct_label  = m.group(3)  # may be None
    origbeam  = m.group(4)
    if ct_label is not None:
        # Normalize internal whitespace: "CT   1" -> "CT 1"
        ct_label = re.sub(r'\s+', ' ', ct_label).strip()
    return session, plan_type, ct_label, origbeam


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


# ============================================================
# Helper: snapshot / diff RD files for rename-after-export
# ============================================================
def snapshot_rd_files(folder):
    return set(
        os.path.basename(p)
        for p in glob.glob(os.path.join(folder, "RD*.dcm"))
    )


def find_new_rd_files(folder, pre_snapshot):
    current   = glob.glob(os.path.join(folder, "RD*.dcm"))
    new_files = [p for p in current
                 if os.path.basename(p) not in pre_snapshot]
    new_files.sort(key=os.path.getmtime)
    return new_files


# ============================================================
# Helper: rename exported RD files to the project naming scheme
# ============================================================
def rename_beam_exports(new_files, export_folder,
                        patient_id, session, plan_type, ct_label,
                        origbeam, beam_set):
    """
    Rename each newly-exported RD*.dcm to:
        dose_{id}_{session}_{plan_type}_{ct_label}_{origbeam}_{segment}.dcm

    {ct_label} has spaces replaced by underscores (e.g. "CT 1" -> "CT_1").
    {segment} is the zero-based beam index within the plan's beam set,
    zero-padded to two digits.
    """
    beams = list(beam_set.Beams)

    if len(new_files) != len(beams):
        print(f"  WARNING: Expected {len(beams)} beam dose file(s), "
              f"found {len(new_files)}. Segment numbering may be mismatched.")

    safe_id       = re.sub(r'[\\/:*?"<>| ]', '_', patient_id)
    safe_session  = re.sub(r'[\\/:*?"<>| ]', '_', session)
    safe_type     = re.sub(r'[\\/:*?"<>| ]', '_', plan_type)
    safe_ct       = re.sub(r'[\\/:*?"<>| ]', '_', ct_label)
    safe_beam     = re.sub(r'[\\/:*?"<>| ]', '_', origbeam)

    final_paths = []
    for segment_idx, src_path in enumerate(new_files):
        desired_name = (
            f"dose_{safe_id}_{safe_session}_{safe_type}_{safe_ct}_"
            f"{safe_beam}_{segment_idx:02d}.dcm"
        )
        desired_path = os.path.join(export_folder, desired_name)

        if os.path.exists(desired_path):
            os.remove(desired_path)

        os.rename(src_path, desired_path)

        rs_beam = beams[segment_idx].Name if segment_idx < len(beams) else "unknown"
        print(f"    segment {segment_idx:02d}  ({rs_beam})  ->  {desired_name}")
        final_paths.append(desired_path)

    return final_paths


# ============================================================
# Resume helpers (key on NEW CT-specific plan names)
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
        f.write(f"# dose_{{id}}_{{session}}_{{plan_type}}_{{ct_label}}_{{origbeam}}_{{segment}}.dcm\n")


# ============================================================
# Helper: extract beam-set configuration from a source beam set
# ============================================================
def get_beam_set_kwargs(src_beam_set):
    """
    Pull machine/modality/technique/etc. from an existing beam set so a new
    beam set can be created with matching configuration. Missing attributes
    are silently omitted (RayStation versions differ slightly).
    """
    kwargs = {}
    try:
        kwargs['MachineName'] = src_beam_set.MachineReference.MachineName
    except Exception:
        pass
    for attr in ('Modality', 'TreatmentTechnique', 'PatientPosition'):
        try:
            val = getattr(src_beam_set, attr, None)
            if val is not None:
                kwargs[attr] = val
        except Exception:
            pass
    try:
        kwargs['NumberOfFractions'] = src_beam_set.FractionationPattern.NumberOfFractions
    except Exception:
        try:
            kwargs['NumberOfFractions'] = src_beam_set.NumberOfFractions
        except Exception:
            pass
    return kwargs


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
        sess, plan_type, ct_label, origbeam = parsed
        if ct_label is not None:
            # CT-specific plan: not a template; skip silently.
            continue
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
            orig_plan_name = plan.Name

            if not plan.BeamSets:
                print(f"  WARNING: '{orig_plan_name}' has no beam sets. Skipping.")
                continue

            original_beam_set = plan.BeamSets[0]

            # --------------------------------------------------------
            # For each target examination, create a CT-specific copy
            # and compute/export dose on that copy.
            # --------------------------------------------------------
            for exam_label in EXAMINATION_LABELS:
                new_plan_name     = f"{sess} {plan_type} {exam_label} {origbeam}"
                new_beam_set_name = new_plan_name  # keep matched for clarity

                if new_plan_name in completed_plans:
                    print(f"\n  SKIPPING '{new_plan_name}' (already exported)")
                    continue

                print(f"\n  --- template='{orig_plan_name}'  ->  '{new_plan_name}' ---")

                # ---- Create new plan on the target examination ------
                existing_plan = None
                for p in case.TreatmentPlans:
                    if p.Name == new_plan_name:
                        existing_plan = p
                        break

                if existing_plan is not None:
                    print(f"    Reusing existing plan '{new_plan_name}'.")
                    new_plan = existing_plan
                else:
                    print(f"    Creating plan on examination '{exam_label}' ...")
                    new_plan = case.AddNewPlan(
                        PlanName=new_plan_name,
                        ExaminationName=exam_label,
                        PlannedBy="",
                        Comment=f"Auto-generated from '{orig_plan_name}'",
                        AllowDuplicateNames=False
                    )

                # ---- Add beam set on the target examination ---------
                if new_plan.BeamSets and len(list(new_plan.BeamSets)) > 0:
                    new_bs = new_plan.BeamSets[0]
                    print(f"    Reusing existing beam set on '{new_plan_name}'.")
                else:
                    bs_kwargs = get_beam_set_kwargs(original_beam_set)
                    print(f"    Adding beam set on '{exam_label}' "
                          f"(kwargs={list(bs_kwargs.keys())}) ...")
                    new_bs = new_plan.AddNewBeamSet(
                        Name=new_beam_set_name,
                        ExaminationName=exam_label,
                        CreateSetupBeams=False,
                        **bs_kwargs
                    )

                    # ---- Copy beams from the original beam set ------
                    print(f"    Copying beams from original beam set ...")
                    new_bs.CopyBeamsFromBeamSet(BeamSetToCopyFrom=original_beam_set)

                examination = case.Examinations[exam_label]
                # case.PatientModel.RegionsOfInterest['Body'].CreateExternalGeometry(Examination=examination, ThresholdLevel=-250)

                # ---- Set default dose grid --------------------------
                vx = DOSE_VOXEL_SIZE
                print(f"    Setting dose grid: {vx['x']} x {vx['y']} x {vx['z']} cm ...")
                new_bs.SetDefaultDoseGrid(VoxelSize=DOSE_VOXEL_SIZE)
                new_bs.FractionDose.UpdateDoseGridStructures()

                # ---- Compute dose (per-beam) ------------------------
                print(f"    Computing dose ({DOSE_ALGORITHM}) ...")
                new_bs.ComputeDose(
                    ComputeBeamDoses=True,
                    DoseAlgorithm=DOSE_ALGORITHM,
                    ForceRecompute=True
                )
                print(f"    Dose computation complete.")

                # ---- Save before export -----------------------------
                patient.Save()
                print(f"    Patient saved.")

                # ---- Snapshot, export, rename -----------------------
                beam_set_id  = new_bs.BeamSetIdentifier()
                pre_snapshot = snapshot_rd_files(export_folder)

                print(f"    Exporting per-beam doses ...")
                case.ScriptableDicomExport(
                    ExportFolderPath=export_folder,
                    BeamSets=[beam_set_id],
                    PhysicalBeamDosesForBeamSets=[beam_set_id],
                    IgnorePreConditionWarnings=True
                )

                new_files = find_new_rd_files(export_folder, pre_snapshot)
                print(f"    {len(new_files)} new RD file(s) detected.")

                if not new_files:
                    print("    WARNING: No new RD files found. Export may have failed.")
                    continue

                final_paths = rename_beam_exports(
                    new_files, export_folder,
                    patient_id, sess, plan_type, exam_label, origbeam, new_bs
                )

                log_plan_completion(progress_log, new_plan_name, final_paths)
                print(f"    Logged completion for '{new_plan_name}'.")

    # --------------------------------------------------------
    print(f"\n{'='*60}")
    print("All sessions processed successfully!")
    print(f"{'='*60}")

except Exception as e:
    import traceback
    print(f"\nERROR: {str(e)}")
    print(traceback.format_exc())
