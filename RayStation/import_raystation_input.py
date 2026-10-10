from connect import *
import os


# ============================================================
# Configuration
# ============================================================
# Root of the exploded-segment files (Step 0.6 output):
#   Raystation_Input/{patient_id}/{session}/
INPUT_ROOT = "F:/ETHOS_Simulations/Raystation_Input"

# Session folders to import. One new case is created per folder, named after it.
SESSIONS = ["Session_1", "Session_2", "Session_3", "Session_4", "Session_5", "Session_6", "Session_7", "Session_8", "Session_9", "Session_10"]


# ============================================================
# Main
# ============================================================
patient    = get_current("Patient")
patient_db = get_current("PatientDB")
patient_id = patient.PatientID

existing_cases = [c.CaseName for c in patient.Cases]

for session in SESSIONS:
    folder = os.path.join(INPUT_ROOT, patient_id, session)
    print(f"[IMPORT] {session}: {folder}")

    if not os.path.isdir(folder):
        print(f"  WARNING: folder not found, skipping")
        continue
    if session in existing_cases:
        print(f"  Case '{session}' already exists, skipping")
        continue

    # Find every series in the folder belonging to the current patient
    studies = patient_db.QueryStudiesFromPath(
        Path=folder, SearchCriterias={'PatientID': patient_id})
    series = []
    for study in studies:
        series += patient_db.QuerySeriesFromPath(Path=folder, SearchCriterias=study)
    print(f"  Found {len(studies)} studies, {len(series)} series")

    if not series:
        print(f"  WARNING: nothing to import for patient {patient_id}, skipping")
        continue

    # CaseName=None makes RayStation create a new case (CaseName is the
    # target of an existing case). Rename the new case afterwards.
    patient.Save()
    cases_before = [c.CaseName for c in patient.Cases]
    warnings = patient.ImportDataFromPath(
        Path=folder, SeriesOrInstances=series, CaseName=None)
    print(f"  Imported. Warnings: {warnings}")

    new_cases = [c for c in patient.Cases if c.CaseName not in cases_before]
    if len(new_cases) != 1:
        print(f"  WARNING: expected 1 new case, found {len(new_cases)}; not renaming")
    else:
        print(f"  Renaming new case '{new_cases[0].CaseName}' -> '{session}'")
        new_cases[0].CaseName = session

    patient.Save()
    existing_cases.append(session)

print("[IMPORT] Done.")
