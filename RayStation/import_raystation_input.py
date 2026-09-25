from connect import *
import os


# ============================================================
# Configuration
# ============================================================
# Root of the exploded-segment files (Step 0.6 output):
#   Raystation_Input/{patient_id}/{session}/
INPUT_ROOT = "F:/ETHOS_Simulations/Raystation_Input"

# Session folders to import. One new case is created per folder, named after it.
SESSIONS = ["Session_1", "Session_2"]


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

    # Import everything into a new case named after the session folder
    patient.Save()
    warnings = patient.ImportDataFromPath(
        Path=folder, SeriesOrInstances=series, CaseName=session)
    print(f"  Imported. Warnings: {warnings}")

    patient.Save()
    existing_cases.append(session)

print("[IMPORT] Done.")
