from connect import *
import numpy as np
import json
import os

# ---------------------------------------------------------------
# Grab the DoseValues object from the evaluation in the state tree
# ---------------------------------------------------------------
case = get_current("Case")

dose_eval   = case.TreatmentDelivery.FractionEvaluations[0] \
                  .DoseOnExaminations[0].DoseEvaluations[0]
dose_values = dose_eval.DoseValues

# ---------------------------------------------------------------
# Pull geometry from InDoseGrid
# ---------------------------------------------------------------
grid = dose_values.InDoseGrid

nx = int(grid.NrVoxels.x)
ny = int(grid.NrVoxels.y)
nz = int(grid.NrVoxels.z)

vx = (float(grid.VoxelSize.x),
      float(grid.VoxelSize.y),
      float(grid.VoxelSize.z))

corner = (float(grid.Corner.x),
          float(grid.Corner.y),
          float(grid.Corner.z))

print(f"Grid: {nx} x {ny} x {nz}  voxel={vx} cm  corner={corner} cm")

# ---------------------------------------------------------------
# Pull the flat dose array and reshape (z, y, x)
# ---------------------------------------------------------------
flat = dose_values.DoseData                       # CLR Array[Single], length = nx*ny*nz
expected = nx * ny * nz

dose = np.fromiter(flat, dtype=np.float32, count=expected).reshape(nz, ny, nx)

print(f"Dose array shape  : {dose.shape}")
print(f"Dose min / max    : {dose.min():.4f} / {dose.max():.4f} (Gy, eval units)")

# ---------------------------------------------------------------
# Save: array + sidecar JSON with geometry
# ---------------------------------------------------------------
out_dir = r"\\Client\C$\Users\YourUser\Desktop\dose_export"  # adjust for your VDI mapping
os.makedirs(out_dir, exist_ok=True)

np.savez_compressed(
    os.path.join(out_dir, "evaluation_dose.npz"),
    dose=dose,
    voxel_size_cm=np.array(vx, dtype=np.float32),
    corner_cm=np.array(corner, dtype=np.float32),
    shape_zyx=np.array(dose.shape, dtype=np.int32),
)

with open(os.path.join(out_dir, "evaluation_dose.json"), "w") as f:
    json.dump({
        "shape_zyx":     list(dose.shape),
        "voxel_size_cm": list(vx),
        "corner_cm":     list(corner),
        "patient_id":    get_current("Patient").PatientID,
        "case_name":     case.CaseName,
    }, f, indent=2)

print(f"Saved to: {out_dir}")