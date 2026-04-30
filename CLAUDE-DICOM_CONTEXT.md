# DICOM_CONTEXT.md — DICOM Chains, MLC Correction & Segment Explosion

> Referenced when working on: `step0_sort_dicom.m`, `step05_fix_mlc_gaps.m`, `step06_explode_segments.m`

## DICOM Reference Chain (Step 0)

Files are matched via the following chain:

```
CT SeriesInstanceUID
  → RTSTRUCT.ReferencedFrameOfReferenceSequence
    → RTPLAN.ReferencedStructureSetSequence
      → RTDOSE.ReferencedRTPlanSequence
```

Step 0 walks this chain to associate each SCT series with its corresponding RT files.

## MLC Gap Correction — Halcyon (Step 0.5)

- Halcyon uses a **dual-layer MLC** (proximal + distal).
- Minimum enforced gap: **0.5 mm**, with **0.4 mm expansion per side**.
- Valid leaf range: **[−140, 140] mm**.

## Beam Segment Explosion (Step 0.6)

`step06_explode_segments` reads the MLC-corrected RTPLAN and writes one file per beam:

- Output filename pattern: `RTPLAN_{patient}_{session}_{reference|adapted}_B<N>.dcm`
- Output location: `Raystation_Input/[PatientID]/[Session]/`
- Each file is a single-segment 2-CP plan; imported into RayStation for per-segment dose calculation.

Output struct fields:

| Field | Content |
|-------|---------|
| `.reference` | Cell array of file paths for reference plan beams |
| `.adapted` | Cell array of file paths for adapted plan beams |
| `.all` | Cell array of all generated file paths |

## Gotchas

- **`dicomCollection()` fragility:** Can fail with malformed DICOM headers. Fall back to manual `dicominfo()` iteration if `dicomCollection` errors out.
