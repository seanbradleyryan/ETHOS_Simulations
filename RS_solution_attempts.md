# RS_solution_attempts.md

**Goal:** `calc_beam_plan_dose.py` — produce and export a dose from two identical plans copied from an external source using RayStation's dose calculations.

---

## Attempt 1: Copy beamset via alternative plan

**Approach:** Create an alternative plan and copy the beamset over from the imported plan.

**Failure:** The action is not scriptable in RayStation. When recording the action via the RayStation script recorder, it throws an error. The action works interactively in the GUI but cannot be automated.

---

## Attempt 2: ComputeDoseOnAdditionalSets

**Approach:** Use `ComputeDoseOnAdditionalSets` to compute dose on the copied plan.

**Failure:** This creates an *evaluation dose*, not a *nominal dose*. Evaluation doses are not exportable from RayStation.

---

## Attempt 3: Direct export of evaluation dose

**Approach:** Access the evaluation dose object directly and export it manually via script.

**Failure:** Performance is prohibitive — approximately 10 seconds per dose, with ~6000 doses to process. Not a viable path at scale.

---

## Attempt 4: New plan on a new image set with beamset copy

**Approach:** Create a new plan on a different image set and copy the beamset over from the source plan.

**Failure:** RayStation throws an error preventing the image set from being changed on this type of plan. The plan type is locked to its original image set.

---

## Open Questions / Next Directions

- Is there a RayStation API method that computes and commits a nominal dose directly onto a copied or newly created plan?
- Can a plan be serialized (e.g., exported as RTPLAN DICOM) and re-imported to a fresh case where dose can be computed normally?
- Is there a way to promote an evaluation dose to a nominal dose via script?
- Are there RayStation version-specific APIs (e.g., `CreateDoseSpecificationPoint`, `SetDoseGrid`, `ComputeDose`) that apply here?
