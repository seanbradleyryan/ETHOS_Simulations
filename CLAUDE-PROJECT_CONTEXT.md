## Dose–CT Coupling (do not treat anatomy and dose as separable)

The dose grid is inherently bound to the CT geometry it was computed on: Dose1 belongs to CT1, Dose2 to CT2. There is no physically valid state with one dose grid on a different CT (e.g., "Dose3 on CT1" is meaningless). Across the two timepoints, anatomy and dose change together by construction — a geometry change IS a dosimetric change, because the beam deposits into a different configuration.

Acoustic waves propagate fast relative to anatomical motion, so the medium is frozen during a single delivery; there is no intra-delivery anatomy change to model.

Consequence: a reconstruction difference driven by anatomy is NOT a false positive. It is the true positive — the system correctly reporting that the patient was dosed differently. Do not frame anatomy-driven reconstruction differences as artifacts to subtract, and do not propose studies that decouple dose from anatomy (e.g., same dose on two CTs) except as an attribution diagnostic, which is not a current priority.

True vs. false positive criterion (deposition vs. sensing path):
- Change at/upstream of dose deposition (target moves, organ fills, beam off-target) → real dose change → TRUE POSITIVE. Not separable from anatomy; do not try to separate.
- Change in the sensing path only, downstream of deposition (air bubble at the transducer face, coupling change, gas pocket between source and array but outside the beam) → no dose change → FALSE POSITIVE → the only class worth a specificity/robustness test.

These are different perturbations on different operators: detection sensitivity (MDD) perturbs the source (p0/dose); the sensing-path robustness test perturbs the medium (acoustic properties along the path).