# ETHOS IRAI Pipeline — Academic Audit Checklist

Audit of the live path: RayStation export → Steps 1.4/1.5 → `pipeline_simulate` → `run_single_field_simulation` → Step 2.5 → noise null → `study_pass_rates_allsegments`.
Items marked **(verify)** could not be confirmed without MATLAB / RayStation / data.

**Bottom line:** the CT_1-vs-CT_3 difference is real in the code, but part of it may be built in rather than measured:
(1) truth information enters the recons (mask + own-truth gain), (2) recon1 has zero model error while recon3 always carries medium mismatch, (3) error bars treat correlated segments as independent.

---

## Suggested order

- [ ] 1. Cheap geometry checks: 1.1, 1.2/1.3, 1.4, 1.5
- [ ] 2. Rerun Steps 2 + 2.5 with truth mask off and recon3 scaled by CT_1 gain (change hash by hand — 1.9; set `metrics_overwrite = true` — 1.8)
- [ ] 3. Oracle CT_3 run (§3) and CT_1 model-error run (2.2)
- [ ] 4. Revisit null + statistics (2.3, 2.4), noise model (2.5), Grüneisen table (2.6)

---

## 1. Coding errors / inconsistencies

### Could be affecting current results

- [ ] **1.1 RayStation may export the wrong plan's dose (verify)** — `RayStation/calc_beam_plan_doses.py:390-406`
  Takes the first dose record on "CT 1"/"CT 3" and reads `DoseEvaluations[0].BeamDoses` without checking it belongs to the current `beam_set`. If RayStation appends one evaluation per plan, every plan after the first exports the first plan's doses.
  - *Check:* CT_3 NPZs from two different original beams actually differ and match RayStation's display.
  - *Fix:* select the evaluation whose beam set matches the current plan.

- [ ] **1.2 CT_3 dose may be misaligned with CBCT3 (verify)** — `pipeline/step15_process_doses.m:659-669`
  If origin differs but dims match, only a warning; dose saved un-resampled and treated as on the reference grid. RS script uses `AllowGridExpansion=True` (`calc_beam_plan_doses.py:325`), which can change the grid.
  - *Check:* grep compress log for "Origin mismatch"; compare `corner_cm` of a CT_1 vs CT_3 NPZ for the same plan.
  - *Fix:* resample whenever origin OR spacing OR dims differ.

- [ ] **1.3 Dose resampler swaps x and y** — `pipeline/step15_process_doses.m:1604-1617`
  `ndgrid` gives x-origin/spacing to rows (Y) and y to columns. Triggers whenever dims differ.

- [ ] **1.4 CBCT labels vs dose labels come from different sources (verify)** — `step15_process_doses.m:1896-1913`
  CBCTs labelled by acquisition time; doses by RayStation exam name. Nothing confirms RS "CT 1" is the earlier CBCT.
  - *Check:* SeriesInstanceUID of RS exam "CT 1" == `CBCT1_resampled.series_uid`.
  - *Fix:* match by UID/exam name, not datetime.

- [ ] **1.5 Possible half-voxel offset (verify)** — `pipeline/step14_npz_to_mat.m:132`
  RayStation `DoseGrid.Corner` treated as voxel centre (DICOM IPP). If it is the outer corner, every dose is offset 1.25 mm/axis from CBCT, masks and Γρ map.
  - *Check:* compare against an exported RTDOSE's ImagePositionPatient.

- [ ] **1.6 Saved recons are max-normalized, not Gy** — `pipeline_simulate.m:85-86`
  Pins for `correction_factor`/`normalize` commented out → `normalize=true` from `utils/get_default_config.m:154` → each recon divided by its own max (`run_single_field_simulation.m:1340`) then summed into `total_recon` (`pipeline_simulate.m:533`). Contradicts the comment above the pins and `CLAUDE.md` ("Gy").

- [ ] **1.7 Totals mix both CTs and both plans** — `step15_process_doses.m:737, 803`
  `total_rs_dose` / `total_recon` sum reference+adapted × CT_1+CT_3. Used as the Step 2.5 noise-floor truth (`step25_segment_metrics.m:703`), `load_recon_dose_data` 'total' mode, Step 3, and sensor placement.

- [ ] **1.8 Step 2.5 silently reuses stale results** — `pipeline/step25_segment_metrics.m:319-324`
  Skips any segment with `segment_metrics`, without checking gamma criteria / cutoff / normalize / plan type. Summary then labels old numbers with new criteria. Same for noise-ensemble cache (`utils/noise_ensemble_error_bars.m:100-121`, keyed only by sim hash).
  - *Fix:* store and compare the criteria; recompute on mismatch.

- [ ] **1.9 Config hash misses recon-changing fields** — `utils/compute_sim_config_hash.m:84-113`
  Missing: `mask_recon_to_dose_region`, `Nt_scaling`, tissue tables (Γ/c/ρ), `elements_per_side`, `element_pitch_mm`, `element_size_mm`, `sensor_standoff_mm`, bone/side options, `uniform_*`.
  Lists `force_uniform_speed` but engine reads `force_uniform_sound_speed` (`run_single_field_simulation.m:221`).
  Input data not hashed — re-running Step 1.5 does not invalidate recons or `sensor_mask_<hash>.mat`.

- [ ] **1.10 Segment gets whole-beam MU, from one plan file** — `step15_process_doses.m:704, 1263-1269`
  One RTPLAN used for reference and adapted; each segment gets original beam's full MU, while real segment MU = Δweight × beam MU (`step06_explode_segments.m:268`). Pulse count inflated ~100×. Cancels in current legacy noise mode; breaks fixed `noise_amp_Pa` mode; saved `noise_stats.snr` non-physical.
  - *Fix:* compute segment MU from CumulativeMetersetWeight of the correct plan type.

- [ ] **1.11 ETHOS truth resampled by array size only** — `utils/load_recon_dose_data.m:603-646`, `pipeline/step3_analysis.m:411-437`
  `imresize3` to target dims ignoring origin/spacing; takes first RTDOSE*.dcm in folder. ETHOS-vs-RS comparisons invalid (`IncludeEthos` defaults true).

- [ ] **1.12 Noise floor uses a different engine and medium** — `utils/noise_ensemble_error_bars.m`
  Own copy of forward/TR; forces coupling bath (line 294, 822-845) while `pipeline_simulate.m:296` doesn't (contrary to "MANDATORY" note in `CLAUDE-SIMULATION_CONTEXT.md`); smooths sensor data across unordered channels (line 457); always Gaussian kernel.

### Latent (harmless now, bite on config change)

- [ ] **1.13 x/y spacing transposition** — `run_single_field_simulation.m:264-267, 798`; CalcGamma `width`; `step15_process_doses.m:1213, 1748` (PixelSpacing(1) is row spacing). Harmless while voxels are isotropic.
- [ ] **1.14 Iterative TR residual skips pulse/filter/deconv chain** — `run_single_field_simulation.m:1173-1179`. Only 1 iteration today (comments say 5/30).
- [ ] **1.15 `metrics_plan_type='any'` pairs by segment number only** — `step25_segment_metrics.m:447-448`; can pair reference-CT_1 with adapted-CT_3.
- [ ] **1.16 Air density clamped to 30 kg/m³** — `utils/create_acoustic_medium.m:194` (comment says 1; table says 1.2).
- [ ] **1.17 Run bookkeeping** — manifest keyed by sorted index (`pipeline_simulate.m:1383`), can defer forever after files are added; lock orphaned if worker dies between mkdir and status write (`pipeline_simulate.m:1244-1265`).
- [ ] **1.18 Sensor only checked against CT_1 body** — `utils/determine_sensor_mask_lateral.m:525-533`; CT_3 skin moving >5 mm outward puts sensor points in tissue.

### Docs vs code

- [ ] Active sensor is `determine_sensor_mask_lateral` (right flank, 32×4.35 mm = 139 mm) per `get_default_config.m:57-68`; comments/CLAUDE docs say anterior 10×10 cm.
- [ ] Anterior `determine_sensor_mask.m` still has the rotation bug (290° workaround).
- [ ] TR iterations: 1 actual vs "5"/"30" in comments.
- [ ] Unused local `define_tissue_tables` in `pipeline_simulate.m:720` — remove to avoid drift.
- [ ] Methods section must describe the config that actually ran.

---

## 2. Scientific framework issues

- [ ] **2.1 Truth information leaks into recons (highest priority)**
  - [ ] Recon masked to >1% of the field's *own true dose* and own CBCT body — `run_single_field_simulation.m:407-408, 1317-1325` (for blind CT_3 = CT_3 truth footprint; also deletes blind artefacts).
  - [ ] recon3 scaled by LS gain fitted to rs_CT3 — `step25_segment_metrics.m:480-485` (absorbs global amplitude change and mismatch amplitude error).
  - [ ] Sensor placed using total dose including CT_3 truths (minor).
  - *Fix:* use a-priori info only — mask both with planned (CT_1) footprint or not at all; scale recon3 with CT_1 gain or a single calibration. Test: does correlation of recon differential with true change survive?

- [ ] **2.2 Lopsided comparison ("inverse crime")**
  recon1 simulated and reconstructed on identical model/grid/medium; recon3 always carries medium mismatch → recon3 scores worse even with no dose change.
  - *Fix:* give CT_1 realistic model error (perturbed sound speed a few %, finer forward grid, sensor position error), reconstruct on nominal CT_1.

- [ ] **2.3 Null doesn't test the claim**
  Step 2.5 null = noise-only recon vs summed 4-plan truth, LS-scaled and masked to truth footprint (`noise_ensemble_error_bars.m:147-151, 810-814`). The 17–20% floor mostly reflects the mask.
  - *Fix:* per-segment null = spread of truth1-vs-recon1 pass rate over independent noise draws (as in `study_pass_rates_individual.m`), plus anatomy-only null from §3. Actually test against it.

- [ ] **2.4 Statistics understate uncertainty**
  SE = std/√n over segments (`studies/study_pass_rates_allsegments.m:829-832`); consecutive sliding-window segments are highly correlated.
  - *Fix:* beam as unit, block bootstrap, or mixed-effects; paired test with clustering. One patient/session = feasibility only.

- [ ] **2.5 Noise model not physical**
  - [ ] Legacy mode fixes SNR = 8 per field relative to own peak (`get_default_config.m:119, 137`; `run_single_field_simulation.m:1493-1508`).
  - [ ] No averaging over pulses (√N).
  - [ ] Wiener deconvolution amplifies noise up to 1/(2√λ) = 50× near ~100 kHz; "SNR" defined per sample at FS ≈ 4 MHz.
  - [ ] Noise per voxel, not per element; no element averaging in TR path.
  - *Fix:* tie to transducer NEP (Pa/√Hz), element area, pulse count; report SNR in recon domain.

- [ ] **2.6 Grüneisen values** — `utils/define_tissue_tables.m:26, 35`
  Water 0.11, fat 0.7, soft 1.0, bone 1.0; literature ≈ soft ~0.2 (37 °C), fat ~0.7–0.9, water 0.11 (20 °C)/~0.2 (37 °C). Blind fat↔soft misconversion 1.6× here vs ~3.5× (other direction) with literature values → "anatomy error" magnitude depends on this table. threshold_1 inconsistent (muscle 0.2, soft 1.0).
  - *Fix:* cite sources or sensitivity analysis.

- [ ] **2.7 CBCT-derived medium**
  - [ ] Hard HU thresholds after point resampling (`step15_process_doses.m:1924-1947`) → CBCT noise near −50 HU flips fat/soft between CT_1 and CT_3. Measure flip rate in unchanged region (paraspinal muscle); area-average HU first or use continuous HU→c/ρ map.
  - [ ] No coupling bath in pipeline → couch voxels get HU-derived properties (possibly "bone" 3200 m/s).
  - [ ] Confirm RayStation uses a commissioned CBCT HU→density table for both truths.
  - [ ] Check CBCT FOV truncation (SI extent) differs between CT_1 and CT_3.

- [ ] **2.8 Resolution vs 3%/3mm gamma**
  σ = 4 µs pulse (FWHM 9.4 µs) + λ=1e-4 → almost no signal >~130 kHz; 0.35 MHz filter (`run_single_field_simulation.m:974`) cuts DC ~8×, 100 kHz ~4×. Usable λ ≈ 12–30 mm, cm-scale PSF. Grid supports ≤290 kHz in tissue, ~69 kHz in air.
  - *Fix:* measure PSF; compare against PSF-blurred truth or justify criterion; justify pulse width vs measured linac pulse; grid convergence check.

- [ ] **2.9 Gamma setup choices**
  - [ ] `restrict=1` axis-only search (`step25_segment_metrics.m:592-593`) overestimates gamma for oblique gradients.
  - [ ] Eval mask from reference only (`step25_segment_metrics.m:350-354`) → extra dose where truth1 <10% never evaluated. Use union mask / both directions.
  - [ ] Add beam-level and composite comparisons (sum segment recons) and target/OAR dose metrics.

- [ ] **2.10 Monte Carlo noise in truths**
  truth1-vs-truth3 includes independent MC noise on thin segments (`calc_beam_plan_doses.py:15`). Quantify with repeat calc on same CT, or tighten uncertainty.

- [ ] **2.11 Minor physics**
  - [ ] Manual Dirichlet TR re-applies attenuation instead of compensating (k-Wave built-in TR compensates) — small <150 kHz except through bone.
  - [ ] Single `alpha_power` for all tissues.
  - [ ] MLC gap fix (step05) → truth is a modified plan; not comparable to clinical ETHOS dose.
  - [ ] State whether probe is couch-mounted or skin-mounted (currently fixed in grid coords).

---

## 3. Known "anatomy error" — make it measurable

- [ ] **Oracle run:** re-run CT_3 fields with `blind_recon_ct3 = false` (forward + recon on CT_3; own hash automatically). Use the same noise seed as the blind run.
- [ ] **Decompose** (TR is linear, shared noise draw):
  `recon3_blind − recon1 = (recon3_blind − recon3_oracle) + (recon3_oracle − recon1)`
  first term = acoustic model (fidelity) error; second = dose change as seen by the system. Respects the dose–CT coupling rule.
- [ ] **(Optional) false-positive class:** reconstruct CT_3 data on hybrid medium = CT_1 with CT_3 properties inside beam region (truth1 ∪ truth3 ≥ 10%). Residual vs oracle = sensing-path-only error.
- [ ] State in paper: if an in-session CBCT near delivery is available at reconstruction time, reconstructing only on CT_1 is conservative.
