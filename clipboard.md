Using scripting environment 'CPython 3.11 (64-bit)'
Running 'F:\ETHOS_Simulations\time_beam_dose_export.py'.
Script started: 20 May 2026, 11:15:34 (hr:min:sec)
Current patient is patient with UUID: fcaf0ce4-8f13-40e4-8827-c7e49d5e7fd6
No machine currently loaded

Connecting to RayStation. (Session id = 16484_2)
======================================================================
RayStation per-beam dose readout/save timing
======================================================================
  Plan        : IMRT Pancreas RT Intent Revision 16
patient.Cases['Case 1'].TreatmentPlans['IMRT Pancreas RT Intent Revision 16'].BeamSets['Total Fraction'].BeamSetIdentifier()
  BeamSet id  : IMRT Pancreas RT Intent Revision 16:Total Fraction
  Test beam   : 0
  N_REPEATS   : 5
  RS version  : unknown
  CWD         : F:\ETHOS_Simulations
  BeamDose src: FractionDose.BeamDoses[0]
  Grid        : 141 x 96 x 141  (1908576 voxels, 7.28 MB float32)
  Voxel (cm)  : x=0.2500 y=0.2500 z=0.2500
  Corner (cm) : x=-21.5361 y=-9.3315 z=-16.0174

[M1] Grid metadata fetch ...
[M2] DoseData handle fetch (no iteration) ...
[M3a] Readout via np.array(list(flat)) ...
[M3b] Readout via np.fromiter(flat, count=N) ...

ERROR: setting an array element with a sequence.
TypeError: only length-1 arrays can be converted to Python scalars

The above exception was the direct cause of the following exception:

Traceback (most recent call last):
  File "F:\ETHOS_Simulations\time_beam_dose_export.py", line 383, in <module>
    m3b_times = measure_m3b_fromiter(bd, nz, ny, nx, N_REPEATS)
                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "F:\ETHOS_Simulations\time_beam_dose_export.py", line 175, in measure_m3b_fromiter
    arr  = np.fromiter(flat, dtype=np.float32, count=total).reshape(nz, ny, nx)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
ValueError: setting an array element with a sequence.


Script completed: 20 May 2026, 11:15:39 (hr:min:sec)