Using scripting environment 'CPython 3.11 (64-bit)'
Running 'F:\ETHOS_Simulations\time_beam_dose_export.py'.
Script started: 20 May 2026, 11:24:56 (hr:min:sec)
Current patient is patient with UUID: fcaf0ce4-8f13-40e4-8827-c7e49d5e7fd6
No machine currently loaded

Connecting to RayStation. (Session id = 16484_3)
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
[M3b2] Readout via np.fromiter(chain.from_iterable(flat), count=N) ...
[M3c] Readout via np.array(tuple(flat)) ...
[M3d] Readout via bulk accessors ...
  Fastest readout: M3a list (0.3338 s)
  Building reference dose array for save tests ...
    built via M3a list
[M4] Save destinations (.npy durable + .npz compressed) ...
  Fastest save: M4d_server_C (.npy) at 0.0161 s -> C:/Temp/RayStationTiming_tmp
[M5] patient.Save() overhead ...
patient.Save()
patient.Save()
patient.Save()
patient.Save()
patient.Save()
[M6] ScriptableDicomExport to fastest save dir ...
patient.Cases['Case 1'].TreatmentPlans['IMRT Pancreas RT Intent Revision 16'].BeamSets['Total Fraction'].BeamSetIdentifier()
patient.Cases['Case 1'].ScriptableDicomExport(AnonymizationSettings=None, Connection=None, ExportFolderPath="C:/Temp/RayStationTiming_tmp", Examinations=[], RtStructureSetsForExaminations=[], RtStructureSetsReferencedFromBeamSets=[], RtStructureSetsWithDicomUIDs=[], BeamSets=["IMRT Pancreas RT Intent Revision 16:Total Fraction"], RtRadiationSetForBeamSets=[], RtRadiationsForBeamSets=[], PhysicalBeamSetDoseForBeamSets=["IMRT Pancreas RT Intent Revision 16:Total Fraction"], EffectiveBeamSetDoseForBeamSets=[], PhysicalBeamDosesForBeamSets=[], EffectiveBeamDosesForBeamSets=[], SpatialRegistrationForExaminations=[], DeformableSpatialRegistrationsForExaminations=[], TreatmentBeamDrrImages=[], SetupBeamDrrImages=[], DicomFilter=None, IgnorePreConditionWarnings=True, RayGatewayTitle=None, TransferSyntaxOverride=None, ExportAsBdspDose=False)
  ...Complete
patient.Cases['Case 1'].ScriptableDicomExport(AnonymizationSettings=None, Connection=None, ExportFolderPath="C:/Temp/RayStationTiming_tmp", Examinations=[], RtStructureSetsForExaminations=[], RtStructureSetsReferencedFromBeamSets=[], RtStructureSetsWithDicomUIDs=[], BeamSets=["IMRT Pancreas RT Intent Revision 16:Total Fraction"], RtRadiationSetForBeamSets=[], RtRadiationsForBeamSets=[], PhysicalBeamSetDoseForBeamSets=["IMRT Pancreas RT Intent Revision 16:Total Fraction"], EffectiveBeamSetDoseForBeamSets=[], PhysicalBeamDosesForBeamSets=[], EffectiveBeamDosesForBeamSets=[], SpatialRegistrationForExaminations=[], DeformableSpatialRegistrationsForExaminations=[], TreatmentBeamDrrImages=[], SetupBeamDrrImages=[], DicomFilter=None, IgnorePreConditionWarnings=True, RayGatewayTitle=None, TransferSyntaxOverride=None, ExportAsBdspDose=False)
  ...Complete
patient.Cases['Case 1'].ScriptableDicomExport(AnonymizationSettings=None, Connection=None, ExportFolderPath="C:/Temp/RayStationTiming_tmp", Examinations=[], RtStructureSetsForExaminations=[], RtStructureSetsReferencedFromBeamSets=[], RtStructureSetsWithDicomUIDs=[], BeamSets=["IMRT Pancreas RT Intent Revision 16:Total Fraction"], RtRadiationSetForBeamSets=[], RtRadiationsForBeamSets=[], PhysicalBeamSetDoseForBeamSets=["IMRT Pancreas RT Intent Revision 16:Total Fraction"], EffectiveBeamSetDoseForBeamSets=[], PhysicalBeamDosesForBeamSets=[], EffectiveBeamDosesForBeamSets=[], SpatialRegistrationForExaminations=[], DeformableSpatialRegistrationsForExaminations=[], TreatmentBeamDrrImages=[], SetupBeamDrrImages=[], DicomFilter=None, IgnorePreConditionWarnings=True, RayGatewayTitle=None, TransferSyntaxOverride=None, ExportAsBdspDose=False)
  ...Complete
patient.Cases['Case 1'].ScriptableDicomExport(AnonymizationSettings=None, Connection=None, ExportFolderPath="C:/Temp/RayStationTiming_tmp", Examinations=[], RtStructureSetsForExaminations=[], RtStructureSetsReferencedFromBeamSets=[], RtStructureSetsWithDicomUIDs=[], BeamSets=["IMRT Pancreas RT Intent Revision 16:Total Fraction"], RtRadiationSetForBeamSets=[], RtRadiationsForBeamSets=[], PhysicalBeamSetDoseForBeamSets=["IMRT Pancreas RT Intent Revision 16:Total Fraction"], EffectiveBeamSetDoseForBeamSets=[], PhysicalBeamDosesForBeamSets=[], EffectiveBeamDosesForBeamSets=[], SpatialRegistrationForExaminations=[], DeformableSpatialRegistrationsForExaminations=[], TreatmentBeamDrrImages=[], SetupBeamDrrImages=[], DicomFilter=None, IgnorePreConditionWarnings=True, RayGatewayTitle=None, TransferSyntaxOverride=None, ExportAsBdspDose=False)
  ...Complete
patient.Cases['Case 1'].ScriptableDicomExport(AnonymizationSettings=None, Connection=None, ExportFolderPath="C:/Temp/RayStationTiming_tmp", Examinations=[], RtStructureSetsForExaminations=[], RtStructureSetsReferencedFromBeamSets=[], RtStructureSetsWithDicomUIDs=[], BeamSets=["IMRT Pancreas RT Intent Revision 16:Total Fraction"], RtRadiationSetForBeamSets=[], RtRadiationsForBeamSets=[], PhysicalBeamSetDoseForBeamSets=["IMRT Pancreas RT Intent Revision 16:Total Fraction"], EffectiveBeamSetDoseForBeamSets=[], PhysicalBeamDosesForBeamSets=[], EffectiveBeamDosesForBeamSets=[], SpatialRegistrationForExaminations=[], DeformableSpatialRegistrationsForExaminations=[], TreatmentBeamDrrImages=[], SetupBeamDrrImages=[], DicomFilter=None, IgnorePreConditionWarnings=True, RayGatewayTitle=None, TransferSyntaxOverride=None, ExportAsBdspDose=False)
  ...Complete

======================================================================
SUMMARY
======================================================================
  Plan          : IMRT Pancreas RT Intent Revision 16
patient.Cases['Case 1'].TreatmentPlans['IMRT Pancreas RT Intent Revision 16'].BeamSets['Total Fraction'].BeamSetIdentifier()
  BeamSet id    : IMRT Pancreas RT Intent Revision 16:Total Fraction
  Test beam     : 0
  Grid          : 141 x 96 x 141  (1908576 voxels, 7.28 MB float32)
  Voxel (cm)    : x=0.2500 y=0.2500 z=0.2500
  N_REPEATS     : 5
  RS version    : unknown
  BeamDose src  : FractionDose.BeamDoses[0]

  Measurement                                      |       Mean +/- Std (s) | Notes
  ------------------------------------------------------------------------------------------------
  M1  Grid metadata fetch                          |     0.0049 +/-  0.0002 | 
  M2  DoseData handle fetch (no iteration)         |     0.3335 +/-  0.0151 | 
  M3a Readout: np.array(list(flat))                |     0.3338 +/-  0.0093 | ok
  M3b Readout: np.fromiter(flat, count)            |                    n/a | failed: setting an array element with a sequence.
  M3b2 Readout: np.fromiter(chain, count)          |                    n/a | failed: setting an array element with a sequence.
  M3c Readout: np.array(tuple(flat))               |     0.3953 +/-  0.1364 | ok
  M3d Readout: DoseValues.DoseDataAsArray          |                    n/a | not available in this build (Object has no member 'DoseDataAsArray'.)
  M3d Readout: DoseValues.GetDoseValues()          |                    n/a | not available in this build (Object has no member 'GetDoseValues'.)
  M3d Readout: GetDoseValues()                     |                    n/a | not available in this build (Object has no member 'GetDoseValues'.)
  M4  M4a_cwd .npy                                 |     0.0526 +/-  0.0181 | 7.28 MB on disk  (F:\ETHOS_Simulations)
  M4  M4a_cwd .npz                                 |     0.1990 +/-  0.0046 | 4.13 MB on disk  (F:\ETHOS_Simulations)
  M4  M4b_server_F .npy                            |     0.0362 +/-  0.0006 | 7.28 MB on disk  (F:/RayStationTiming_tmp)
  M4  M4b_server_F .npz                            |     0.1988 +/-  0.0034 | 4.13 MB on disk  (F:/RayStationTiming_tmp)
  M4  M4c_citrix_C .npy                            |    23.4181 +/-  0.7182 | 7.28 MB on disk  (//Client/C$/Users/80030361/Documents/RayStationTiming_tmp)
  M4  M4c_citrix_C .npz                            |    20.2263 +/-  0.7742 | 4.13 MB on disk  (//Client/C$/Users/80030361/Documents/RayStationTiming_tmp)
  M4  M4d_server_C .npy                            |     0.0161 +/-  0.0018 | 7.28 MB on disk  (C:/Temp/RayStationTiming_tmp)
  M4  M4d_server_C .npz                            |     0.1724 +/-  0.0016 | 4.13 MB on disk  (C:/Temp/RayStationTiming_tmp)
  M5  patient.Save()                               |     0.0793 +/-  0.0236 | 
  M6  ScriptableDicomExport (1 beam)               |     3.5515 +/-  0.1504 | ok

======================================================================
FOOTER
======================================================================
  Best readout      : M3a list  (0.3338 s vs M3a 0.3338 s, 1.00x)
  Best save dest    : M4d_server_C (.npy)  0.0161 s vs Citrix .npz 20.2263 s, 1256.23x
  6000-dose run     : current 34.27 h  ->  optimized 0.58 h  (58.76x)
patient.Cases['Case 1'].TreatmentPlans['IMRT Pancreas RT Intent Revision 16'].BeamSets['Total Fraction'].BeamSetIdentifier()

  CSV written: F:\ETHOS_Simulations\timing_results_20260520_112907.csv
  Cleanup: removed 2/42 scratch file(s).

Script completed: 20 May 2026, 11:29:07 (hr:min:sec)