ScriptableDicomExport(..)
Exports specified DICOM datasets to either disk or SCP.
Parameters:
AnonymizationSettings: dict[unicode, Object] - Settings to be used if exported data is anonymized.

Parameters:
* Anonymize - Anonymize all exported datasets (bool).
* AnonymizedName - Patients name to set in anonymized datasets.
* AnonymizedID - Patient ID to set in anonymized datasets.
* RetainDates - Retain Longitudinal Temporal Information with Full Dates according to Dicom anonymization profile (bool).
* RetainDeviceIdentity - Retain Device Identity according to Dicom anonymization profile (bool).
* RetainInstitutionIdentity - Retain Institution Identity according to Dicom anonymization profile (bool).
* RetainUIDs - Retain UIDs according to Dicom anonymization profile (bool).
* RetainSafePrivateAttributes - Retain Safe Private Attributes (private attributes known by RayStation) (bool).
Connection: dict[unicode, Object] - Connection parameters when exporting to DICOM PACS.

A connection can be specified in one of two ways:
1. By setting parameters Node, Port, CallingAE and CalledAE
2. By setting only the parameter Title

If only Title is set, the other parameters will be read from the configured DICOM application entity (in Clinic Settings) with the given title.

Parameters:
* Node - IP address or hostname of remote SCP.
* Port - Port of remote SCP.
* CallingAE - AE title of local SCU. Leave empty to use hostname.
* CalledAE - AE title of remote SCP.
* Title - Title of configured DICOM application entity in Clinic Settings.
ExportFolderPath: unicode - Export target folder. Only used for file exports. Leave empty for SCP export
Examinations: list[unicode] - List of examination that shall be exported. Specified by Examination names
Argument snippet: Examinations = [examination.Name]
RtStructureSetsForExaminations: list[unicode] - List of examination names for which the structure set shall be exported.
Specified by examination names
Argument snippet: RtStructureSetsForExaminations = [examination.Name]
RtStructureSetsReferencedFromBeamSets: list[unicode] - List of beamset identifiers for which the referenced structure set that shall be exported
The identifier shall be specified as PlanName:DicomPlanLabel (ex. "Plan 1:BS 1" for plan "Plan 1" with beam set "BS 1")
Argument snippet: RtStructureSetsReferencedFromBeamSets = ["%s:%s"%(plan.Name, beam_set.DicomPlanLabel)]
Alternative snippet: RtStructureSetsReferencedFromBeamSets = [beam_set.BeamSetIdentifier()]
RtStructureSetsWithDicomUIDs: list[unicode] - List of DicomUIDs of SubStructureSets for which the RT Structure Set shall be exported.
Specified by ModificationInfo.DicomUID
Argument snippet: RtStructureSetsWithDicomUIDs = [structureSet.SubStructureSets[0].ModificationInfo.DicomUID]
BeamSets: list[unicode] - List of beamset identifiers that shall be exported
The identifier shall be specified as PlanName:DicomPlanLabel (ex. "Plan 1:BS 1" for plan "Plan 1" with beam set "BS 1")
Argument snippet: BeamSets = ["%s:%s"%(plan.Name, beam_set.DicomPlanLabel)]
Alternative snippet: BeamSets = [beam_set.BeamSetIdentifier()]
RtRadiationSetForBeamSets: list[unicode] - List of beamset identifiers for which the RT Radiation Set that shall be exported. Only CyberKnife beamsets are supported.
The identifier shall be specified as PlanName:DicomPlanLabel (ex. "Plan 1:BS 1" for plan "Plan 1" with beam set "BS 1")
Argument snippet: RtRadiationSetForBeamSets = ["%s:%s"%(plan.Name, beam_set.DicomPlanLabel)]
Alternative snippet: RtRadiationSetForBeamSets = [beam_set.BeamSetIdentifier()]
RtRadiationsForBeamSets: list[unicode] - List of beamset identifiers for which the RT Radiations that shall be exported. Only CyberKnife beamsets are supported.
The identifier shall be specified as PlanName:DicomPlanLabel (ex. "Plan 1:BS 1" for plan "Plan 1" with beam set "BS 1")
Argument snippet: RtRadiationsForBeamSets = ["%s:%s"%(plan.Name, beam_set.DicomPlanLabel)]
Alternative snippet: RtRadiationsForBeamSets = [beam_set.BeamSetIdentifier()]
PhysicalBeamSetDoseForBeamSets: list[unicode] - List of beamset identifiers for which the physical beam set dose shall be exported
The identifier shall be specified as PlanName:DicomPlanLabel (ex. "Plan 1:BS 1" for plan "Plan 1" with beam set "BS 1")
Argument snippet: PhysicalBeamSetDoseForBeamSets = ["%s:%s"%(plan.Name, beam_set.DicomPlanLabel)]
Alternative snippet: PhysicalBeamSetDoseForBeamSets = [beam_set.BeamSetIdentifier()]
EffectiveBeamSetDoseForBeamSets: list[unicode] - List of beamset identifiers for which the effective beam set dose shall be exported
The identifier shall be specified as PlanName:DicomPlanLabel (ex. "Plan 1:BS 1" for plan "Plan 1" with beam set "BS 1")
Argument snippet: EffectiveBeamSetDoseForBeamSets = ["%s:%s"%(plan.Name, beam_set.DicomPlanLabel)]
Alternative snippet: EffectiveBeamSetDoseForBeamSets = [beam_set.BeamSetIdentifier()]
PhysicalBeamDosesForBeamSets: list[unicode] - List of beamset identifiers for which all physical beam doses shall be exported
The identifier shall be specified as PlanName:DicomPlanLabel (ex. "Plan 1:BS 1" for plan "Plan 1" with beam set "BS 1")
Argument snippet: PhysicalBeamDosesForBeamSets = ["%s:%s"%(plan.Name, beam_set.DicomPlanLabel)]
Alternative snippet: PhysicalBeamDosesForBeamSets = [beam_set.BeamSetIdentifier()]
EffectiveBeamDosesForBeamSets: list[unicode] - List of beamset identifiers for which all effective beam doses shall be exported
The identifier shall be specified as PlanName:DicomPlanLabel (ex. "Plan 1:BS 1" for plan "Plan 1" with beam set "BS 1")
Argument snippet: EffectiveBeamDosesForBeamSets = ["%s:%s"%(plan.Name, beam_set.DicomPlanLabel)]
Alternative snippet: EffectiveBeamDosesForBeamSets = [beam_set.BeamSetIdentifier()]
SpatialRegistrationForExaminations: list[unicode] - List of examination pairs for which the registration object shall be exported.
The pair is specified as fromExaminationName:toExaminationName
Argument snippet: SpatialRegistrationForExaminations = ["%s:%s"%(fromExamination.Name, toExamination.Name)]
DeformableSpatialRegistrationsForExaminations: list[unicode] - List of examination pairs for which the registration object shall be exported.
The pair is specified as structureRegistrationGroupName:fromExaminationName:toExaminationName
Argument snippet: DeformableSpatialRegistrationsForExaminations = ["%s:%s:%s"%(case.PatientModel.StructureRegistrationGroups[0].Name, fromExamination.Name, toExamination.Name)]
TreatmentBeamDrrImages: list[unicode] - List of beamset identifiers for which all treatment beam DRRs shall be exported
The identifier shall be specified as PlanName:DicomPlanLabel (ex. "Plan 1:BS 1" for plan "Plan 1" with beam set "BS 1")
Argument snippet: TreatmentBeamDrrImages = ["%s:%s"%(plan.Name, beam_set.DicomPlanLabel)]
Alternative snippet: TreatmentBeamDrrImages = [beam_set.BeamSetIdentifier()]
If you want to specify a single beam or specific DrrSetting other then Default, the identifier shall be specified as PlanName:DicomPlanLabel:BeamName:DrrSettingName
(ex. "Plan 1:BS 1:B 1:DRR 1" for plan "Plan 1", with beam set "BS 1", and beam "B 1" and DrrSetting "DRR 1"
The last two argument can be omitted if wanted.
Not specifying beam will take all beams in the beam set
Not specifying DrrSetting will use the setting named "Default"
Argument snippet: TreatmentBeamDrrImages = ["%s:%s:%s:%s"%(plan.Name, beam_set.DicomPlanLabel, "", "")] # all beams with Default DrrSetting
Argument snippet: TreatmentBeamDrrImages = ["%s:%s:%s:%s"%(plan.Name, beam_set.DicomPlanLabel, "", "DRR 1")] # all beams with DrrSetting named "DRR 1"
Argument snippet: TreatmentBeamDrrImages = ["%s:%s:%s:%s"%(plan.Name, beam_set.DicomPlanLabel, beam.Name, "")] # only the selected beam with Default DrrSetting
Argument snippet: TreatmentBeamDrrImages = ["%s:%s:%s:%s"%(plan.Name, beam_set.DicomPlanLabel, beam.Name, "DRR 1")] # only the selected beam with DrrSetting named "DRR 1"
SetupBeamDrrImages: list[unicode] - List of beamset identifiers for which all setup beam DRRs shall be exported
The identifier shall be specified as PlanName:DicomPlanLabel (ex. "Plan 1:BS 1" for plan "Plan 1" with beam set "BS 1")
Argument snippet: SetupBeamDrrImages = ["%s:%s"%(plan.Name, beam_set.DicomPlanLabel)]
Alternative snippet: SetupBeamDrrImages = [beam_set.BeamSetIdentifier()]
If you want to specify a single beam or specific DrrSetting other then Default, the identifier shall be specified as PlanName:DicomPlanLabel:BeamName:DrrSettingName
(ex. "Plan 1:BS 1:B 1:DRR 1" for plan "Plan 1", with beam set "BS 1", and beam "B 1" and DrrSetting "DRR 1"
The last two argument can be omitted if wanted.
Not specifying beam will take all beams in the beam set
Not specifying DrrSetting will use the setting named "Default"
Argument snippet: SetupBeamDrrImages = ["%s:%s:%s:%s"%(plan.Name, beam_set.DicomPlanLabel, "", "")] # all beams with Default DrrSetting
Argument snippet: SetupBeamDrrImages = ["%s:%s:%s:%s"%(plan.Name, beam_set.DicomPlanLabel, "", "DRR 1")] # all beams with DrrSetting named "DRR 1"
Argument snippet: SetupBeamDrrImages = ["%s:%s:%s:%s"%(plan.Name, beam_set.DicomPlanLabel, beam.Name, "")] # only the selected beam with Default DrrSetting
Argument snippet: SetupBeamDrrImages = ["%s:%s:%s:%s"%(plan.Name, beam_set.DicomPlanLabel, beam.Name, "DRR 1")] # only the selected beam with DrrSetting named "DRR 1"
DicomFilter: unicode - Dicom filter to use during export. Specified by filter name.
IgnorePreConditionWarnings: bool - Switch for disabeling warnings. For any clinical export, warnings must be handled by first exporting with this argument set to False.
Use a try - except pattern to catch all warnings. After the warnings been handled, the export can bu run again with this attribute
set to True.
Code snippet:
try:
case.ScriptableDicomExport(...)
except SystemError as error:
HandleWarnings(error)
case.ScriptableDicomExport(... IgnorePreConditionWarnings=True)
RayGatewayTitle: unicode - RayGateway title for RayGateway export.
TransferSyntaxOverride: unicode - Transfer syntax override, for file export only. Overrides default setting in Clinic Settings.
Optional, shall be set to either "Implicit" or "Explicit" if set.
ExportAsBdspDose: bool - Export the attribute Beam Dose (300A,0084) as Beam Dose Specification Point (BDSP) dose.
This will override the 'Export Beam Dose (300A, 0084) as nominal contribution' setting on the machine in RayPhysics.