function sim_dir = get_simulation_dir(working_dir, patient_id, session)
%GET_SIMULATION_DIR Folder holding the simulation data of one patient/session.
%
%   sim_dir = get_simulation_dir(working_dir, patient_id, session)
%
%   PURPOSE:
%   Single definition of where pipeline_simulate's data lives:
%       <working_dir>/RayStationFiles/<patient_id>/<session>/processed/simulated_doses/
%   processed/ itself holds only the pipeline_compress (Step 1.5) outputs.
%   Every simulation output, lock and cache goes here, keyed by the 8-char
%   config hash (compute_sim_config_hash; gruneisen_method is in the hash, so
%   all methods share this one folder):
%   <base>_recon_<hash>.mat, total_recon_dose_<hash>.mat, sensor_mask_<hash>.mat,
%   cache_manifest_<hash>.mat, segment_metrics_summary_<hash>.mat,
%   full_field_dose_CT_*_<hash>.mat, noise_ensemble_*_<hash>.mat,
%   adjoint_power_cache_*.mat, the *.lock / *.metrics_lock folders, plus
%   config_registry.json (hash -> config), simulation_log.txt and field_status/.
%   Analysis outputs (figures, tables) go to AnalysisResults/ instead; see
%   config_label.
%
%   INPUTS:
%       working_dir - Repo root holding RayStationFiles/.
%       patient_id  - Char, e.g. '1194203'.
%       session     - Char, e.g. 'Session_1'.
%
%   OUTPUTS:
%       sim_dir - Char path. Not created here.
%
%   EXAMPLE:
%       sim_dir = get_simulation_dir(get_repo_root(), '1194203', 'Session_1');
%
%   See also: config_label, compute_sim_config_hash, pipeline_simulate

    sim_dir = fullfile(working_dir, 'RayStationFiles', char(patient_id), char(session), ...
        'processed', 'simulated_doses');
end
