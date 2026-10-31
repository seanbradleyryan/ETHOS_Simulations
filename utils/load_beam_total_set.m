function beam_set = load_beam_total_set(patient_id, session, config, beam, ct1_str, ct3_str)
%LOAD_BEAM_TOTAL_SET Load one beam's summed RS dose + beam recon on both CTs (beam mode).
%
%   beam_set = load_beam_total_set(patient_id, session, config, beam, ct1_str, ct3_str)
%
%   PURPOSE:
%   Beam-mode counterpart of load_recon_dose_data(Mode='set', Beam=beam). When
%   pipeline_simulate runs with reconstruction_mode = 'beam', Step 2 simulates
%   each beam's summed dose instead of each segment. This loads, for one beam,
%   the RS truth (total_dose_<CT>_B<n>.mat from pipeline_compress) and the beam
%   recon (total_dose_<CT>_B<n>_recon_<hash>.mat) on each of the two CTs. The
%   output has the same shape as load_recon_dose_data's 'set' output, so Step 2.5
%   and Step 3 use it unchanged. Each CT gets one field with seg_num 0 (= the
%   whole beam).
%
%   INPUTS:
%       patient_id - Char/string patient identifier (e.g. '1194203').
%       session    - Char/string session name (e.g. 'Session_1').
%       config     - Struct with:
%           .working_dir      - Base directory path.
%           .config_hash      - 8-char sim hash. When absent/empty the hash is
%                               computed from config (needs the full sim CONFIG).
%       beam       - Beam number.
%       ct1_str    - Reference CT label (e.g. 'CT_1').
%       ct3_str    - Second CT label (e.g. 'CT_3').
%
%   OUTPUTS:
%       beam_set - Struct with:
%           .config_hash
%           .metadata - processed/metadata.mat (spacing, dimensions, ...)
%           .fields   - 1x2 struct array (ct1_str, ct3_str), each with
%                       .recon_dose, .rs_dose, .rtplan (.seg_num = 0, .ct_label),
%                       .source_mat_filename, .recon_file
%
%   ALGORITHM:
%       1. Resolve the hash (config.config_hash, else compute_sim_config_hash).
%       2. For each CT: load the beam recon and the beam's summed RS dose
%          (densified when saved sparse).
%       3. Load metadata.mat for the grid spacing.
%   A missing file errors; callers skip that beam.
%
%   EXAMPLE:
%       config.working_dir      = get_repo_root();
%       config.config_hash      = 'a1b2c3d4';
%       B = load_beam_total_set('1194203', 'Session_1', config, 6, 'CT_1', 'CT_3');
%
%   DEPENDENCIES:
%       - compute_sim_config_hash, get_simulation_dir
%       - Beam totals from step15_process_doses; beam recons from pipeline_simulate.
%
%   See also: load_recon_dose_data, step25_segment_metrics, step3_analysis

    if ~ischar(patient_id) && ~isstring(patient_id)
        error('load_beam_total_set:InvalidInput', ...
            'patient_id must be a string or character array.');
    end
    patient_id = char(patient_id);
    if ~ischar(session) && ~isstring(session)
        error('load_beam_total_set:InvalidInput', ...
            'session must be a string or character array.');
    end
    session = char(session);

    if isfield(config, 'config_hash') && ~isempty(config.config_hash)
        hash8 = config.config_hash;
    else
        hash8 = compute_sim_config_hash(config);
    end
    processed_dir = fullfile(config.working_dir, 'RayStationFiles', ...
        patient_id, session, 'processed');
    sim_dir = get_simulation_dir(config.working_dir, patient_id, session, hash8);

    ct_labels = {ct1_str, ct3_str};
    for k = 1:2
        base       = sprintf('total_dose_%s_B%d', ct_labels{k}, beam);
        recon_file = fullfile(sim_dir, sprintf('%s_recon_%s.mat', base, hash8));
        R = load(recon_file, 'recon_dose');
        T = load(fullfile(processed_dir, [base, '.mat']));
        if isfield(T, 'beam_total_sparse')
            rs_dose = reshape(full(T.beam_total_sparse), T.beam_total_dims);
        else
            rs_dose = T.beam_total;
        end
        ct_fields(k) = struct('recon_dose', R.recon_dose, 'rs_dose', rs_dose, ...
            'rtplan', struct('seg_num', 0, 'ct_label', ct_labels{k}), ...
            'source_mat_filename', [base, '.mat'], 'recon_file', recon_file); %#ok<AGROW>
    end

    L = load(fullfile(processed_dir, 'metadata.mat'), 'metadata');
    beam_set = struct();
    beam_set.config_hash = hash8;
    beam_set.metadata    = L.metadata;
    beam_set.fields      = ct_fields;
end
