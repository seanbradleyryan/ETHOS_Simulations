function results = step205_run_full_field_simulation(patient_id, session, config)
%STEP205_RUN_FULL_FIELD_SIMULATION Simulate each CT's total RS dose as one field.
%
%   results = step205_run_full_field_simulation(patient_id, session, config)
%   results = step205_run_full_field_simulation()   % standalone: 1194203, Session_1
%
%   PURPOSE:
%   Step 2 simulates every beam/segment field on its own. This step simulates
%   the whole plan in ONE k-Wave run per CT: the per-CT total RayStation dose
%   written by step15_process_doses (pipeline_compress.m) is passed to
%   run_single_field_simulation as if it were a single field, so it goes
%   through the same sensor placement, p0 = D*Gamma*rho, forward propagation,
%   reconstruction and pressure -> dose conversion. The CT_1 / CT_3 split
%   matches Step 2 (run_field_dispatch in pipeline_simulate):
%       CT_1 total -> full access: forward and reconstruction on CT_1
%       CT_3 total -> blind: forward on CT_3, reconstruction on CT_1
%                     (only when config.blind_recon_ct3 is true, as in Step 2)
%   The plan sensor mask Step 2 saved (sensor_mask_<hash>.mat) is reused. The
%   four Step 2.5 comparisons (gamma + local SSIM) are then run on the totals.
%
%   The analysis should compare these outputs with the SUM of the per-field
%   Step 2 recons (total_recon_dose_<hash>.mat).
%
%   CACHE: outputs are keyed on the sim config hash, like Step 2. A CT whose
%   full_field_dose_<CT>_<hash>.mat exists is not simulated again. When both
%   exist and the CT_1 file already holds full_field_metrics, the stored result
%   is returned. config.step205_overwrite = true redoes everything.
%   A <sim_dir>/full_field_<hash>.lock folder stops two pipeline instances from
%   running this step at once. A crashed run can leave it behind: delete it.
%
%   INPUTS:
%       patient_id - Char/string patient ID (e.g. '1194203').
%       session    - Char/string session name (e.g. 'Session_2').
%       config     - Simulation CONFIG, as built in pipeline_simulate
%                    (defaults in []):
%           .working_dir            REQUIRED. Root holding RayStationFiles/ and
%                                   SimulationResults/.
%           .gruneisen_method       REQUIRED. SimulationResults subfolder.
%           .config_hash            [computed] 8-char sim hash. Computed with
%                                   compute_sim_config_hash(config) when absent.
%           .sensor_placement_method  A determine_sensor_mask method needs the
%                                   Step 2 sensor_mask_<hash>.mat on disk.
%           .blind_recon_ct3        [false] Reconstruct the CT_3 total on CT_1.
%           .gamma_dose_pct         [3.0]   Gamma dose criterion (%).
%           .gamma_dist_mm          [3.0]   Gamma distance criterion (mm).
%           .gamma_dose_cutoff_pct  [10.0]  Eval mask: voxels >= this % of the
%                                   reference max.
%           .metrics_normalize      [true]  Scale each recon by its least-squares
%                                   gain to its own-CT truth before comparing.
%           .step205_overwrite      [false] Simulate and compare again even if
%                                   the outputs are on disk.
%           .Nt_scaling             Ignored: forced to 2 in this step.
%           Plus every field run_single_field_simulation and
%           create_acoustic_medium read.
%
%   OUTPUTS:
%       results - Struct:
%           .patient_id, .session, .config_hash
%           .status       'complete' | 'cached' | 'locked'
%           .recon_files  {CT_1 file, CT_3 file}
%           .metrics      full_field_metrics struct ([] when locked)
%
%       Written to SimulationResults/<pid>/<session>/<method>/ (next to the
%       Step 2 recons):
%           full_field_dose_CT_1_<hash>.mat  recon_dose, config_hash, ct_label,
%                                            blind_recon, meterset_MU,
%                                            noise_stats, full_field_metrics
%           full_field_dose_CT_3_<hash>.mat  same, without full_field_metrics
%       full_field_metrics has the layout of Step 2.5's segment_metrics:
%       .comparison(d) with .name, .mask_idx, .gamma_vals, .gamma_pass_rate,
%       .ssim_vals, .ssim_mean, plus .vol_size, .spacing and the LS gains.
%       Re-expand a map with: M = nan(m.vol_size); M(c.mask_idx) = c.gamma_vals;
%
%   ALGORITHM:
%       1. Return the stored result if both recons and the metrics are on disk.
%       2. Claim the lock folder.
%       3. Load the per-CT totals (total_dose_CT_1.mat / total_dose_CT_3.mat).
%       4. For each CT not on disk yet: load the CBCTs, build the media with
%          create_acoustic_medium, sum the metersets of that CT's fields (sets
%          the pulse count), run run_single_field_simulation with the
%          precomputed sensor, save the recon.
%       5. Four comparisons (truth1_vs_truth3, truth1_vs_recon1,
%          truth1_vs_recon3, truth3_vs_recon3): global gamma + local SSIM over
%          the reference's 10% mask. Saved into the CT_1 file.
%       6. Release the lock.
%
%   EXAMPLE:
%       addpath(genpath('utils')); addpath('pipeline');
%       CONFIG = get_default_config();
%       CONFIG.blind_recon_ct3             = true;    % match pipeline_simulate
%       CONFIG.blind_recon_reference_label = 'CT_1';  % so the hash matches
%       r = step205_run_full_field_simulation('1194203', 'Session_2', CONFIG);
%
%   DEPENDENCIES:
%       - run_single_field_simulation, create_acoustic_medium, k-Wave
%       - list_processed_field_doses, compute_sim_config_hash, least_squares_gain
%       - CalcGamma, compute_local_ssim (Image Processing Toolbox)
%       - Step 1.5 outputs in RayStationFiles/<pid>/<session>/processed/
%       - Step 2 sensor_mask_<hash>.mat (for determine_sensor_mask methods)
%
%   See also: run_single_field_simulation, pipeline_simulate,
%             step25_segment_metrics, step3_analysis

    %% ======================== STANDALONE DEFAULTS ========================
    % Called with no arguments (e.g. the editor's Run button): patient 1194203,
    % Session_1, and the CONFIG pipeline_simulate builds. blind_recon_ct3 and
    % blind_recon_reference_label are in the config hash, so they must match
    % pipeline_simulate for the Step 2 sensor_mask_<hash>.mat to be found.
    if nargin == 0
        repo_root = fileparts(fileparts(mfilename('fullpath')));
        addpath(genpath(fullfile(repo_root, 'utils')));
        patient_id = '1194203';
        session    = 'Session_1';
        config     = get_default_config();
        config.blind_recon_ct3             = true;
        config.blind_recon_reference_label = 'CT_1';
        config.plot_results                = false;   % no engine figures, as in pipeline_simulate
    end

    %% ======================== INPUT VALIDATION ========================
    if ~ischar(patient_id) && ~isstring(patient_id)
        error('step205_run_full_field_simulation:InvalidInput', ...
            'patient_id must be a string or character array.');
    end
    patient_id = char(patient_id);
    if ~ischar(session) && ~isstring(session)
        error('step205_run_full_field_simulation:InvalidInput', ...
            'session must be a string or character array.');
    end
    session = char(session);
    if ~isstruct(config) || ~isfield(config, 'working_dir') || ~isfield(config, 'gruneisen_method')
        error('step205_run_full_field_simulation:MissingConfig', ...
            'config must be a struct containing working_dir and gruneisen_method.');
    end

    %% ======================== CONFIG DEFAULTS ========================
    if ~isfield(config, 'config_hash') || isempty(config.config_hash)
        config.config_hash = compute_sim_config_hash(config);
    end
    config = default_field(config, 'blind_recon_ct3',       false);
    config = default_field(config, 'gamma_dose_pct',        3.0);
    config = default_field(config, 'gamma_dist_mm',         3.0);
    config = default_field(config, 'gamma_dose_cutoff_pct', 10.0);
    config = default_field(config, 'metrics_normalize',     true);
    config = default_field(config, 'step205_overwrite',     false);
    % Forced for this step, whatever the caller set. Nt_scaling is not in the
    % config hash, so recons already on disk under this hash are still reused.
    config.Nt_scaling = 2;

    hash8         = config.config_hash;
    ct_labels     = {'CT_1', 'CT_3'};
    sim_dir       = fullfile(config.working_dir, 'SimulationResults', patient_id, ...
        session, config.gruneisen_method);
    processed_dir = fullfile(config.working_dir, 'RayStationFiles', patient_id, ...
        session, 'processed');
    recon_files   = {fullfile(sim_dir, sprintf('full_field_dose_CT_1_%s.mat', hash8)), ...
                     fullfile(sim_dir, sprintf('full_field_dose_CT_3_%s.mat', hash8))};

    results             = struct();
    results.patient_id  = patient_id;
    results.session     = session;
    results.config_hash = hash8;
    results.recon_files = recon_files;
    results.metrics     = [];

    fprintf('\n=========================================================\n');
    fprintf('  [STEP 2.05] Full-Field Simulation (total RS dose per CT)\n');
    fprintf('  Patient %s | %s | hash %s | blind CT_3 = %d\n', ...
        patient_id, session, hash8, config.blind_recon_ct3);
    fprintf('=========================================================\n');

    %% ======================== CACHE CHECK ========================
    need_sim = ~[isfile(recon_files{1}), isfile(recon_files{2})];
    if config.step205_overwrite
        need_sim(:) = true;
    end
    if ~any(need_sim) && ~isempty(whos('-file', recon_files{1}, 'full_field_metrics'))
        L = load(recon_files{1}, 'full_field_metrics');
        results.metrics = L.full_field_metrics;
        results.status  = 'cached';
        fprintf('[STEP 2.05] Already computed for hash %s; loaded from:\n  %s\n', ...
            hash8, recon_files{1});
        print_metrics(results.metrics);
        return;
    end

    %% ======================== LOCK ========================
    % Atomic folder lock (same mkdir trick as Steps 2 / 2.5) so two pipeline
    % instances never simulate or write these files at the same time.
    if ~exist(sim_dir, 'dir'), mkdir(sim_dir); end
    lock_dir = fullfile(sim_dir, sprintf('full_field_%s.lock', hash8));
    [ok, ~, msgid] = mkdir(lock_dir);
    if ~ok || strcmpi(msgid, 'MATLAB:MKDIR:DirectoryExists')
        fprintf(['[STEP 2.05] Lock held by another instance; skipping.\n' ...
                 '           If no other run is active, delete: %s\n'], lock_dir);
        results.status = 'locked';
        return;
    end

    try
        %% ======================== PER-CT TOTAL RS DOSES ========================
        % Used both as the simulation input and as the truth for the metrics.
        [field_index, dose_metadata] = list_processed_field_doses(patient_id, session, config);
        rs = struct();
        for k = 1:2
            lbl   = ct_labels{k};
            fpath = fullfile(processed_dir, sprintf('total_dose_%s.mat', lbl));
            if ~isfile(fpath)
                error('step205_run_full_field_simulation:NoTotalDose', ...
                    'Per-CT total dose not found: %s\nRun pipeline_compress (Step 1.5) first.', fpath);
            end
            L = load(fpath);
            if isfield(L, 'ct_total_sparse')
                rs.(lbl) = double(reshape(full(L.ct_total_sparse), L.ct_total_dims));
            else
                rs.(lbl) = double(L.ct_total);
            end
            fprintf('  %s total RS dose: max %.4f Gy\n', lbl, max(rs.(lbl)(:)));
        end

        %% ======================== SIMULATION INPUTS ========================
        if any(need_sim)
            % Plan sensor mask saved by Step 2 (the one every field used).
            % Other placement methods build their mask inline, as in Step 2.
            sensor_file        = fullfile(sim_dir, sprintf('sensor_mask_%s.mat', hash8));
            precomputed_sensor = [];
            if isfile(sensor_file)
                L = load(sensor_file, 'precomputed_sensor');
                precomputed_sensor = L.precomputed_sensor;
                fprintf('  Using precomputed plan sensor: %s\n', sensor_file);
            elseif ismember(config.sensor_placement_method, ...
                    {'determine_sensor_mask', 'determine_sensor_mask_lateral'})
                error('step205_run_full_field_simulation:NoSensorMask', ...
                    ['Precomputed sensor mask not found: %s\n' ...
                     'Run Step 2 (pipeline_simulate) with this CONFIG first.'], sensor_file);
            end

            beam_metadata = [];
            if isfield(dose_metadata, 'beam_metadata')
                beam_metadata = dose_metadata.beam_metadata;
            end

            % Both CBCTs and media: CT_1 is also the CT_3 blind recon geometry.
            L = load(fullfile(processed_dir, 'CBCT1_resampled.mat'), 'CBCT1_resampled');
            cbct.CT_1 = L.CBCT1_resampled;
            L = load(fullfile(processed_dir, 'CBCT3_resampled.mat'), 'CBCT3_resampled');
            cbct.CT_3 = L.CBCT3_resampled;
            medium.CT_1 = create_acoustic_medium(cbct.CT_1, config);
            medium.CT_3 = create_acoustic_medium(cbct.CT_3, config);

            % Pulse count: a CT's total dose is delivered over every pulse of
            % that CT's fields, so use the summed meterset of those fields.
            % Loads one field file at a time.
            fprintf('  Summing metersets over %d field file(s)...\n', numel(field_index));
            meterset = struct('CT_1', 0, 'CT_3', 0);
            for i = 1:numel(field_index)
                L   = load(field_index(i).file, 'field_dose');
                lbl = strrep(L.field_dose.ct_label, '-', '_');
                if isfield(meterset, lbl)
                    meterset.(lbl) = meterset.(lbl) + L.field_dose.meterset;
                end
            end
        end

        %% ======================== SIMULATIONS ========================
        recon = struct();
        for k = 1:2
            lbl = ct_labels{k};
            if ~need_sim(k)
                L = load(recon_files{k}, 'recon_dose');
                recon.(lbl) = L.recon_dose;
                fprintf('  %s: cached recon loaded (%s)\n', lbl, recon_files{k});
                continue;
            end

            % CT_3 -> blind (forward on CT_3, reconstruct on CT_1), as in Step 2.
            blind_recon  = strcmp(lbl, 'CT_3') && config.blind_recon_ct3;
            medium_recon = [];
            if blind_recon
                medium_recon = medium.CT_1;
            end
            sim_config             = config;
            sim_config.blind_recon = blind_recon;

            fd              = struct();
            fd.dose_Gy      = rs.(lbl);
            fd.meterset     = meterset.(lbl);
            fd.gantry_angle = 0;   % whole plan: no single gantry angle

            fprintf('\n  %s total: %.1f MU, blind = %d. Simulating...\n', ...
                lbl, fd.meterset, blind_recon);
            [recon_dose, sim_results] = run_single_field_simulation(fd, cbct.(lbl), ...
                medium.(lbl), beam_metadata, sim_config, precomputed_sensor, medium_recon);
            recon_dose = gather(recon_dose);

            % Same noise diagnostics Step 2 saves with each field recon.
            noise_stats = struct('sensor_signal_peak', NaN, 'noise_amp', NaN, 'snr', NaN);
            names = fieldnames(noise_stats);
            for n = 1:numel(names)
                if isfield(sim_results, names{n}) && ~isempty(sim_results.(names{n}))
                    noise_stats.(names{n}) = double(gather(sim_results.(names{n})));
                end
            end

            config_hash = hash8;            %#ok<NASGU> saved variable
            ct_label    = lbl;              %#ok<NASGU> saved variable
            meterset_MU = meterset.(lbl);   %#ok<NASGU> saved variable
            save(recon_files{k}, 'recon_dose', 'config_hash', 'ct_label', ...
                'blind_recon', 'meterset_MU', 'noise_stats', '-v7.3');
            recon.(lbl) = recon_dose;
            fprintf('  Saved %s (max %.4f Gy)\n', recon_files{k}, max(recon_dose(:)));
        end
        clear medium cbct

        %% ======================== FOUR STEP 2.5 COMPARISONS ========================
        gain_CT1 = 1;
        gain_CT3 = 1;
        if config.metrics_normalize
            gain_CT1 = least_squares_gain(rs.CT_1, recon.CT_1);
            gain_CT3 = least_squares_gain(rs.CT_3, recon.CT_3);
        end
        vols = struct('rs_CT1', rs.CT_1, 'rs_CT3', rs.CT_3, ...
            'recon_CT1', double(recon.CT_1) * gain_CT1, ...
            'recon_CT3', double(recon.CT_3) * gain_CT3);

        % (name, reference, target). The reference builds the eval mask.
        comparisons = { ...
            'truth1_vs_truth3', 'rs_CT1', 'rs_CT3'; ...
            'truth1_vs_recon1', 'rs_CT1', 'recon_CT1'; ...
            'truth1_vs_recon3', 'rs_CT1', 'recon_CT3'; ...
            'truth3_vs_recon3', 'rs_CT3', 'recon_CT3'};
        spacing   = dose_metadata.spacing(:)';
        cutoff_fr = config.gamma_dose_cutoff_pct / 100;

        comp = struct('name', {}, 'ref_field', {}, 'tgt_field', {}, 'mask_idx', {}, ...
            'gamma_vals', {}, 'gamma_pass_rate', {}, 'ssim_vals', {}, 'ssim_mean', {});
        for d = 1:size(comparisons, 1)
            ref  = vols.(comparisons{d, 2});
            tgt  = vols.(comparisons{d, 3});
            mask = ref >= cutoff_fr * max(ref(:));
            fprintf('  [%d/4] %s: gamma + SSIM over %d voxels...\n', ...
                d, comparisons{d, 1}, nnz(mask));

            % Global gamma, same CalcGamma call as Step 2.5.
            ref_struct = struct('start', [0, 0, 0], 'width', spacing, 'data', ref);
            tgt_struct = struct('start', [0, 0, 0], 'width', spacing, 'data', tgt);
            gmap = CalcGamma(ref_struct, tgt_struct, config.gamma_dose_pct, ...
                config.gamma_dist_mm, 'local', 0, 'limit', 2 * config.gamma_dist_mm, ...
                'restrict', 1, 'cpu', 1);
            g_vals = single(gmap(mask));

            [smap, s_frac] = compute_local_ssim(ref, tgt, mask);
            s_vals = single([]);
            if ~isempty(smap)
                s_vals = single(smap(mask));
            end

            comp(d) = struct('name', comparisons{d, 1}, ...
                'ref_field', comparisons{d, 2}, 'tgt_field', comparisons{d, 3}, ...
                'mask_idx', uint32(find(mask)), 'gamma_vals', g_vals, ...
                'gamma_pass_rate', 100 * mean(double(g_vals) <= 1), ...
                'ssim_vals', s_vals, 'ssim_mean', 100 * s_frac);
        end

        full_field_metrics                       = struct();
        full_field_metrics.config_hash           = hash8;
        full_field_metrics.spacing               = spacing;
        full_field_metrics.vol_size              = size(rs.CT_1);
        full_field_metrics.normalize             = config.metrics_normalize;
        full_field_metrics.recon_ct1_gain        = gain_CT1;
        full_field_metrics.recon_ct3_gain        = gain_CT3;
        full_field_metrics.gamma_dose_pct        = config.gamma_dose_pct;
        full_field_metrics.gamma_dist_mm         = config.gamma_dist_mm;
        full_field_metrics.gamma_dose_cutoff_pct = config.gamma_dose_cutoff_pct;
        full_field_metrics.comparison            = comp;
        save(recon_files{1}, 'full_field_metrics', '-append');
    catch ME
        if isfolder(lock_dir), rmdir(lock_dir, 's'); end
        rethrow(ME);
    end
    if isfolder(lock_dir), rmdir(lock_dir, 's'); end

    results.metrics = full_field_metrics;
    results.status  = 'complete';
    fprintf('\n[STEP 2.05] Complete. Saved in %s\n', sim_dir);
    print_metrics(full_field_metrics);
end


function config = default_field(config, name, value)
    if ~isfield(config, name) || isempty(config.(name))
        config.(name) = value;
    end
end


function print_metrics(m)
    fprintf('  %-18s %9s %9s\n', 'Comparison', 'Gamma %', 'SSIM %');
    for d = 1:numel(m.comparison)
        fprintf('  %-18s %9.1f %9.1f\n', m.comparison(d).name, ...
            m.comparison(d).gamma_pass_rate, m.comparison(d).ssim_mean);
    end
end
