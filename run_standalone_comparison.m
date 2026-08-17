%% =========================================================================
%  RUN_STANDALONE_COMPARISON.m   (wrapper)
%  Two-dose k-Wave comparison: a listed dose (A) and its counterpart on the
%  other CT image (B). Each dose is reconstructed through the shared engine via
%  run_standalone_field -- the reference-CT dose full-access, the other-CT dose
%  BLIND (forward on its own CT, time-reversal on the reference CT). The two
%  reconstructions are then gamma-compared and saved for analysis.
%
%  Thin driver: all physics lives in run_single_field_simulation; setup + medium
%  build live in run_standalone_field; plots/gamma use the utils/ helpers.
%  =========================================================================

clear; clc; close all;

% Ensure the moved helper functions in utils/ are on the path (run from root).
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'utils')));
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'pipeline')));  % step* functions live here

%% ========================= CONFIGURATION ================================

CONFIG = get_default_config();
CONFIG.mask_recon_to_dose_region     = false;    % zero recon dose outside the >1% dose mask


% --- Comparison-specific settings ---
CONFIG.ct_pair      = [1, 3];   % the two CT images paired for comparison
CONFIG.reference_ct = 1;        % geometry we (blindly) reconstruct on

CONFIG.num_tr_iter = 1; 
CONFIG.conv_noise_level = .01

% Noise-only null hypothesis: instead of a single noise-only run (which jitters
% ~18-20% from the random noise draw), we run an ENSEMBLE of noise-only
% reconstructions via noise_ensemble_error_bars and use its mean +/- std gamma
% pass rate as the null band / error bar. The ensemble is cached (keyed on the
% sensor + noise + recon config, NOT the dose), so it is computed once per
% session and reused for every beam/segment. Set false to skip it entirely.
CONFIG.include_noise_only     = true;
CONFIG.noise_ensemble_minutes = 30;    % wall-time budget for a fresh ensemble

% --- Example overrides ---
% CONFIG.dose_filename         = 'dose_1194203_Session_1_reference_CT_1_B15_112.mat';
% CONFIG.reconstruction_method = 'tr';

%% ===================== PLOT SELECTION ===================================
%  Enable/disable each figure drawn at the end of this script. These only take
%  effect when CONFIG.plot_results is true (the master switch). Edit freely --
%  this PLOTS struct is local to this script.
PLOTS = struct();
PLOTS.dose_panels       = true;   % 3-view panels: the two reconstructed doses
PLOTS.convergence       = false;   % TR max-pressure / relative-change history
PLOTS.gamma_maps        = false;   % per-pair gamma + error axial maps
PLOTS.noise_only_panels = false;   % 3-view panels: blind recon vs noise-only recon

%% ===================== RESOLVE DOSE PAIR & CBCT PATHS ====================
%  Resolve the listed dose (A), then derive its counterpart (B) on the other CT
%  by swapping the _CT_k token. Each dose is paired with its own CBCT geometry.

if ~isempty(CONFIG.dose_file_override)
    dose_filepath_A = CONFIG.dose_file_override;
    [processed_dir, baseA, extA] = fileparts(dose_filepath_A);
    dose_basename_A = [baseA, extA];
else
    processed_dir   = fullfile(CONFIG.working_dir, 'RayStationFiles', ...
        CONFIG.patient_id, CONFIG.session, 'processed');
    dose_basename_A = CONFIG.dose_filename;
    dose_filepath_A = fullfile(processed_dir, dose_basename_A);
end

ct_tok = regexp(dose_basename_A, '_CT_(\d+)_', 'tokens', 'once');
if isempty(ct_tok)
    error('run_standalone_comparison:NoCTtoken', ...
        'Dose filename "%s" has no _CT_k token; cannot locate its counterpart.', dose_basename_A);
end
ct_a = str2double(ct_tok{1});
if ~ismember(ct_a, CONFIG.ct_pair)
    error('run_standalone_comparison:CTnotInPair', ...
        'Dose CT index %d is not in CONFIG.ct_pair = [%s].', ct_a, num2str(CONFIG.ct_pair));
end
ct_b_all = CONFIG.ct_pair(CONFIG.ct_pair ~= ct_a);
ct_b     = ct_b_all(1);
ct_list  = [ct_a, ct_b];

if ~ismember(CONFIG.reference_ct, CONFIG.ct_pair)
    error('run_standalone_comparison:BadReferenceCT', ...
        'CONFIG.reference_ct = %d must be one of CONFIG.ct_pair = [%s].', ...
        CONFIG.reference_ct, num2str(CONFIG.ct_pair));
end
recon_cbct_filepath = fullfile(processed_dir, sprintf('CBCT%d_resampled.mat', CONFIG.reference_ct));

dose_basename_B = regexprep(dose_basename_A, '_CT_\d+_', sprintf('_CT_%d_', ct_b));
dose_filepath_B = fullfile(processed_dir, dose_basename_B);

if ~isempty(CONFIG.cbct_file_override)
    cbct_filepath_A = CONFIG.cbct_file_override;
else
    cbct_filepath_A = fullfile(processed_dir, sprintf('CBCT%d_resampled.mat', ct_a));
end
cbct_filepath_B = fullfile(processed_dir, sprintf('CBCT%d_resampled.mat', ct_b));

dose_filepath_list = {dose_filepath_A, dose_filepath_B};
cbct_filepath_list = {cbct_filepath_A, cbct_filepath_B};
label_list         = {sprintf('CT_%d (listed)', ct_a), sprintf('CT_%d (counterpart)', ct_b)};

fprintf('=========================================================\n');
fprintf('  Standalone k-Wave Comparison (wrapper)\n');
fprintf('  Patient: %s / %s\n', CONFIG.patient_id, CONFIG.session);
fprintf('  Dose A:  %s\n', dose_filepath_list{1});
fprintf('  Dose B:  %s\n', dose_filepath_list{2});
fprintf('  Reference CT (recon geometry): CT_%d\n', CONFIG.reference_ct);
fprintf('=========================================================\n');

%% ===================== PER-DOSE RECONSTRUCTION (RUN TWICE) ===============

RES = struct([]);
for di = 1:2
    ct_this = ct_list(di);
    blind   = (ct_this ~= CONFIG.reference_ct);

    fprintf('\n########## DOSE %d/2: %s ##########\n', di, label_list{di});
    if blind
        fprintf('   [BLIND] Forward on CT_%d, reconstruction on reference CT_%d.\n', ...
            ct_this, CONFIG.reference_ct);
    else
        fprintf('   [FULL ACCESS] Forward and reconstruction on CT_%d.\n', ct_this);
    end

    cfg = CONFIG;
    cfg.dose_file_override = dose_filepath_list{di};
    cfg.cbct_file_override = cbct_filepath_list{di};
    cfg.blind_recon        = blind;
    if blind
        cfg.recon_cbct_file_override = recon_cbct_filepath;
    end

    out = run_standalone_field(cfg);
    sr  = out.sim_results;

    % Capture the reference-CT (full-access) geometry so the noise-only null
    % ensemble can build its forward bundle once on the reconstruction geometry.
    if ~blind
        sct_ref    = out.sct;
        beam_meta  = out.beam_metadata;
        gantry_ref = out.field_dose.gantry_angle;
    end

    RES(di).recon_dose        = out.recon_dose;
    RES(di).doseGrid          = out.doseGrid;
    RES(di).spacing_mm        = out.spacing;
    RES(di).sensor_mask       = out.sensor_mask;
    RES(di).density           = out.density;
    RES(di).p0                = out.p0;
    RES(di).reconPressure     = out.reconPressure;
    RES(di).gridSize          = out.gridSize;
    RES(di).fwd_time          = getf(sr, 'forward_time_s',   0);
    RES(di).tr_time           = getf(sr, 'tr_time_s',        0);
    RES(di).conv_max_pressure = getf(sr, 'conv_max_pressure', []);
    RES(di).conv_rel_change   = getf(sr, 'conv_rel_change',   []);
    RES(di).num_iters_done    = getf(sr, 'num_iters_done',    1);
    RES(di).label             = label_list{di};
end

%% ===================== EXTRACT / MAP ====================================

recon_A    = RES(1).recon_dose;
recon_B    = RES(2).recon_dose;
spacing_mm = RES(1).spacing_mm;

if ~isequal(RES(1).gridSize, RES(2).gridSize)
    error('run_standalone_comparison:GridMismatch', ...
        'The two doses are on different grids ([%s] vs [%s]).', ...
        num2str(RES(1).gridSize), num2str(RES(2).gridSize));
end

idx1 = find(ct_list == CONFIG.reference_ct, 1);   % reference-CT (full access)
idx2 = find(ct_list ~= CONFIG.reference_ct, 1);   % other-CT (blind)
if isempty(idx2), tmp = setdiff([1,2], idx1); idx2 = tmp(1); end

dose1  = double(RES(idx1).doseGrid);   recon1 = double(RES(idx1).recon_dose);
dose2  = double(RES(idx2).doseGrid);   recon2 = double(RES(idx2).recon_dose);
label1 = RES(idx1).label;              label2 = RES(idx2).label;
sensor1 = RES(idx1).sensor_mask;       sensor2 = RES(idx2).sensor_mask;

%% ===================== NOISE-ONLY NULL HYPOTHESIS (ENSEMBLE) ============
%  Replaces the old single noise-only run. noise_ensemble_error_bars runs the
%  forward simulation ONCE (on the reference/reconstruction geometry) and then
%  redraws electronic noise + reconstructs many times, returning the mean +/- std
%  gamma pass rate vs the CT-ref truth (dose1) -- the null band / error bar. It is
%  CACHED, keyed on the sensor/noise/recon config (NOT the dose), so it is computed
%  once per session and reused for every beam/segment.
%    >>> NOTE: this cache assumes the ultrasound array sits in the SAME place for
%    every beam/segment. Reexamine -- and re-key the cache on the sensor mask --
%    if sensor placement ever becomes per-beam (see noise_ensemble_error_bars). <<<
%
%  No single noise recon volume is produced anymore, so the per-volume noise gamma
%  pairs and the noise plot below are skipped; the null band alone is used.

null_ensemble = [];
recon_noise   = [];   % kept [] -- guards the legacy per-volume noise blocks below
if CONFIG.include_noise_only
    fprintf('\n########## NOISE-ONLY NULL ENSEMBLE: CT_%d geometry ##########\n', CONFIG.reference_ct);
    null_ensemble = noise_ensemble_error_bars(CONFIG, dose1, spacing_mm, ...
        sct_ref, gantry_ref, beam_meta, 'TimeBudgetMin', CONFIG.noise_ensemble_minutes);
end

%% ===================== NORMALIZE ========================================
% Least-squares normalization: scale EVERY recon to the SAME target, the CT-ref
% truth (dose1), so all recons and the null share one common scale. (The null
% ensemble in noise_ensemble_error_bars likewise normalizes to dose1.)
if isfield(CONFIG, 'normalize') && CONFIG.normalize
    g1 = least_squares_gain(dose1, recon1);
    g2 = least_squares_gain(dose1, recon2);
    recon1 = recon1 * g1;
    recon2 = recon2 * g2;
    recon_A = recon_A * least_squares_gain(dose1, recon_A);
    recon_B = recon_B * least_squares_gain(dose1, recon_B);
    if ~isempty(recon_noise)
        recon_noise = recon_noise * least_squares_gain(dose1, recon_noise);
    end
    fprintf('\n[NORM] LSQ gains (all to CT %d truth): recon(%s)=%.4g, recon(%s)=%.4g\n', ...
        CONFIG.reference_ct, label1, g1, label2, g2);
end

%% ========================= GAMMA ANALYSIS ==============================
%  3%/3mm comparisons in two groups. Each row:
%  {name, group, reference(truth), target, title, sensor}.
%    HEADLINE    -- everything referenced to the CT ref truth:
%                     CT ref truth  vs  CT blind truth
%                     CT ref truth  vs  CT ref   recon
%                     CT ref truth  vs  CT blind recon
%                     CT ref truth  vs  noise-only recon
%    SUBHEADLINE -- everything referenced to the CT blind truth:
%                     CT blind truth vs  CT blind recon
%                     CT blind truth vs  noise-only recon
%  ct_ref = full-access recon geometry (dose1/recon1); ct_blind = blind
%  geometry (dose2/recon2). recon_noise is on the blind geometry.

ct_ref   = CONFIG.reference_ct;
blind_ids = ct_list(ct_list ~= ct_ref);
ct_blind  = blind_ids(1);

gamma_criteria = {3, 3, '3%/3mm'};

% --- Headline: referenced to CT ref truth ---
gamma_pairs = {
    'headline_truth_vs_truth',      'Headline', dose1, dose2, ...
        sprintf('CT %d truth  vs  CT %d truth', ct_ref, ct_blind), sensor1;
    'headline_truth_vs_refrecon',   'Headline', dose1, recon1, ...
        sprintf('CT %d truth  vs  CT %d recon', ct_ref, ct_ref),   sensor1;
    'headline_truth_vs_blindrecon', 'Headline', dose1, recon2, ...
        sprintf('CT %d truth  vs  CT %d recon', ct_ref, ct_blind), sensor2;
};
if ~isempty(recon_noise)
    gamma_pairs = [gamma_pairs; {
        'headline_truth_vs_noise',  'Headline', dose1, recon_noise, ...
            sprintf('CT %d truth  vs  noise-only recon', ct_ref),  sensor2;
    }];
end

% --- Subheadline: referenced to CT blind truth ---
gamma_pairs = [gamma_pairs; {
    'sub_truth_vs_blindrecon',      'Subheadline', dose2, recon2, ...
        sprintf('CT %d truth  vs  CT %d recon', ct_blind, ct_blind), sensor2;
}];
if ~isempty(recon_noise)
    gamma_pairs = [gamma_pairs; {
        'sub_truth_vs_noise',       'Subheadline', dose2, recon_noise, ...
            sprintf('CT %d truth  vs  noise-only recon', ct_blind), sensor2;
    }];
end

gamma_results = struct();
if exist('CalcGamma', 'file') == 2
    for pp = 1:size(gamma_pairs, 1)
        pair_name  = gamma_pairs{pp, 1};
        ref_vol    = gamma_pairs{pp, 3};
        tgt_vol    = gamma_pairs{pp, 4};
        fprintf('\n[Gamma] %s\n', gamma_pairs{pp, 5});
        gr = compute_gamma(ref_vol, tgt_vol, spacing_mm, 'Criteria', gamma_criteria);
        gr.name        = pair_name;
        gr.group       = gamma_pairs{pp, 2};
        gr.title       = gamma_pairs{pp, 5};
        gr.ref         = ref_vol;
        gr.tgt         = tgt_vol;
        gr.sensor_mask = gamma_pairs{pp, 6};
        gamma_results.(pair_name) = gr;
    end

    fprintf('\n  ===== Gamma Pass Rates (3%%/3mm, 10%% cutoff on reference) =====\n');
    pair_fn = fieldnames(gamma_results);
    for grp = {'Headline', 'Subheadline'}
        fprintf('\n  %s:\n', grp{1});
        for pp = 1:numel(pair_fn)
            gr = gamma_results.(pair_fn{pp});
            if strcmp(gr.group, grp{1})
                fprintf('    %-36s %8.2f%%\n', gr.title, gr.pass_rates(1));
            end
        end
    end
    if ~isempty(null_ensemble)
        if null_ensemble.from_cache, cache_tag = ', cached'; else, cache_tag = ''; end
        fprintf('\n  Noise-only null (ensemble, ref = CT %d truth):\n', ct_ref);
        fprintf('    %-36s %8.2f%% +/- %.2f%%  (n=%d%s)\n', 'noise-only null band', ...
            null_ensemble.mean_pass_rate, null_ensemble.std_pass_rate, ...
            null_ensemble.num_samples, cache_tag);
        gamma_results.null_ensemble = null_ensemble;   % mean/std null band for error bars
    end
else
    warning('CalcGamma not found. Skipping gamma analysis.');
    gamma_results = [];
end

%% ========================= SAVE RESULTS ================================

if CONFIG.save_results
    output_fname = sprintf('standalone_comparison_%s_%s_%s_%d.mat', ...
        CONFIG.reconstruction_method, CONFIG.sensor_placement_method, ...
        CONFIG.gruneisen_method, CONFIG.num_time_reversal_iter);

    results = struct();
    results.recon_dose_A    = recon_A;
    results.recon_dose_B    = recon_B;
    results.label_A         = RES(1).label;
    results.label_B         = RES(2).label;
    results.original_dose_A = RES(1).doseGrid;
    results.original_dose_B = RES(2).doseGrid;
    results.p0_A            = RES(1).p0;
    results.p0_B            = RES(2).p0;
    results.reconPressure_A = RES(1).reconPressure;
    results.reconPressure_B = RES(2).reconPressure;
    results.sensor_mask_A   = RES(1).sensor_mask;
    results.sensor_mask_B   = RES(2).sensor_mask;
    if ~isempty(recon_noise)
        results.recon_dose_noise = recon_noise;   % noise-only control (blind geometry)
    end
    if ~isempty(null_ensemble)
        results.null_ensemble = null_ensemble;    % mean/std/pass_rates null band
    end
    results.config          = CONFIG;
    results.spacing_mm      = spacing_mm;
    results.grid_size       = RES(1).gridSize;
    results.fwd_time_sec    = [RES(1).fwd_time, RES(2).fwd_time];
    results.tr_time_sec     = [RES(1).tr_time,  RES(2).tr_time];
    if ~isempty(gamma_results)
        results.gamma = gamma_results;
    end
    save(output_fname, '-struct', 'results', '-v7.3');
    fprintf('\nResults saved to: %s\n', output_fname);
end

%% ========================= POST-SIMULATION PLOTS =======================

%  Each enabled config gets ONE window whose plots are collected as labeled
%  tabs (uitab) of a single uitabgroup, so related views stay together.

if CONFIG.plot_results
    if PLOTS.dose_panels
        figDose = figure('Name', 'Dose Panels', 'Color', 'w', ...
            'NumberTitle', 'off', 'Position', [50, 50, 1380, 720]);
        tgDose  = uitabgroup(figDose);
        tabDose = uitab(tgDose, 'Title', 'Reconstructed Doses');
        plot_dose_panels(recon_A, recon_B, RES(1).sensor_mask, RES(1).density, spacing_mm, ...
            'Dose Comparison: Two Reconstructed Doses', CONFIG.viz_smooth_sigma, ...
            {RES(1).label, RES(2).label}, tabDose);
    end

    if PLOTS.convergence && ...
            any(strcmpi(CONFIG.reconstruction_method, {'tr', 'hybrid'})) && ~isempty(RES(1).conv_max_pressure)
        figConv = figure('Name', 'Convergence', 'Color', 'w', ...
            'NumberTitle', 'off', 'Position', [150, 520, 760, 430]);
        tgConv  = uitabgroup(figConv);
        tabConv = uitab(tgConv, 'Title', 'p_0 Convergence');
        plot_convergence_history(RES(1).conv_max_pressure, RES(1).conv_rel_change, ...
            RES(1).num_iters_done, CONFIG.convergence_tol, max(RES(1).p0(:)), tabConv);
    end

    if PLOTS.gamma_maps && ~isempty(gamma_results) && isstruct(gamma_results)
        figGamma = figure('Name', 'Gamma Maps', 'Color', 'w', ...
            'NumberTitle', 'off', 'Position', [50, 100, 1450, 400]);
        tgGamma  = uitabgroup(figGamma);
        pair_fn  = fieldnames(gamma_results);
        for k = 1:numel(pair_fn)
            gr = gamma_results.(pair_fn{k});
            [~, max_dose_idx] = max(gr.ref(:));
            [~, ~, cz_gamma]  = ind2sub(size(gr.ref), max_dose_idx);
            tabGamma = uitab(tgGamma, 'Title', gr.title);
            plot_gamma_and_error_axial(gr, gr.ref, gr.tgt, gr.sensor_mask, cz_gamma, gr.title, tabGamma);
        end
    end

    if PLOTS.noise_only_panels && ~isempty(recon_noise)
        figNoise = figure('Name', 'Noise-only Control', 'Color', 'w', ...
            'NumberTitle', 'off', 'Position', [50, 50, 1380, 720]);
        tgNoise  = uitabgroup(figNoise);
        tabNoise = uitab(tgNoise, 'Title', sprintf('%s recon vs noise-only', label2));
        plot_dose_panels(recon2, recon_noise, sensor2, RES(idx2).density, spacing_mm, ...
            sprintf('Noise-only control: %s recon vs noise-only recon', label2), ...
            CONFIG.viz_smooth_sigma, {sprintf('%s recon', label2), 'Noise-only recon'}, tabNoise);
    end
end

fprintf('\nStandalone comparison complete.\n');


%% ========================= LOCAL HELPER =================================

function v = getf(s, f, d)
%GETF Struct field with default.
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
