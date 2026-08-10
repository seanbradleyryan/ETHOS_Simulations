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

%% ========================= CONFIGURATION ================================

CONFIG = get_default_config();
CONFIG.mask_recon_to_dose_region     = false;    % zero recon dose outside the >1% dose mask


% --- Comparison-specific settings ---
CONFIG.ct_pair      = [1, 3];   % the two CT images paired for comparison
CONFIG.reference_ct = 1;        % geometry we (blindly) reconstruct on

CONFIG.num_tr_iter = 1; 

% --- Example overrides ---
% CONFIG.dose_filename         = 'dose_1194203_Session_1_reference_CT_1_B15_112.mat';
% CONFIG.reconstruction_method = 'tr';

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

%% ===================== EXTRACT / MAP / NORMALIZE =========================

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

% Least-squares normalization: scale each recon to its OWN-CT truth.
if isfield(CONFIG, 'normalize') && CONFIG.normalize
    g1 = least_squares_gain(dose1, recon1);
    g2 = least_squares_gain(dose2, recon2);
    recon1 = recon1 * g1;
    recon2 = recon2 * g2;
    recon_A = recon_A * least_squares_gain(double(RES(1).doseGrid), recon_A);
    recon_B = recon_B * least_squares_gain(double(RES(2).doseGrid), recon_B);
    fprintf('\n[NORM] LSQ gains: recon(%s)=%.4g, recon(%s)=%.4g\n', label1, g1, label2, g2);
end

%% ========================= GAMMA ANALYSIS ==============================
%  Four 3%/3mm comparisons. Each pair: {name, reference(truth), target, title, sensor}.

gamma_criteria = {3, 3, '3%/3mm'};
gamma_pairs = {
    'dose1_vs_dose2',  dose1,  dose2,  sprintf('%s truth  vs  %s truth', label1, label2), sensor1;
    'recon1_vs_dose1', dose1,  recon1, sprintf('%s recon  vs  %s truth', label1, label1), sensor1;
    'recon2_vs_dose2', dose2,  recon2, sprintf('%s recon  vs  %s truth', label2, label2), sensor2;
    'recon2_vs_dose1', dose1,  recon2, sprintf('%s recon  vs  %s truth', label2, label1), sensor2;
};

gamma_results = struct();
if exist('CalcGamma', 'file') == 2
    for pp = 1:size(gamma_pairs, 1)
        pair_name  = gamma_pairs{pp, 1};
        ref_vol    = gamma_pairs{pp, 2};
        tgt_vol    = gamma_pairs{pp, 3};
        fprintf('\n[Gamma] %s\n', gamma_pairs{pp, 4});
        gr = compute_gamma(ref_vol, tgt_vol, spacing_mm, 'Criteria', gamma_criteria);
        gr.name  = pair_name;
        gr.title = gamma_pairs{pp, 4};
        gr.ref   = ref_vol;
        gr.tgt   = tgt_vol;
        gr.sensor_mask = gamma_pairs{pp, 5};
        gamma_results.(pair_name) = gr;
    end

    fprintf('\n  ------ Gamma Pass Rates (10%% cutoff on reference) ------\n');
    pair_fn = fieldnames(gamma_results);
    for pp = 1:numel(pair_fn)
        gr = gamma_results.(pair_fn{pp});
        fprintf('  %-34s %8.2f%%\n', gr.title, gr.pass_rates(1));
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

if CONFIG.plot_results
    plot_dose_panels(recon_A, recon_B, RES(1).sensor_mask, RES(1).density, spacing_mm, ...
        'Dose Comparison: Two Reconstructed Doses', CONFIG.viz_smooth_sigma, ...
        {RES(1).label, RES(2).label});

    if any(strcmpi(CONFIG.reconstruction_method, {'tr', 'hybrid'})) && ~isempty(RES(1).conv_max_pressure)
        plot_convergence_history(RES(1).conv_max_pressure, RES(1).conv_rel_change, ...
            RES(1).num_iters_done, CONFIG.convergence_tol, max(RES(1).p0(:)));
    end

    if ~isempty(gamma_results) && isstruct(gamma_results)
        pair_fn = fieldnames(gamma_results);
        for k = 1:numel(pair_fn)
            gr = gamma_results.(pair_fn{k});
            [~, max_dose_idx] = max(gr.ref(:));
            [~, ~, cz_gamma]  = ind2sub(size(gr.ref), max_dose_idx);
            plot_gamma_and_error_axial(gr, gr.ref, gr.tgt, gr.sensor_mask, cz_gamma, gr.title);
        end
    end
end

fprintf('\nStandalone comparison complete.\n');


%% ========================= LOCAL HELPER =================================

function v = getf(s, f, d)
%GETF Struct field with default.
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
