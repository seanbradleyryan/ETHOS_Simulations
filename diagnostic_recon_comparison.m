%% =========================================================================
%  DIAGNOSTIC_RECON_COMPARISON.m
%  Reconstruction-method comparison on one field
%
%  Runs like run_standalone_simulation, once per method in
%  CONFIG.recon_methods: each run goes through run_standalone_field (load
%  dose+CBCT -> build medium -> run_single_field_simulation). Every method
%  sees the SAME measured data (same dose, sensor and noise draw) and runs the
%  SAME number of iterations (no early stop). Each recon is then normalized
%  and gamma-compared to the truth dose, and the pass rates are tabulated.
%
%  Note: run_standalone_field re-runs the forward simulation for every
%  method; rng(CONFIG.rng_seed) before each run makes the noise identical.
%  =========================================================================

clear; clc; close all;

% Ensure the moved helper functions in utils/ are on the path (run from root).
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'utils')));
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'pipeline')));  % step* functions live here

%% ========================= CONFIGURATION ================================

CONFIG = get_default_config();

% --- Comparison settings ---
CONFIG.recon_methods  = {'tr', 'adjoint', 'adjoint_tv'};   % any of 'tr','das','ubp','adjoint','adjoint_tv'
CONFIG.num_iter       = 5;     % iterations for every iterative method
CONFIG.rng_seed       = 1;     % same electronic-noise draw for every method
CONFIG.gamma_criteria = {3, 3, '3%/3mm'; 5, 5, '5%/5mm'};

% Same iteration count for every method: tol 0 disables the early stop.
CONFIG.num_time_reversal_iter = CONFIG.num_iter;
CONFIG.adjoint_num_iter       = CONFIG.num_iter;
CONFIG.convergence_tol        = 0;

% --- adjoint_tv knobs (see run_single_field_simulation) ---
CONFIG.adjoint_tv_lambda     = 0.01;
CONFIG.adjoint_tv_power_iter = 5;
CONFIG.adjoint_tv_inner_iter = 10;

% --- Example overrides ---
% CONFIG.dose_filename           = 'dose_1194203_Session_1_reference_CT_1_B15_112.mat';
% CONFIG.sensor_placement_method = 'spherical';

% Gamma runs at the end; fail now rather than after the k-Wave runs.
if exist('CalcGamma', 'file') ~= 2
    error('diagnostic_recon_comparison:NoCalcGamma', 'CalcGamma must be on the MATLAB path.');
end

numMethods = numel(CONFIG.recon_methods);

fprintf('=========================================================\n');
fprintf('  Reconstruction method comparison\n');
fprintf('=========================================================\n');
fprintf('  Patient:      %s / %s\n', CONFIG.patient_id, CONFIG.session);
fprintf('  Dose file:    %s\n', CONFIG.dose_filename);
fprintf('  Sensor:       %s\n', CONFIG.sensor_placement_method);
fprintf('  Tissue model: %s\n', CONFIG.gruneisen_method);
fprintf('  Methods:      %s (%d iterations each)\n', strjoin(CONFIG.recon_methods, ', '), CONFIG.num_iter);
fprintf('=========================================================\n\n');

%% ========================= RUN EACH METHOD ==============================

RES = struct([]);
for mi = 1:numMethods
    method = CONFIG.recon_methods{mi};
    fprintf('\n########## METHOD %d/%d: %s ##########\n', mi, numMethods, method);

    cfg = CONFIG;
    cfg.reconstruction_method = method;
    rng(CONFIG.rng_seed);   % identical noise -> identical measured data
    out = run_standalone_field(cfg);
    sr  = out.sim_results;

    RES(mi).method         = method;
    RES(mi).recon_dose     = out.recon_dose;
    RES(mi).tr_time        = getf(sr, 'tr_time_s', 0);
    RES(mi).num_iters_done = getf(sr, 'num_iters_done', 1);
end

% Truth and geometry are the same for every method: take them from the last run.
doseGrid    = out.doseGrid;
spacing_mm  = out.spacing;
sensor_mask = out.sensor_mask;
density     = out.density;

%% ========================= NORMALIZE + GAMMA ============================
%  'max' : divide both doses by their own maxima.
%  'lsq' : scale each recon by the least-squares gain to the (unscaled) truth.

norm_method = 'none';
if CONFIG.normalize
    norm_method = lower(CONFIG.normalize_method);
end
if strcmp(norm_method, 'max')
    doseGrid = doseGrid / max(doseGrid(:));
end

for mi = 1:numMethods
    recon = RES(mi).recon_dose;
    switch norm_method
        case 'max'
            if max(recon(:)) > 0, recon = recon / max(recon(:)); end
        case 'lsq'
            recon = recon * least_squares_gain(doseGrid, recon);
        case 'none'
        otherwise
            error('diagnostic_recon_comparison:BadNormalizeMethod', ...
                'Unknown CONFIG.normalize_method "%s" (use ''max'' or ''lsq'').', norm_method);
    end
    RES(mi).recon_dose = recon;

    fprintf('\n[Gamma] %s\n', RES(mi).method);
    RES(mi).gamma = compute_gamma(doseGrid, recon, spacing_mm, 'Criteria', CONFIG.gamma_criteria);
end

%% ========================= SUMMARY TABLE ================================

critLabels = CONFIG.gamma_criteria(:, 3)';
passRates  = zeros(numMethods, numel(critLabels));
fprintf('\n========= GAMMA PASS RATES (normalization: %s) =========\n', norm_method);
fprintf('  %-12s %6s %12s', 'Method', 'Iters', 'Recon (s)');
fprintf(' %10s', critLabels{:});
fprintf('\n');
for mi = 1:numMethods
    passRates(mi, :) = RES(mi).gamma.pass_rates(:)';
    fprintf('  %-12s %6d %12.1f', RES(mi).method, RES(mi).num_iters_done, RES(mi).tr_time);
    fprintf(' %9.2f%%', passRates(mi, :));
    fprintf('\n');
end
fprintf('=========================================================\n');

%% ========================= SAVE RESULTS ================================

if CONFIG.save_results
    output_fname = sprintf('recon_comparison_%s_%s_%d.mat', ...
        CONFIG.sensor_placement_method, CONFIG.gruneisen_method, CONFIG.num_iter);
    results = struct();
    results.methods       = CONFIG.recon_methods;
    results.res           = RES;
    results.pass_rates    = passRates;
    results.original_dose = doseGrid;
    results.sensor_mask   = sensor_mask;
    results.spacing_mm    = spacing_mm;
    results.config        = CONFIG;
    save(output_fname, '-struct', 'results', '-v7.3');
    fprintf('\nResults saved to: %s\n', output_fname);
end

%% ========================= PLOTS =======================================

if CONFIG.plot_results
    % Pass rate per method, one bar per criterion; red line = 90% pass rate.
    figure('Name', 'Gamma Pass Rates', 'Color', 'w', 'NumberTitle', 'off');
    bar(categorical(CONFIG.recon_methods, CONFIG.recon_methods), passRates);
    yline(90, 'r-', 'LineWidth', 1.5);
    ylim([0, 100]); grid on;
    ylabel('Gamma pass rate (%)');
    legend([critLabels, {'90%'}], 'Location', 'southoutside', 'Orientation', 'horizontal');
    title(sprintf('Truth vs recon, %d iterations each', CONFIG.num_iter));

    % Dose panels and gamma maps: one window each, one tab per method.
    [~, max_dose_idx] = max(doseGrid(:));
    [~, ~, cz_gamma]  = ind2sub(size(doseGrid), max_dose_idx);
    figDose  = figure('Name', 'Dose Panels', 'Color', 'w', 'NumberTitle', 'off', 'Position', [50, 50, 1380, 720]);
    tgDose   = uitabgroup(figDose);
    figGamma = figure('Name', 'Gamma Maps', 'Color', 'w', 'NumberTitle', 'off', 'Position', [50, 100, 1450, 400]);
    tgGamma  = uitabgroup(figGamma);
    for mi = 1:numMethods
        plot_dose_panels(doseGrid, RES(mi).recon_dose, sensor_mask, density, spacing_mm, ...
            sprintf('Truth vs %s', RES(mi).method), CONFIG.viz_smooth_sigma, ...
            {'Truth', RES(mi).method}, uitab(tgDose, 'Title', RES(mi).method));
        plot_gamma_and_error_axial(RES(mi).gamma, doseGrid, RES(mi).recon_dose, sensor_mask, ...
            cz_gamma, RES(mi).method, uitab(tgGamma, 'Title', RES(mi).method));
    end
end

fprintf('\nReconstruction comparison complete.\n');


%% ========================= LOCAL HELPER =================================

function v = getf(s, f, d)
%GETF Struct field with default.
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
