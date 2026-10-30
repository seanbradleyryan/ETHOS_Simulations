%% =========================================================================
%  DIAGNOSTIC_RECON_COMPARISON.m
%  Reconstruction-method comparison on one field
%
%  Runs like run_standalone_simulation, double loop: once per method in
%  CONFIG.recon_methods and, for the adjoint_tv* methods, once per TV weight
%  in CONFIG.tv_lambdas (other methods run once). Each run goes through
%  run_standalone_field (load dose+CBCT -> build medium ->
%  run_single_field_simulation). Every run sees the SAME measured data (same
%  dose, sensor and noise draw) and runs the SAME number of iterations (no
%  early stop). Each recon is then normalized and gamma-compared to the truth
%  dose, the pass rates are tabulated, and the best lambda per TV method is
%  printed.
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
CONFIG.recon_methods  = {
    'tr', 
    %'adjoint', 
    'adjoint_tv'};   % any of 'tr','das','ubp','adjoint','adjoint_tv','adjoint_tv_tr'
CONFIG.num_iter       = 1;     % iterations for every iterative method
CONFIG.rng_seed       = 42;     % same electronic-noise draw for every method
CONFIG.gamma_criteria = {3, 3, '3%/3mm'};
CONFIG.Nt_scaling = 3; 
CONFIG.conv_noise_level = .125; 
CONFIG.conv_use_deconv = false; 
CONFIG.gruneisen_method = 'continuous_1'


% CONFIG.sensor_placement_method = 'box'; 
% CONFIG.gruneisen_method = 'uniform'; 
% 50, 75

%CONFIG.num_iter = 1; 
%CONFIG.Nt_scaling = 3; 
%33, 30.8, 30

% CONFIG.adjoint_tv_lambda = 0; 




% Same iteration count for every method: tol 0 disables the early stop.
CONFIG.num_time_reversal_iter = CONFIG.num_iter;
CONFIG.adjoint_num_iter       = CONFIG.num_iter;

CONFIG.convergence_tol        = 0;

% --- adjoint_tv knobs (see run_single_field_simulation) ---
CONFIG.tv_lambdas            = [0, 0.001, 0.01, 0.1];   % adjoint_tv_lambda sweep (adjoint_tv* methods only)
CONFIG.adjoint_tv_force_power_recalc = false;   % true = ignore processed/simulated_doses/adjoint_power_cache_<hash>.mat, recompute + overwrite
CONFIG.adjoint_tv_power_iter = 5;
CONFIG.adjoint_tv_inner_iter = CONFIG.num_iter; 

% --- Example overrides ---
% CONFIG.dose_filename           = 'dose_1194203_Session_1_reference_CT_1_B15_112.mat';
% CONFIG.sensor_placement_method = 'spherical';

% Gamma runs at the end; fail now rather than after the k-Wave runs.
if exist('CalcGamma', 'file') ~= 2
    error('diagnostic_recon_comparison:NoCalcGamma', 'CalcGamma must be on the MATLAB path.');
end

numMethods = numel(CONFIG.recon_methods);
isTV       = startsWith(CONFIG.recon_methods, 'adjoint_tv');   % only these use adjoint_tv_lambda
numRuns    = sum(~isTV) + sum(isTV) * numel(CONFIG.tv_lambdas);

fprintf('=========================================================\n');
fprintf('  Reconstruction method comparison\n');
fprintf('=========================================================\n');
fprintf('  Patient:      %s / %s\n', CONFIG.patient_id, CONFIG.session);
fprintf('  Dose file:    %s\n', CONFIG.dose_filename);
fprintf('  Sensor:       %s\n', CONFIG.sensor_placement_method);
fprintf('  Tissue model: %s\n', CONFIG.gruneisen_method);
fprintf('  Methods:      %s (%d iterations each)\n', strjoin(CONFIG.recon_methods, ', '), CONFIG.num_iter);
fprintf('  TV lambdas:   %s (%d runs total)\n', mat2str(CONFIG.tv_lambdas), numRuns);
fprintf('=========================================================\n\n');

%% ========================= RUN EACH METHOD x LAMBDA =====================

RES = struct([]);
ri  = 0;   % run index into RES
for mi = 1:numMethods
    method = CONFIG.recon_methods{mi};
    if isTV(mi)
        lambdas = CONFIG.tv_lambdas;
    else
        lambdas = NaN;   % lambda unused: run once
    end

    for li = 1:numel(lambdas)
        ri  = ri + 1;
        cfg = CONFIG;
        cfg.reconstruction_method = method;
        if isnan(lambdas(li))
            label = method;
        else
            cfg.adjoint_tv_lambda = lambdas(li);
            label = sprintf('%s, lambda=%g', method, lambdas(li));
        end
        fprintf('\n########## RUN %d/%d: %s ##########\n', ri, numRuns, label);

        rng(CONFIG.rng_seed);   % identical noise -> identical measured data
        out = run_standalone_field(cfg);
        sr  = out.sim_results;

        RES(ri).method         = method;
        RES(ri).tv_lambda      = lambdas(li);
        RES(ri).label          = label;
        RES(ri).recon_dose     = out.recon_dose;
        RES(ri).tr_time        = getf(sr, 'tr_time_s', 0);
        RES(ri).num_iters_done = getf(sr, 'num_iters_done', 1);
    end
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

for ri = 1:numRuns
    recon = RES(ri).recon_dose;
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
    RES(ri).recon_dose = recon;

    fprintf('\n[Gamma] %s\n', RES(ri).label);
    RES(ri).gamma = compute_gamma(doseGrid, recon, spacing_mm, 'Criteria', CONFIG.gamma_criteria);
end

%% ========================= SUMMARY TABLE ================================

critLabels = CONFIG.gamma_criteria(:, 3)';
passRates  = zeros(numRuns, numel(critLabels));
fprintf('\n========= GAMMA PASS RATES (normalization: %s) =========\n', norm_method);
fprintf('  %-28s %6s %12s', 'Method', 'Iters', 'Recon (s)');
fprintf(' %10s', critLabels{:});
fprintf('\n');
for ri = 1:numRuns
    passRates(ri, :) = RES(ri).gamma.pass_rates(:)';
    fprintf('  %-28s %6d %12.1f', RES(ri).label, RES(ri).num_iters_done, RES(ri).tr_time);
    fprintf(' %9.2f%%', passRates(ri, :));
    fprintf('\n');
end

% Best lambda per TV method, ranked by the first gamma criterion.
for mi = find(isTV(:)')
    runIdx = find(strcmp({RES.method}, CONFIG.recon_methods{mi}));
    [bestRate, b] = max(passRates(runIdx, 1));
    fprintf('  Best lambda for %s: %g (%s pass rate %.2f%%)\n', ...
        CONFIG.recon_methods{mi}, RES(runIdx(b)).tv_lambda, critLabels{1}, bestRate);
end
fprintf('=========================================================\n');

%% ========================= SAVE RESULTS ================================

if CONFIG.save_results
    output_fname = sprintf('recon_comparison_%s_%s_%d.mat', ...
        CONFIG.sensor_placement_method, CONFIG.gruneisen_method, CONFIG.num_iter);
    results = struct();
    results.methods       = CONFIG.recon_methods;
    results.tv_lambdas    = CONFIG.tv_lambdas;
    results.labels        = {RES.label};   % one per row of pass_rates
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
    % Pass rate per run (method x lambda), one bar per criterion; red line = 90% pass rate.
    runLabels = {RES.label};
    figure('Name', 'Gamma Pass Rates', 'Color', 'w', 'NumberTitle', 'off');
    bar(categorical(runLabels, runLabels), passRates);
    yline(90, 'r-', 'LineWidth', 1.5);
    ylim([0, 100]); grid on;
    ylabel('Gamma pass rate (%)');
    legend([critLabels, {'90%'}], 'Location', 'southoutside', 'Orientation', 'horizontal');
    title(sprintf('Truth vs recon, %d iterations each', CONFIG.num_iter));

    % Dose panels and gamma maps: one window each, one tab per run.
    [~, max_dose_idx] = max(doseGrid(:));
    [~, ~, cz_gamma]  = ind2sub(size(doseGrid), max_dose_idx);
    figDose  = figure('Name', 'Dose Panels', 'Color', 'w', 'NumberTitle', 'off', 'Position', [50, 50, 1380, 720]);
    tgDose   = uitabgroup(figDose);
    figGamma = figure('Name', 'Gamma Maps', 'Color', 'w', 'NumberTitle', 'off', 'Position', [50, 100, 1450, 400]);
    tgGamma  = uitabgroup(figGamma);
    for ri = 1:numRuns
        plot_dose_panels(doseGrid, RES(ri).recon_dose, sensor_mask, density, spacing_mm, ...
            sprintf('Truth vs %s', RES(ri).label), CONFIG.viz_smooth_sigma, ...
            {'Truth', RES(ri).label}, uitab(tgDose, 'Title', RES(ri).label));
        plot_gamma_and_error_axial(RES(ri).gamma, doseGrid, RES(ri).recon_dose, sensor_mask, ...
            cz_gamma, RES(ri).label, uitab(tgGamma, 'Title', RES(ri).label));
    end
end

fprintf('\nReconstruction comparison complete.\n');


%% ========================= LOCAL HELPER =================================

function v = getf(s, f, d)
%GETF Struct field with default.
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
