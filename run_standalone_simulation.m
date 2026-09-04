%% =========================================================================
%  RUN_STANDALONE_SIMULATION.m   (wrapper)
%  Standalone k-Wave Photoacoustic Forward + Time-Reversal Simulation
%
%  Thin driver: it sets up a CONFIG (get_default_config + overrides), runs ONE
%  field through the shared simulation engine via run_standalone_field
%  (load dose+CBCT -> build medium -> run_single_field_simulation), then does
%  the standalone analysis (normalize -> gamma -> log -> save -> plots) using
%  the shared helpers in utils/. All physics lives in run_single_field_simulation.
%  =========================================================================

clear; clc; close all;

% Ensure the moved helper functions in utils/ are on the path (run from root).
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'utils')));
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'pipeline')));  % step* functions live here

%% ========================= CONFIGURATION ================================
%  Defaults come from get_default_config (the single source of truth).
%  Override only what you want to change for this run below.

CONFIG = get_default_config();
CONFIG.sensor_placement_method = 'spherical';
CONFIG.reconstruction_method = 'ubp'; 
CONFIG.gruneisen_method = 'uniform';


% --- Example overrides (uncomment / edit as needed) ---
% CONFIG.dose_filename          = 'dose_1194203_Session_1_reference_CT_1_B15_112.mat';
% CONFIG.cbct_filename          = 'CBCT1_resampled.mat';
% CONFIG.reconstruction_method  = 'tr';
% CONFIG.num_time_reversal_iter = 10;
% CONFIG.downscale_factor       = 1;

%% ========================= RUN SIMULATION ===============================

fprintf('=========================================================\n');
fprintf('  Standalone k-Wave Photoacoustic Simulation (wrapper)\n');
fprintf('=========================================================\n');
fprintf('  Patient:      %s / %s\n', CONFIG.patient_id, CONFIG.session);
fprintf('  Dose file:    %s\n', CONFIG.dose_filename);
fprintf('  Sensor:       %s\n', CONFIG.sensor_placement_method);
fprintf('  Tissue model: %s\n', CONFIG.gruneisen_method);
fprintf('  Recon method: %s\n', CONFIG.reconstruction_method);
fprintf('=========================================================\n\n');

out = run_standalone_field(CONFIG);

recon_dose    = out.recon_dose;
doseGrid      = out.doseGrid;
reconPressure = out.reconPressure;
p0            = out.p0;
sensor_mask   = out.sensor_mask;
density       = out.density;
spacing_mm    = out.spacing;
gridSize      = out.gridSize;
sr            = out.sim_results;

conv_max_pressure = getf(sr, 'conv_max_pressure', []);
conv_rel_change   = getf(sr, 'conv_rel_change',   []);
num_iters_done    = getf(sr, 'num_iters_done',    1);
fwd_time          = getf(sr, 'forward_time_s',    0);
tr_time           = getf(sr, 'tr_time_s',         0);

%% ========================= NORMALIZE DOSES ==============================
%  'max' : divide both doses by their own maxima.
%  'lsq' : scale the recon by the least-squares gain to the (unscaled) truth.

if isfield(CONFIG, 'normalize') && CONFIG.normalize
    norm_method = 'max';
    if isfield(CONFIG, 'normalize_method') && ~isempty(CONFIG.normalize_method)
        norm_method = lower(CONFIG.normalize_method);
    end
    switch norm_method
        case 'max'
            dg_max = max(doseGrid(:));  rd_max = max(recon_dose(:));
            if dg_max > 0, doseGrid   = doseGrid   / dg_max; end
            if rd_max > 0, recon_dose = recon_dose / rd_max; end
            fprintf('[NORM] max-normalized both doses (dg_max=%.4f, rd_max=%.4f).\n', dg_max, rd_max);
        case 'lsq'
            g = least_squares_gain(doseGrid, recon_dose);
            recon_dose = recon_dose * g;
            fprintf('[NORM] least-squares gain (recon -> truth): g = %.4f\n', g);
        otherwise
            error('run_standalone_simulation:BadNormalizeMethod', ...
                'Unknown CONFIG.normalize_method "%s" (use ''max'' or ''lsq'').', norm_method);
    end
end

doseThreshold = 0.01 * max(doseGrid(:));

%% ========================= RESULTS SUMMARY ==============================

fprintf('\n========= RESULTS =========\n');
fprintf('  Original dose:      [%.6f, %.4f]\n', min(doseGrid(:)),   max(doseGrid(:)));
fprintf('  Reconstructed dose: [%.6f, %.4f]\n', min(recon_dose(:)), max(recon_dose(:)));
dose_region = doseGrid > doseThreshold;
if any(dose_region(:))
    abs_error = abs(recon_dose(dose_region) - doseGrid(dose_region));
    rel_error = abs_error ./ max(doseGrid(dose_region), 1e-10) * 100;
    fprintf('  Mean abs err: %.6f\n', mean(abs_error));
    fprintf('  Mean rel err: %.2f%%\n', mean(rel_error));
    fprintf('  Max  rel err: %.2f%%\n', max(rel_error));
end
fprintf('===========================\n');

%% ========================= GAMMA ANALYSIS ==============================

gamma_results = compute_gamma(doseGrid, recon_dose, spacing_mm);
append_gamma_log(CONFIG, gamma_results, CONFIG.dose_filename);

%% ========================= SAVE RESULTS ================================

if CONFIG.save_results
    output_fname = sprintf('standalone_results_%s_%s_%s_%d.mat', ...
        CONFIG.reconstruction_method, CONFIG.sensor_placement_method, ...
        CONFIG.gruneisen_method, CONFIG.num_time_reversal_iter);

    results = struct();
    results.recon_dose    = recon_dose;
    results.original_dose = doseGrid;
    results.p0            = p0;
    results.reconPressure = reconPressure;
    results.sensor_mask   = sensor_mask;
    results.config        = CONFIG;
    results.spacing_mm    = spacing_mm;
    results.grid_size     = gridSize;
    results.fwd_time_sec  = fwd_time;
    results.tr_time_sec   = tr_time;
    if ~isempty(gamma_results)
        results.gamma = gamma_results;
    end
    save(output_fname, '-struct', 'results', '-v7.3');
    fprintf('\nResults saved to: %s\n', output_fname);
end

%% ========================= POST-SIMULATION PLOTS =======================

if CONFIG.plot_results
    dose_mask_vis = double(doseGrid) >= 0.10 * max(double(doseGrid(:)));
    plot_sensor_dose_planes(dose_mask_vis, sensor_mask, spacing_mm, density, CONFIG);

    plot_dose_panels(doseGrid, recon_dose, sensor_mask, density, spacing_mm, ...
        'Dose Comparison: Original vs Reconstructed', CONFIG.viz_smooth_sigma);

    if any(strcmpi(CONFIG.reconstruction_method, {'tr', 'hybrid'})) && ~isempty(conv_max_pressure)
        plot_convergence_history(conv_max_pressure, conv_rel_change, ...
            num_iters_done, CONFIG.convergence_tol, max(p0(:)));
    end

    if ~isempty(gamma_results) && isfield(gamma_results, 'maps')
        [~, max_dose_idx] = max(doseGrid(:));
        [~, ~, cz_gamma]  = ind2sub(gridSize, max_dose_idx);
        plot_gamma_and_error_axial(gamma_results, doseGrid, recon_dose, sensor_mask, cz_gamma);
    end
end

fprintf('\nStandalone simulation complete.\n');


%% ========================= LOCAL HELPER =================================

function v = getf(s, f, d)
%GETF Struct field with default.
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
