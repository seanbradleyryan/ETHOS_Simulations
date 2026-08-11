%% =========================================================================
%  RUN_AIR_SOUND_SPEED_COMPARISON.m
%  Compare k-Wave photoacoustic reconstruction with air sound speed = 343 m/s
%  (physically correct, slow) vs 1480 m/s (water speed, fast) under otherwise
%  identical conditions, and gamma-analyze the difference.
%
%  Motivation: simTime = 2.5*gridDiag/minC and minC is set by the slowest
%  tissue. Real air (343 m/s) inflates Nt; using 1480 m/s raises minC (~1450)
%  and cuts Nt ~4x. This script validates that the speedup does not materially
%  change the reconstructed dose.
%
%  Modeled on run_standalone_simulation.m. Noise is disabled for both runs so
%  any gamma difference is attributable to the air-speed change alone.
%% =========================================================================

clear; clc; close all;

% Ensure the moved helper functions in utils/ are on the path (run from root).
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'utils')));

%% ========================= CONFIGURATION ================================

% Shared simulation defaults come from get_default_config (single source of
% truth). This comparison overrides only what it needs; every other field uses
% the default, which already matches this script's previous hard-coded values.
CONFIG = get_default_config();

CONFIG.num_time_reversal_iter = 10;  % more TR iterations than the default (1)

% Pulse convolution / noise DISABLED so any gamma difference between the two
% runs is attributable to the air sound-speed change alone.
CONFIG.convolution_kernel     = 0;   % 0 => frequency response only, no noise

CONFIG.output_file            = 'air_sound_speed_comparison_results.mat';

% --- Air sound speeds to sweep (m/s). First entry is the physically-correct
%     reference (343); second is the fast approximation (1480). ---
CONFIG.air_sound_speeds = [343, 1480];

%% ========================= RESOLVE FILE PATHS ============================

if ~isempty(CONFIG.dose_file_override)
    dose_filepath = CONFIG.dose_file_override;
else
    processed_dir = fullfile(CONFIG.working_dir, 'RayStationFiles', ...
        CONFIG.patient_id, CONFIG.session, 'processed');
    dose_filepath = fullfile(processed_dir, CONFIG.dose_filename);
end

if ~isempty(CONFIG.cbct_file_override)
    cbct_filepath = CONFIG.cbct_file_override;
else
    if ~exist('processed_dir', 'var')
        processed_dir = fullfile(CONFIG.working_dir, 'RayStationFiles', ...
            CONFIG.patient_id, CONFIG.session, 'processed');
    end
    cbct_filepath = fullfile(processed_dir, CONFIG.cbct_filename);
end

%% ========================= LOAD PLAN BEAM METADATA =======================

if ~isfield(CONFIG, 'beam_metadata') || isempty(CONFIG.beam_metadata)
    if exist('processed_dir', 'var')
        metadata_filepath = fullfile(processed_dir, 'metadata.mat');
    else
        metadata_filepath = '';
    end

    if ~isempty(metadata_filepath) && isfile(metadata_filepath)
        try
            md = load(metadata_filepath, 'metadata');
            if isfield(md, 'metadata') && isfield(md.metadata, 'beam_metadata') ...
                    && ~isempty(md.metadata.beam_metadata)
                CONFIG.beam_metadata = md.metadata.beam_metadata;
                fprintf('  Loaded beam_metadata for %d beams from %s\n', ...
                    length(CONFIG.beam_metadata), metadata_filepath);
            else
                fprintf('  [WARN] metadata.mat present but no beam_metadata field.\n');
            end
        catch ME
            fprintf('  [WARN] Failed to load %s: %s\n', metadata_filepath, ME.message);
        end
    else
        fprintf('  [WARN] No metadata.mat in processed_dir; sensor exclusion zone will be empty.\n');
    end
end

%% ========================= PRINT CONFIGURATION ===========================

fprintf('=========================================================\n');
fprintf('  Air Sound-Speed Comparison  (343 vs 1480 m/s)\n');
fprintf('=========================================================\n');
fprintf('  Patient:         %s / %s\n', CONFIG.patient_id, CONFIG.session);
fprintf('  Dose file:       %s\n', dose_filepath);
fprintf('  CBCT file:       %s\n', cbct_filepath);
fprintf('  Sensor:          %s\n', CONFIG.sensor_placement_method);
fprintf('  Tissue model:    %s\n', CONFIG.gruneisen_method);
fprintf('  Recon method:    %s\n', CONFIG.reconstruction_method);
fprintf('  Noise:           DISABLED (convolution_kernel = 0)\n');
fprintf('  Air speeds:      [%s] m/s\n', num2str(CONFIG.air_sound_speeds));
fprintf('=========================================================\n\n');

%% ========================= LOAD DATA ====================================

fprintf('[LOAD] Loading dose data...\n');
if ~isfile(dose_filepath)
    error('Dose file not found: %s', dose_filepath);
end
dose_data = load(dose_filepath);
dose_fields = fieldnames(dose_data);

step15_spacing_mm = [];
fd_gantry_angle   = [];

if isfield(dose_data, 'field_dose')
    fd = dose_data.field_dose;
    if ~isfield(fd, 'dose_Gy')
        error('field_dose struct missing dose_Gy field.');
    end
    if (isfield(fd, 'is_sparse') && fd.is_sparse) || issparse(fd.dose_Gy)
        if ~isfield(fd, 'dose_dims')
            error('field_dose.dose_dims missing - cannot reconstruct sparse dose.');
        end
        doseGrid = reshape(full(fd.dose_Gy), fd.dose_dims);
        fprintf('       Loaded: field_dose.dose_Gy (sparse -> [%d x %d x %d])\n', fd.dose_dims);
    else
        doseGrid = double(fd.dose_Gy);
        fprintf('       Loaded: field_dose.dose_Gy (dense)\n');
    end
    if isfield(fd, 'spacing') && ~isempty(fd.spacing)
        step15_spacing_mm = fd.spacing(:)';
    end
    if isfield(fd, 'meterset') && ~isempty(fd.meterset) && fd.meterset > 0
        if CONFIG.meterset ~= fd.meterset
            fprintf('       [INFO] Overriding CONFIG.meterset: %.2f -> %.2f MU\n', ...
                CONFIG.meterset, fd.meterset);
            CONFIG.meterset = fd.meterset;
        end
    end
    if isfield(fd, 'gantry_angle')
        fd_gantry_angle = fd.gantry_angle;
        fprintf('       Gantry angle: %.1f deg\n', fd_gantry_angle);
    end
elseif isfield(dose_data, 'total_rs_dose_sparse')
    if ~isfield(dose_data, 'total_rs_dose_dims')
        error('total_rs_dose_dims missing - cannot reconstruct sparse total dose.');
    end
    doseGrid = reshape(full(dose_data.total_rs_dose_sparse), dose_data.total_rs_dose_dims);
    fprintf('       Loaded: total_rs_dose_sparse ([%d x %d x %d])\n', dose_data.total_rs_dose_dims);
elseif isfield(dose_data, 'total_rs_dose')
    doseGrid = dose_data.total_rs_dose;
elseif isfield(dose_data, 'dose_Gy')
    doseGrid = dose_data.dose_Gy;
elseif length(dose_fields) == 1
    doseGrid = dose_data.(dose_fields{1});
else
    error('Cannot auto-detect dose variable. Fields found: %s', strjoin(dose_fields, ', '));
end
doseGrid = double(doseGrid);

if ~isnumeric(doseGrid) || ndims(doseGrid) ~= 3
    error('Dose data must be a 3D numeric array.');
end

gridSize = size(doseGrid);
Nx = gridSize(1); Ny = gridSize(2); Nz = gridSize(3);
fprintf('       Grid size: [%d x %d x %d]\n', Nx, Ny, Nz);
fprintf('       Dose range: [%.6f, %.4f] Gy\n', min(doseGrid(:)), max(doseGrid(:)));

fprintf('[LOAD] Loading CBCT data...\n');
if ~isfile(cbct_filepath)
    error('CBCT file not found: %s', cbct_filepath);
end
cbct_data = load(cbct_filepath);
if isfield(cbct_data, 'CBCT1_resampled')
    sct = cbct_data.CBCT1_resampled;
elseif isfield(cbct_data, 'CBCT3_resampled')
    sct = cbct_data.CBCT3_resampled;
else
    error('CBCT1_resampled / CBCT3_resampled variable not found in %s', cbct_filepath);
end
if ~isfield(sct, 'cubeHU')
    error('CBCT resampled struct missing required field: cubeHU');
end

if isfield(sct, 'spacing') && ~isempty(sct.spacing)
    spacing_mm = sct.spacing(:)';
elseif ~isempty(step15_spacing_mm)
    spacing_mm = step15_spacing_mm;
    fprintf('       [INFO] Using spacing from field_dose file: [%.3f %.3f %.3f] mm\n', spacing_mm);
else
    error('CBCT resampled struct missing required field: spacing');
end

sctSize = size(sct.cubeHU);
if ~isequal(gridSize, sctSize)
    error(['Dose grid [%d %d %d] does not match CBCT grid [%d %d %d].'], ...
        Nx, Ny, Nz, sctSize(1), sctSize(2), sctSize(3));
end

fprintf('       Spacing: [%.2f, %.2f, %.2f] mm\n', spacing_mm);
fprintf('       HU range: [%.0f, %.0f]\n', min(sct.cubeHU(:)), max(sct.cubeHU(:)));

% Mask to body once (shared by both runs)
if isfield(sct, 'bodyMask')
    doseGrid = doseGrid .* double(sct.bodyMask);
end

%% ========================= GRID DOWNSCALING (shared) ====================

if CONFIG.downscale_factor ~= 1
    df     = CONFIG.downscale_factor;
    new_Nx = max(1, round(Nx / df));
    new_Ny = max(1, round(Ny / df));
    new_Nz = max(1, round(Nz / df));

    fprintf('[DS]  Downscaling: [%d x %d x %d] -> [%d x %d x %d]\n', ...
        Nx, Ny, Nz, new_Nx, new_Ny, new_Nz);

    doseGrid   = max(imresize3(doseGrid, [new_Nx, new_Ny, new_Nz]), 0);
    sct.cubeHU = imresize3(sct.cubeHU, [new_Nx, new_Ny, new_Nz]);
    if isfield(sct, 'bodyMask')
        sct.bodyMask = imresize3(single(sct.bodyMask), [new_Nx, new_Ny, new_Nz], 'nearest') > 0.5;
    end
    if isfield(sct, 'couchMask')
        sct.couchMask = imresize3(single(sct.couchMask), [new_Nx, new_Ny, new_Nz], 'nearest') > 0.5;
    end
    if isfield(sct, 'cubeDensity')
        sct.cubeDensity = imresize3(sct.cubeDensity, [new_Nx, new_Ny, new_Nz]);
    end

    spacing_mm = spacing_mm .* ([Nx, Ny, Nz] ./ [new_Nx, new_Ny, new_Nz]);
    sct.spacing = spacing_mm;
    Nx = new_Nx; Ny = new_Ny; Nz = new_Nz;
    gridSize = [Nx, Ny, Nz];
end

%% ========================= RUN BOTH AIR SPEEDS ==========================

air_speeds = CONFIG.air_sound_speeds;
n_runs     = numel(air_speeds);

for k = 1:n_runs
    air_c = air_speeds(k);
    fprintf('\n=========================================================\n');
    fprintf('  RUN %d/%d  -  air sound speed = %d m/s\n', k, n_runs, air_c);
    fprintf('=========================================================\n');

    cfg = CONFIG;
    cfg.tissue_tables = define_tissue_tables();
    % Standalone uniform-medium mapping (mirrors run_standalone_simulation).
    cfg.tissue_tables.uniform = struct( ...
        'density',     cfg.uniform_density, ...
        'sound_speed', cfg.uniform_sound_speed, ...
        'alpha_coeff', cfg.uniform_alpha_coeff, ...
        'alpha_power', 1.1, ...
        'gruneisen',   cfg.uniform_gruneisen);
    % Inject the swept air sound speed.
    cfg.tissue_tables.threshold_2.air.sound_speed = air_c;

    out = simulate_and_reconstruct(cfg, sct, doseGrid, spacing_mm, fd_gantry_angle);
    out.air_c = air_c;
    res(k) = out;

    fprintf('  [RUN %d] air=%d m/s -> minC=%.0f m/s, Nt=%d, fwd=%.1fs, tr=%.1fs\n', ...
        k, air_c, out.minC, out.Nt, out.fwd_time, out.tr_time);
end

%% ========================= SPEEDUP SUMMARY ==============================

fprintf('\n========= SPEEDUP SUMMARY =========\n');
fprintf('  %-10s  %-8s  %-8s  %-10s  %-10s\n', 'air(m/s)', 'minC', 'Nt', 'fwd(s)', 'tr(s)');
for k = 1:n_runs
    fprintf('  %-10d  %-8.0f  %-8d  %-10.1f  %-10.1f\n', ...
        res(k).air_c, res(k).minC, res(k).Nt, res(k).fwd_time, res(k).tr_time);
end
if n_runs == 2
    dt_ratio = (res(1).fwd_time + res(1).tr_time) / max(res(2).fwd_time + res(2).tr_time, eps);
    fprintf('  Total-time speedup (343 -> 1480): %.2fx\n', dt_ratio);
end
fprintf('===================================\n');

%% ========================= IDENTIFY REF / TGT ===========================
% Reference = 343 m/s (physically correct). Target = 1480 m/s (fast).

idx_ref = find(air_speeds == 343, 1);
idx_tgt = find(air_speeds == 1480, 1);
if isempty(idx_ref), idx_ref = 1; end
if isempty(idx_tgt), idx_tgt = n_runs; end

recon_ref  = res(idx_ref).recon_dose;   % 343 m/s (truth)
recon_tgt  = res(idx_tgt).recon_dose;   % 1480 m/s
spacing_mm = res(idx_ref).spacing_mm;
sensor_msk = res(idx_ref).sensor_mask;
density_bg = res(idx_ref).density;
gridSize   = res(idx_ref).gridSize;

if ~isequal(size(recon_ref), size(recon_tgt))
    error('Recon grids differ between runs: [%s] vs [%s].', ...
        num2str(size(recon_ref)), num2str(size(recon_tgt)));
end

%% ========================= NORMALIZE RECONS =============================
% Least-squares relative normalization: scale each recon by the gain that best
% fits it to its own RayStation truth (doseGrid) over the truth's 10% region, so
% the recon-vs-recon gamma is amplitude-consistent. (CONFIG.normalize)
if isfield(CONFIG, 'normalize') && CONFIG.normalize
    fprintf('\n[NORM] Least-squares normalizing each recon to the truth dose.\n');
    g_ref = least_squares_gain(res(idx_ref).doseGrid, recon_ref);
    g_tgt = least_squares_gain(res(idx_tgt).doseGrid, recon_tgt);
    recon_ref = recon_ref * g_ref;
    recon_tgt = recon_tgt * g_tgt;
    fprintf('       Gains: 343 m/s x%.4g, 1480 m/s x%.4g\n', g_ref, g_tgt);
end

%% ========================= GAMMA ANALYSIS ===============================
%  343-recon (reference) vs 1480-recon (target).

gamma_results = struct();

if exist('CalcGamma', 'file') == 2
    fprintf('\n[Gamma] Running gamma analysis (343-recon ref vs 1480-recon tgt)...\n');

    ref_struct.start = [0, 0, 0];
    ref_struct.width = spacing_mm;
    ref_struct.data  = double(recon_ref);

    tgt_struct.start = [0, 0, 0];
    tgt_struct.width = spacing_mm;
    tgt_struct.data  = double(recon_tgt);

    low_dose_cutoff = 0.10 * max(recon_ref(:));
    gamma_eval_mask = recon_ref >= low_dose_cutoff;

    gamma_criteria = {10, 10, '10%/10mm'; 5, 5, '5%/5mm'; 3, 3, '3%/3mm'};
    gamma_maps     = cell(size(gamma_criteria, 1), 1);
    pass_rates     = zeros(size(gamma_criteria, 1), 1);

    for gc = 1:size(gamma_criteria, 1)
        pct_crit = gamma_criteria{gc, 1};
        dta_crit = gamma_criteria{gc, 2};
        lbl      = gamma_criteria{gc, 3};
        fprintf('  [%s] Computing...', lbl);
        try
            gmap = CalcGamma(ref_struct, tgt_struct, pct_crit, dta_crit, ...
                'local', 0, 'limit', dta_crit * 2, 'restrict', 1);
            gamma_maps{gc} = gmap;
            eval_vals      = gmap(gamma_eval_mask);
            pass_rate      = 100 * mean(eval_vals <= 1);
            pass_rates(gc) = pass_rate;
            fprintf('  Pass rate: %.2f%%\n', pass_rate);
        catch ME
            warning('Gamma [%s] failed: %s', lbl, ME.message);
            gamma_maps{gc} = [];
            pass_rates(gc) = NaN;
            fprintf('  FAILED\n');
        end
    end

    gamma_results.maps        = gamma_maps;
    gamma_results.pass_rates  = pass_rates;
    gamma_results.criteria    = gamma_criteria;
    gamma_results.cutoff_Gy   = low_dose_cutoff;
    gamma_results.eval_mask   = gamma_eval_mask;

    fprintf('\n  ------ Gamma Pass Rates (343-recon vs 1480-recon, 10%% cutoff) ------\n');
    for gc = 1:size(gamma_criteria, 1)
        if isnan(pass_rates(gc))
            fprintf('  %-12s  FAILED\n', gamma_criteria{gc, 3});
        else
            fprintf('  %-12s  %.2f%%\n', gamma_criteria{gc, 3}, pass_rates(gc));
        end
    end
else
    warning('CalcGamma not found. Skipping gamma analysis.');
    gamma_results = [];
end

%% ========================= SAVE RESULTS =================================

if CONFIG.save_results
    results = struct();
    results.recon_343      = recon_ref;
    results.recon_1480     = recon_tgt;
    results.air_speeds     = air_speeds;
    results.minC           = [res.minC];
    results.Nt             = [res.Nt];
    results.fwd_time_sec   = [res.fwd_time];
    results.tr_time_sec    = [res.tr_time];
    results.spacing_mm     = spacing_mm;
    results.grid_size      = gridSize;
    results.sensor_mask    = sensor_msk;
    results.config         = CONFIG;
    if ~isempty(gamma_results)
        results.gamma = gamma_results;
    end
    save(CONFIG.output_file, '-struct', 'results', '-v7.3');
    fprintf('\nResults saved to: %s\n', CONFIG.output_file);
end

%% ========================= VISUALIZATION ================================

if CONFIG.plot_results
    % Shared sensor/dose geometry (identical across runs).
    dose_mask_vis = double(recon_ref) >= 0.10 * max(double(recon_ref(:)));
    plot_sensor_dose_planes(dose_mask_vis, logical(sensor_msk), spacing_mm, density_bg, CONFIG);

    % Recon comparison: 343 (top) vs 1480 (bottom).
    plot_dose_panels(recon_ref, recon_tgt, sensor_msk, density_bg, spacing_mm, ...
        'Recon Dose: air 343 m/s (top) vs 1480 m/s (bottom)', CONFIG.viz_smooth_sigma);

    % Gamma + absolute difference at the reference max-dose axial slice.
    if ~isempty(gamma_results) && isfield(gamma_results, 'maps')
        [~, max_idx]     = max(recon_ref(:));
        [~, ~, cz_gamma] = ind2sub(gridSize, max_idx);
        plot_gamma_and_error_axial(gamma_results, recon_ref, recon_tgt, sensor_msk, cz_gamma);
    end
end

fprintf('\nAir sound-speed comparison complete.\n');


%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

function g = least_squares_gain(rs_truth, recon)
%LEAST_SQUARES_GAIN Scalar gain aligning a recon to its RS truth (relative norm).
%  g = sum(rs.*recon)/sum(recon.^2) over the truth's 10% low-dose region; the g
%  minimizing ||rs - g*recon||^2 there. Falls back to 1 for empty/zero inputs.
%  (Method mirrored from study_pass_rates_allsegments.m.)
    rs_truth = double(rs_truth);
    recon    = double(recon);
    if max(rs_truth(:)) > 0
        mask = rs_truth >= 0.10 * max(rs_truth(:));
    else
        mask = true(size(rs_truth));
    end
    r = recon(mask);
    denom = sum(r .^ 2);
    if denom > 0
        g = sum(rs_truth(mask) .* r) / denom;
    else
        g = 1;
    end
end


function out = simulate_and_reconstruct(CONFIG, sct, doseGrid, spacing_mm, fd_gantry_angle)
%SIMULATE_AND_RECONSTRUCT  Full forward + reconstruction pipeline for one air
%  sound speed. Ported faithfully from run_standalone_simulation.m (medium ->
%  p0 -> padding -> sensor -> forward -> recon -> pressure->dose). Noise is
%  disabled (CONFIG.convolution_kernel = 0); only the sensor frequency response
%  is applied. Returns recon_dose plus geometry and timing/step diagnostics.

    dx = spacing_mm(1) / 1000;
    dy = spacing_mm(2) / 1000;
    dz = spacing_mm(3) / 1000;

    gridSize = size(doseGrid);
    Nx = gridSize(1); Ny = gridSize(2); Nz = gridSize(3);

    %% ---------------- CREATE ACOUSTIC MEDIUM ----------------
    fprintf('[3/7] Creating acoustic medium (method: %s)...\n', CONFIG.gruneisen_method);
    medium = create_acoustic_medium(sct, CONFIG);

    gs = medium.grid_size;
    if CONFIG.force_uniform_density
        medium.density = ones(gs) * CONFIG.uniform_density;
    end
    if CONFIG.force_uniform_sound_speed
        medium.sound_speed = ones(gs) * CONFIG.uniform_sound_speed;
    end
    if CONFIG.force_uniform_attenuation
        medium.alpha_coeff = ones(gs) * CONFIG.uniform_alpha_coeff;
        medium.alpha_power = 1.1;
    end
    if CONFIG.force_uniform_gruneisen
        medium.gruneisen = ones(gs) * CONFIG.uniform_gruneisen;
    end
    if isfield(sct, 'bodyMask')
        outside_body = ~logical(sct.bodyMask);
        if isfield(sct, 'couchMask')
            outside_body = outside_body | logical(sct.couchMask);
        end
        medium.density(outside_body)     = CONFIG.uniform_density;
        medium.sound_speed(outside_body) = CONFIG.uniform_sound_speed;
        medium.alpha_coeff(outside_body) = CONFIG.uniform_alpha_coeff;
        medium.gruneisen(outside_body)   = CONFIG.uniform_gruneisen;
    end

    fprintf('       Density range:     [%.0f, %.0f] kg/m^3\n', min(medium.density(:)), max(medium.density(:)));
    fprintf('       Sound speed range: [%.0f, %.0f] m/s\n',    min(medium.sound_speed(:)), max(medium.sound_speed(:)));
    fprintf('       Gruneisen range:   [%.4f, %.4f]\n',        min(medium.gruneisen(:)), max(medium.gruneisen(:)));

    %% ---------------- INITIAL PRESSURE p0 ----------------
    if isfield(sct, 'bodyMask')
        doseGrid = doseGrid .* double(sct.bodyMask);
    end

    fprintf('[4/7] Computing initial pressure...\n');
    meterset       = CONFIG.meterset;
    num_pulses     = ceil(meterset / CONFIG.dose_per_pulse_cGy);
    dose_per_pulse = doseGrid / num_pulses;

    p0 = dose_per_pulse .* medium.gruneisen .* medium.density;
    p0 = smooth(p0);

    fprintf('       Meterset: %.2f MU -> %d pulses\n', meterset, num_pulses);
    fprintf('       p0 range: [%.2e, %.2e] Pa\n', min(p0(:)), max(p0(:)));

    doseThreshold = 0.01 * max(doseGrid(:));
    doseMask      = doseGrid > doseThreshold;

    if ~any(doseMask(:)) || max(p0(:)) == 0
        error('simulate_and_reconstruct:NoSignal', ...
            'No significant dose or zero initial pressure.');
    end

    %% ---------------- OPTIMAL GRID PADDING ----------------
    fprintf('[PAD] Computing FFT-optimal padded dimensions...\n');
    Nx_orig = Nx; Ny_orig = Ny; Nz_orig = Nz;
    gridSize_orig = gridSize;
    medium_orig   = medium;

    if CONFIG.use_grid_padding
        Nx_pad = find_optimal_kwave_size(Nx, CONFIG.pml_size);
        Ny_pad = find_optimal_kwave_size(Ny, CONFIG.pml_size);
        Nz_pad = find_optimal_kwave_size(Nz, CONFIG.pml_size);
    else
        Nx_pad = Nx; Ny_pad = Ny; Nz_pad = Nz;
    end

    did_pad = ~isequal([Nx_pad, Ny_pad, Nz_pad], [Nx, Ny, Nz]);
    if did_pad
        fprintf('[PAD] Padding grid: [%d %d %d] -> [%d %d %d]\n', Nx, Ny, Nz, Nx_pad, Ny_pad, Nz_pad);
        density_pad    = ones(Nx_pad, Ny_pad, Nz_pad) * 1000;
        soundSpeed_pad = ones(Nx_pad, Ny_pad, Nz_pad) * 1540;
        alphaCoeff_pad = zeros(Nx_pad, Ny_pad, Nz_pad);
        gruneisen_pad  = zeros(Nx_pad, Ny_pad, Nz_pad);

        density_pad(1:Nx, 1:Ny, 1:Nz)    = medium.density;
        soundSpeed_pad(1:Nx, 1:Ny, 1:Nz) = medium.sound_speed;
        if numel(medium.alpha_coeff) > 1
            alphaCoeff_pad(1:Nx, 1:Ny, 1:Nz) = medium.alpha_coeff;
        else
            alphaCoeff_pad(:) = medium.alpha_coeff;
        end
        gruneisen_pad(1:Nx, 1:Ny, 1:Nz) = medium.gruneisen;

        medium.density     = density_pad;
        medium.sound_speed = soundSpeed_pad;
        medium.alpha_coeff = alphaCoeff_pad;
        medium.gruneisen   = gruneisen_pad;

        p0_pad = zeros(Nx_pad, Ny_pad, Nz_pad);
        p0_pad(1:Nx, 1:Ny, 1:Nz) = p0;
        p0 = p0_pad;

        Nx = Nx_pad; Ny = Ny_pad; Nz = Nz_pad;
        gridSize = [Nx, Ny, Nz];
    else
        fprintf('[PAD] Grid [%d %d %d] already FFT-optimal.\n', Nx, Ny, Nz);
    end

    %% ---------------- SENSOR PLACEMENT ----------------
    fprintf('[5/7] Placing sensor (method: %s)...\n', CONFIG.sensor_placement_method);
    sensor      = struct();
    sensor.mask = zeros(Nx, Ny, Nz);
    sensor_info_orig = struct('num_elements', 0);

    switch CONFIG.sensor_placement_method
        case 'full_plane_anterior'
            sensor.mask(CONFIG.sensor_x_index, :, :) = 1;
        case 'full_plane_lateral'
            sensor.mask(:, CONFIG.sensor_y_index, :) = 1;
        case 'spherical'
            sph_radius  = floor(min([Nx, Ny, Nz]) / 2) - CONFIG.pml_size;
            sensor.mask = makeSphere(Nx, Ny, Nz, sph_radius);
            sph_cx = floor(Nx/2) + 1; sph_cy = floor(Ny/2) + 1; sph_cz = floor(Nz/2) + 1;
            [Xg_sph, Yg_sph, Zg_sph] = ndgrid(1:Nx, 1:Ny, 1:Nz);
            ball_mask = (Xg_sph - sph_cx).^2 + (Yg_sph - sph_cy).^2 + (Zg_sph - sph_cz).^2 <= sph_radius^2;
            p0 = p0 .* ball_mask;
            clear Xg_sph Yg_sph Zg_sph
        case 'box'
            bx_lo = 3; bx_hi_x = Nx - 3; bx_hi_y = Ny - 3; bx_hi_z = Nz - 3;
            if bx_hi_x <= bx_lo || bx_hi_y <= bx_lo || bx_hi_z <= bx_lo
                error('simulate_and_reconstruct:BoxTooSmall', ...
                    'Grid [%d %d %d] too small for box sensor.', Nx, Ny, Nz);
            end
            sensor.mask(bx_lo,   bx_lo:bx_hi_y, bx_lo:bx_hi_z) = 1;
            sensor.mask(bx_hi_x, bx_lo:bx_hi_y, bx_lo:bx_hi_z) = 1;
            sensor.mask(bx_lo:bx_hi_x, bx_lo,   bx_lo:bx_hi_z) = 1;
            sensor.mask(bx_lo:bx_hi_x, bx_hi_y, bx_lo:bx_hi_z) = 1;
            sensor.mask(bx_lo:bx_hi_x, bx_lo:bx_hi_y, bx_lo)   = 1;
            sensor.mask(bx_lo:bx_hi_x, bx_lo:bx_hi_y, bx_hi_z) = 1;
        case 'determine_sensor_mask'
            sct_for_sensor = sct;
            if ~isfield(sct_for_sensor, 'couchMask')
                sct_for_sensor.couchMask = false(size(sct_for_sensor.bodyMask));
            end
            if ~isfield(sct_for_sensor, 'origin')
                sct_for_sensor.origin = [0, 0, 0];
            end
            sct_for_sensor.spacing = spacing_mm;

            field_dose_for_sensor = struct();
            field_dose_for_sensor.dose_Gy     = doseGrid;
            if ~isempty(fd_gantry_angle)
                field_dose_for_sensor.gantry_angle = fd_gantry_angle;
            else
                field_dose_for_sensor.gantry_angle = 0;
            end
            field_dose_for_sensor.origin      = sct_for_sensor.origin;
            field_dose_for_sensor.spacing     = spacing_mm;
            field_dose_for_sensor.dimensions  = [Nx_orig, Ny_orig, Nz_orig];

            beam_meta = [];
            if isfield(CONFIG, 'beam_metadata') && ~isempty(CONFIG.beam_metadata)
                beam_meta = CONFIG.beam_metadata;
            end

            [sensor_mask_orig, sensor_info_orig] = determine_sensor_mask( ...
                sct_for_sensor, field_dose_for_sensor, beam_meta, CONFIG);

            gp = sensor_info_orig.grid_pad;
            if gp.expanded
                fprintf('       [Sensor] Grid expansion: Y(+%d/+%d), X(+%d/+%d), Z(+%d/+%d). Re-padding with water.\n', ...
                    gp.y_pre, gp.y_post, gp.x_pre, gp.x_post, gp.z_pre, gp.z_post);

                density_unp    = medium.density(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
                soundSpeed_unp = medium.sound_speed(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
                if numel(medium.alpha_coeff) > 1
                    alphaCoeff_unp = medium.alpha_coeff(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
                else
                    alphaCoeff_unp = medium.alpha_coeff;
                end
                gruneisen_unp  = medium.gruneisen(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
                p0_unp         = p0(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);

                Nx_exp = Nx_orig + gp.y_pre + gp.y_post;
                Ny_exp = Ny_orig + gp.x_pre + gp.x_post;
                Nz_exp = Nz_orig + gp.z_pre + gp.z_post;

                density_exp    = ones(Nx_exp, Ny_exp, Nz_exp)  * 1000;
                soundSpeed_exp = ones(Nx_exp, Ny_exp, Nz_exp)  * 1540;
                alphaCoeff_exp = zeros(Nx_exp, Ny_exp, Nz_exp);
                gruneisen_exp  = zeros(Nx_exp, Ny_exp, Nz_exp);
                p0_exp         = zeros(Nx_exp, Ny_exp, Nz_exp);

                xr = gp.y_pre + (1:Nx_orig);
                yr = gp.x_pre + (1:Ny_orig);
                zr = gp.z_pre + (1:Nz_orig);

                density_exp(xr, yr, zr)    = density_unp;
                soundSpeed_exp(xr, yr, zr) = soundSpeed_unp;
                if numel(alphaCoeff_unp) > 1
                    alphaCoeff_exp(xr, yr, zr) = alphaCoeff_unp;
                else
                    alphaCoeff_exp(:) = alphaCoeff_unp;
                end
                gruneisen_exp(xr, yr, zr)  = gruneisen_unp;
                p0_exp(xr, yr, zr)         = p0_unp;

                medium.density     = density_exp;
                medium.sound_speed = soundSpeed_exp;
                medium.alpha_coeff = alphaCoeff_exp;
                medium.gruneisen   = gruneisen_exp;
                p0 = p0_exp;

                medium_orig = struct( ...
                    'density',     density_exp, ...
                    'sound_speed', soundSpeed_exp, ...
                    'alpha_coeff', alphaCoeff_exp, ...
                    'gruneisen',   gruneisen_exp);

                doseGrid_exp = zeros(Nx_exp, Ny_exp, Nz_exp);
                doseGrid_exp(xr, yr, zr) = doseGrid;
                doseGrid = doseGrid_exp;

                doseMask_exp = false(Nx_exp, Ny_exp, Nz_exp);
                doseMask_exp(xr, yr, zr) = doseMask;
                doseMask = doseMask_exp;

                if isfield(sct, 'bodyMask') && isequal(size(sct.bodyMask), [Nx_orig, Ny_orig, Nz_orig])
                    body_exp = false(Nx_exp, Ny_exp, Nz_exp);
                    body_exp(xr, yr, zr) = sct.bodyMask;
                    sct.bodyMask = body_exp;
                end
                if isfield(sct, 'couchMask') && ~isempty(sct.couchMask) && ...
                        isequal(size(sct.couchMask), [Nx_orig, Ny_orig, Nz_orig])
                    couch_exp = false(Nx_exp, Ny_exp, Nz_exp);
                    couch_exp(xr, yr, zr) = sct.couchMask;
                    sct.couchMask = couch_exp;
                end

                Nx_orig = Nx_exp; Ny_orig = Ny_exp; Nz_orig = Nz_exp;
                gridSize_orig = [Nx_orig, Ny_orig, Nz_orig];

                if CONFIG.use_grid_padding
                    Nx_pad2 = find_optimal_kwave_size(Nx_orig, CONFIG.pml_size);
                    Ny_pad2 = find_optimal_kwave_size(Ny_orig, CONFIG.pml_size);
                    Nz_pad2 = find_optimal_kwave_size(Nz_orig, CONFIG.pml_size);
                else
                    Nx_pad2 = Nx_orig; Ny_pad2 = Ny_orig; Nz_pad2 = Nz_orig;
                end

                if ~isequal([Nx_pad2, Ny_pad2, Nz_pad2], [Nx_orig, Ny_orig, Nz_orig])
                    density_pad    = ones(Nx_pad2, Ny_pad2, Nz_pad2)  * 1000;
                    soundSpeed_pad = ones(Nx_pad2, Ny_pad2, Nz_pad2)  * 1540;
                    alphaCoeff_pad = zeros(Nx_pad2, Ny_pad2, Nz_pad2);
                    gruneisen_pad  = zeros(Nx_pad2, Ny_pad2, Nz_pad2);
                    p0_pad2        = zeros(Nx_pad2, Ny_pad2, Nz_pad2);

                    density_pad(1:Nx_orig, 1:Ny_orig, 1:Nz_orig)    = medium.density;
                    soundSpeed_pad(1:Nx_orig, 1:Ny_orig, 1:Nz_orig) = medium.sound_speed;
                    if numel(medium.alpha_coeff) > 1
                        alphaCoeff_pad(1:Nx_orig, 1:Ny_orig, 1:Nz_orig) = medium.alpha_coeff;
                    else
                        alphaCoeff_pad(:) = medium.alpha_coeff;
                    end
                    gruneisen_pad(1:Nx_orig, 1:Ny_orig, 1:Nz_orig) = medium.gruneisen;
                    p0_pad2(1:Nx_orig, 1:Ny_orig, 1:Nz_orig)       = p0;

                    medium.density     = density_pad;
                    medium.sound_speed = soundSpeed_pad;
                    medium.alpha_coeff = alphaCoeff_pad;
                    medium.gruneisen   = gruneisen_pad;
                    p0 = p0_pad2;
                end

                Nx = Nx_pad2; Ny = Ny_pad2; Nz = Nz_pad2;
                gridSize = [Nx, Ny, Nz];
                sensor.mask = zeros(Nx, Ny, Nz);
            end

            if isfield(CONFIG, 'plot_exclusion_zone') && CONFIG.plot_exclusion_zone
                plot_exclusion_zone_views(sct, sensor_info_orig, spacing_mm, ...
                    sprintf('Exclusion zone (gantry %.1f deg)', field_dose_for_sensor.gantry_angle));
            end

            m1 = min(Nx, size(sensor_mask_orig, 1));
            m2 = min(Ny, size(sensor_mask_orig, 2));
            m3 = min(Nz, size(sensor_mask_orig, 3));
            sensor.mask(1:m1, 1:m2, 1:m3) = double(sensor_mask_orig(1:m1, 1:m2, 1:m3));
        case 'fixed_anterior'
            fixed_struct = struct();
            if isfield(CONFIG, 'beam_metadata') && ~isempty(CONFIG.beam_metadata)
                fixed_struct.beam_metadata = CONFIG.beam_metadata;
            end
            fixed_struct.total_dose = doseGrid;
            [sensor_mask_orig, ~] = determine_sensor_placement_fixed(CONFIG, sct, fixed_struct);
            m1 = min(Nx, size(sensor_mask_orig, 1));
            m2 = min(Ny, size(sensor_mask_orig, 2));
            m3 = min(Nz, size(sensor_mask_orig, 3));
            sensor.mask(1:m1, 1:m2, 1:m3) = double(sensor_mask_orig(1:m1, 1:m2, 1:m3));
        otherwise
            error('Unknown sensor_placement_method: "%s"', CONFIG.sensor_placement_method);
    end

    numSensorPts = sum(sensor.mask(:));
    fprintf('       Sensor: %d active points\n', numSensorPts);
    if numSensorPts == 0
        error('simulate_and_reconstruct:EmptySensor', 'Sensor mask is empty.');
    end

    %% ---------------- k-WAVE GRID & MEDIUM SETUP ----------------
    fprintf('[6/7] Setting up k-Wave grid...\n');
    kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

    maxC = max(medium.sound_speed(:));
    minC = min(medium.sound_speed(medium.sound_speed > 0));
    dt   = CONFIG.cfl_number * min([dx, dy, dz]) / maxC;

    gridDiag = sqrt((Nx*dx)^2 + (Ny*dy)^2 + (Nz*dz)^2);
    simTime  = 2.5 * gridDiag / minC;
    Nt       = ceil(simTime / dt);

    kgrid.dt = dt;
    kgrid.Nt = Nt;
    fprintf('       dt = %.2e s, Nt = %d, minC = %.0f m/s, T_sim = %.2e s\n', dt, Nt, minC, simTime);

    kmedium             = struct();
    kmedium.density     = medium.density;
    kmedium.sound_speed = medium.sound_speed;
    kmedium.alpha_coeff = medium.alpha_coeff;
    kmedium.alpha_power = 1.1;

    if CONFIG.use_gpu
        try
            gpuDevice;
            dataCast = 'gpuArray-single';
        catch
            dataCast = 'single';
        end
    else
        dataCast = 'single';
    end

    inputArgs = {'Smooth', false, 'PMLInside', false, 'PMLSize', CONFIG.pml_size, ...
                 'DataCast', dataCast, 'PlotSim', false};

    %% ---------------- FORWARD SIMULATION ----------------
    fprintf('\n[7/7] Running k-Wave forward simulation...\n');
    source_fwd    = struct();
    source_fwd.p0 = p0;

    fwd_tic    = tic;
    sensorData = kspaceFirstOrder3D(kgrid, kmedium, source_fwd, sensor, inputArgs{:});
    fwd_time   = toc(fwd_tic);
    fprintf('       Forward complete (%.1f s). Sensor data: [%d x %d]\n', ...
        fwd_time, size(sensorData, 1), size(sensorData, 2));

    sensorData = smooth(sensorData);
    FS         = 1 / kgrid.dt;
    sensorData_measured = sensorData;

    %% ---------------- FREQUENCY RESPONSE (noise disabled) ----------------
    if CONFIG.convolution_kernel > 0
        conv_kernel_sigma  = CONFIG.convolution_kernel;
        conv_noise_level   = CONFIG.conv_noise_level;
        conv_deconv_lambda = CONFIG.conv_deconv_lambda;

        sigma_samples = conv_kernel_sigma / dt;
        kernel_half   = ceil(4 * sigma_samples);
        t_kernel      = (-kernel_half : kernel_half)';
        gauss_kernel  = exp(-t_kernel.^2 / (2 * sigma_samples^2));
        gauss_kernel  = gauss_kernel / sum(gauss_kernel);

        sensorData_cpu = double(gather(sensorData));
        Nt_data        = size(sensorData_cpu, 2);
        H       = fft(gauss_kernel, Nt_data).';
        H_conj  = conj(H);
        H_power = abs(H).^2;

        sensorData_conv  = real(ifft(fft(sensorData_cpu, [], 2) .* H, [], 2));
        sensorData_resp  = gaussianFilter(sensorData_conv, FS, 0.35e6, 100, true);
        noise_amp        = conv_noise_level * max(abs(sensorData_resp(:)));
        sensorData_noisy = sensorData_resp + noise_amp * randn(size(sensorData_resp));
        sensorData_deconv = real(ifft( ...
            fft(sensorData_noisy, [], 2) .* H_conj ./ (H_power + conv_deconv_lambda), [], 2));

        sensorData          = single(sensorData_deconv);
        sensorData_measured = single(sensorData_deconv);
    else
        sensorData          = gaussianFilter(sensorData, FS, 0.35e6, 100, true);
        sensorData_measured = sensorData;
        fprintf('       Pulse convolution disabled; frequency response applied.\n');
    end

    %% ---------------- RECONSTRUCTION ----------------
    das_opts = struct();
    if isfield(CONFIG, 'das_use_elements'), das_opts.use_elements = CONFIG.das_use_elements; end
    if isfield(CONFIG, 'das_envelope'),     das_opts.envelope     = CONFIG.das_envelope;     end
    if isfield(CONFIG, 'das_use_aperture'), das_opts.use_aperture = CONFIG.das_use_aperture; end
    if isfield(CONFIG, 'das_aperture_cos'), das_opts.aperture_cos = CONFIG.das_aperture_cos; end
    if isfield(CONFIG, 'das_depth_weight'), das_opts.depth_weight = CONFIG.das_depth_weight; end
    if isfield(CONFIG, 'das_interp'),       das_opts.interp       = CONFIG.das_interp;       end

    if ~exist('sensor_info_orig', 'var') || ~isstruct(sensor_info_orig)
        sensor_info_orig = struct('num_elements', 0);
    end

    tr_time = 0;
    switch lower(CONFIG.reconstruction_method)
    case 'tr'
        fprintf('       Running iterative time reversal (%d iterations, tol=%.1e)...\n', ...
            CONFIG.num_time_reversal_iter, CONFIG.convergence_tol);

        reconPressure      = zeros(gridSize);
        reconPressure_prev = zeros(gridSize);
        num_iters_done     = 0;

        tr_total_tic = tic;
        for tr_iter = 1:CONFIG.num_time_reversal_iter
            fprintf('       --- TR Iteration %d/%d ---\n', tr_iter, CONFIG.num_time_reversal_iter);

            source_tr        = struct();
            source_tr.p_mask = sensor.mask;
            source_tr.p      = fliplr(sensorData);
            source_tr.p_mode = 'dirichlet';

            sensor_tr        = struct();
            sensor_tr.mask   = ones(Nx, Ny, Nz);
            sensor_tr.record = {'p_final'};

            p0_recon = kspaceFirstOrder3D(kgrid, kmedium, source_tr, sensor_tr, inputArgs{:});
            if isstruct(p0_recon) && isfield(p0_recon, 'p_final')
                reconPressure = reshape(p0_recon.p_final, [Nx, Ny, Nz]);
            else
                reconPressure = reshape(p0_recon, [Nx, Ny, Nz]);
            end
            reconPressure = max(reconPressure, 0);
            num_iters_done = tr_iter;
            fprintf('       Max pressure: %.4e Pa\n', max(reconPressure(:)));

            converged = false;
            if tr_iter > 1
                norm_prev = norm(reconPressure_prev(:));
                if norm_prev > 0
                    rel_change = norm(reconPressure(:) - reconPressure_prev(:)) / norm_prev;
                else
                    rel_change = Inf;
                end
                fprintf('       Rel change: %.4e\n', rel_change);
                if rel_change < CONFIG.convergence_tol
                    fprintf('       *** Converged at iteration %d ***\n', tr_iter);
                    converged = true;
                end
            end
            reconPressure_prev = reconPressure;

            if converged, break; end

            if tr_iter < CONFIG.num_time_reversal_iter
                source_resid    = struct();
                source_resid.p0 = reconPressure;
                sensorDataRecon = kspaceFirstOrder3D(kgrid, kmedium, source_resid, sensor, inputArgs{:});
                sensorData      = sensorData + (sensorData_measured - sensorDataRecon);
            end
        end
        reconPressure = gather(reconPressure) * CONFIG.correction_factor;
        tr_time = toc(tr_total_tic);
        fprintf('       Time reversal complete (%.1f s).\n', tr_time);

    case 'das'
        fprintf('       Running Delay-And-Sum reconstruction...\n');
        das_tic = tic;
        reconPressure = das_reconstruct(sensorData, sensor, sensor_info_orig, medium, ...
                                        Nx, Ny, Nz, dx, dy, dz, dt, das_opts);
        reconPressure = reconPressure * CONFIG.correction_factor;
        tr_time = toc(das_tic);
        fprintf('       DAS complete (%.1f s).\n', tr_time);

    case 'hybrid'
        fprintf('       Running HYBRID reconstruction (DAS seed + TR)...\n');
        hybrid_tic = tic;
        N_iter = CONFIG.num_time_reversal_iter;
        reconPressure = das_reconstruct(sensorData, sensor, sensor_info_orig, medium, ...
                                        Nx, Ny, Nz, dx, dy, dz, dt, das_opts);
        reconPressure_prev = reconPressure;
        if N_iter > 1
            source_resid    = struct();
            source_resid.p0 = reconPressure;
            sensorDataRecon = kspaceFirstOrder3D(kgrid, kmedium, source_resid, sensor, inputArgs{:});
            sensorData      = sensorData + (sensorData_measured - sensorDataRecon);
        end
        for tr_iter = 2:N_iter
            fprintf('       --- TR Iteration %d/%d ---\n', tr_iter, N_iter);
            source_tr        = struct();
            source_tr.p_mask = sensor.mask;
            source_tr.p      = fliplr(sensorData);
            source_tr.p_mode = 'dirichlet';
            sensor_tr        = struct();
            sensor_tr.mask   = ones(Nx, Ny, Nz);
            sensor_tr.record = {'p_final'};
            p0_recon = kspaceFirstOrder3D(kgrid, kmedium, source_tr, sensor_tr, inputArgs{:});
            if isstruct(p0_recon) && isfield(p0_recon, 'p_final')
                reconPressure = reshape(p0_recon.p_final, [Nx, Ny, Nz]);
            else
                reconPressure = reshape(p0_recon, [Nx, Ny, Nz]);
            end
            reconPressure = max(reconPressure, 0);
            norm_prev = norm(reconPressure_prev(:));
            if norm_prev > 0
                rel_change = norm(reconPressure(:) - reconPressure_prev(:)) / norm_prev;
            else
                rel_change = Inf;
            end
            reconPressure_prev = reconPressure;
            if rel_change < CONFIG.convergence_tol
                fprintf('       *** Converged at iteration %d ***\n', tr_iter);
                break;
            end
            if tr_iter < N_iter
                source_resid    = struct();
                source_resid.p0 = reconPressure;
                sensorDataRecon = kspaceFirstOrder3D(kgrid, kmedium, source_resid, sensor, inputArgs{:});
                sensorData      = sensorData + (sensorData_measured - sensorDataRecon);
            end
        end
        reconPressure = gather(reconPressure) * CONFIG.correction_factor;
        tr_time = toc(hybrid_tic);
        fprintf('       Hybrid complete (%.1f s).\n', tr_time);
    otherwise
        error('Unknown reconstruction_method: "%s"', CONFIG.reconstruction_method);
    end

    %% ---------------- CROP TO ORIGINAL SIZE ----------------
    if did_pad
        reconPressure = reconPressure(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
        Nx = Nx_orig; Ny = Ny_orig; Nz = Nz_orig;
        gridSize    = gridSize_orig;
        medium      = medium_orig;
        sensor.mask = sensor.mask(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
    end

    %% ---------------- PRESSURE SCALE CORRECTION ----------------
    if CONFIG.use_pressure_scale_correction
        p0_max_orig = max(p0(1:Nx_orig, 1:Ny_orig, 1:Nz_orig), [], 'all');
        recon_max   = max(reconPressure(:));
        if recon_max > 0
            reconPressure = reconPressure * (p0_max_orig / recon_max);
        end
    end

    %% ---------------- PRESSURE -> DOSE ----------------
    conversionFactor = medium.gruneisen .* medium.density;
    conversionFactor(conversionFactor == 0) = 1;
    reconDosePerPulse = reconPressure ./ conversionFactor;

    body_mask_plot = ones(gridSize);
    if isfield(sct, 'bodyMask') && isequal(size(sct.bodyMask), gridSize)
        body_mask_plot = double(sct.bodyMask);
    end
    if isfield(CONFIG, 'mask_recon_to_dose_region') && ~CONFIG.mask_recon_to_dose_region
        recon_dose = reconDosePerPulse * num_pulses .* body_mask_plot;
    else
        recon_dose = reconDosePerPulse * num_pulses .* double(doseMask) .* body_mask_plot;
    end
    % Zero air-classified voxels (no real PA signal).
    recon_dose(medium.density < 100) = 0;

    if did_pad
        p0 = p0(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
    end

    fprintf('       Reconstructed dose: [%.4f, %.4f] Gy\n', min(recon_dose(:)), max(recon_dose(:)));

    %% ---------------- BUILD OUTPUT ----------------
    out = struct();
    out.recon_dose  = recon_dose;
    out.doseGrid    = doseGrid;
    out.spacing_mm  = spacing_mm;
    out.sensor_mask = sensor.mask;
    out.density     = medium_orig.density;
    out.gridSize    = gridSize;
    out.minC        = minC;
    out.Nt          = Nt;
    out.fwd_time    = fwd_time;
    out.tr_time     = tr_time;
end


function plot_sensor_dose_planes(dose_mask, sensor_mask, spacing_mm, density, config)
%PLOT_SENSOR_DOSE_PLANES  1x3 anatomical view of sensor geometry vs dose mask.

    [Nx3, Ny3, Nz3] = size(dose_mask);
    x_ax = (1:Nx3) * spacing_mm(1);
    y_ax = (1:Ny3) * spacing_mm(2);
    z_ax = (1:Nz3) * spacing_mm(3);

    dose_cor  = squeeze(any(dose_mask,   2));
    dose_sag  = squeeze(any(dose_mask,   1));
    dose_axi  = squeeze(any(dose_mask,   3));
    sens_cor  = squeeze(any(sensor_mask, 2));
    sens_sag  = squeeze(any(sensor_mask, 1));
    sens_axi  = squeeze(any(sensor_mask, 3));

    have_density = ~isempty(density) && isequal(size(density), size(dose_mask));
    wl_center = 1050; wl_width = 350;
    wl_min    = wl_center - wl_width / 2;
    if have_density
        ct_projs = {
            squeeze(mean(double(density), 2))', ...
            squeeze(mean(double(density), 1))', ...
            squeeze(mean(double(density), 3))'
        };
    end

    figure('Name', 'Sensor Placement vs Dose Mask', 'Color', 'w', ...
        'NumberTitle', 'off', 'Position', [80, 80, 1300, 420]);
    sgtitle(sprintf('Sensor Placement vs Recon Dose Mask (\\geq10%% max)   |   Sensor: %s', ...
        config.sensor_placement_method), 'FontWeight', 'bold', 'FontSize', 11);

    view_data = {
        dose_cor', sens_cor', x_ax, z_ax, 'X (mm)', 'Z (mm)', 'Coronal  (mean-proj along Y)';
        dose_sag', sens_sag', y_ax, z_ax, 'Y (mm)', 'Z (mm)', 'Sagittal  (mean-proj along X)';
        dose_axi', sens_axi', x_ax, y_ax, 'X (mm)', 'Y (mm)', 'Axial  (mean-proj along Z)';
    };

    dose_color   = [0.20, 0.50, 0.90];
    sensor_color = [0.90, 0.10, 0.10];

    for col = 1:3
        ax  = subplot(1, 3, col);
        d2d = double(view_data{col, 1});
        s2d = double(view_data{col, 2});
        xv  = view_data{col, 3};
        yv  = view_data{col, 4};

        if have_density
            dn     = (ct_projs{col} - wl_min) / wl_width;
            dn     = max(0, min(1, dn));
            bg_rgb = repmat(dn, [1, 1, 3]);
        else
            bg_rgb = ones([size(d2d), 3]);
        end
        image(ax, xv, yv, bg_rgb);
        hold(ax, 'on');

        if any(d2d(:))
            dose_rgb = repmat(reshape(dose_color, [1,1,3]), [size(d2d), 1]);
            h_dose   = image(ax, xv, yv, dose_rgb);
            h_dose.AlphaData = d2d * 0.45;
        end
        if any(s2d(:))
            sens_rgb = repmat(reshape(sensor_color, [1,1,3]), [size(s2d), 1]);
            h_sens   = image(ax, xv, yv, sens_rgb);
            h_sens.AlphaData = s2d * 0.85;
        end

        hold(ax, 'off');
        axis(ax, 'xy'); axis(ax, 'image');
        xlabel(ax, view_data{col, 5}); ylabel(ax, view_data{col, 6});
        title(ax, view_data{col, 7}, 'FontSize', 10);

        hold(ax, 'on');
        p1 = patch(ax, NaN, NaN, dose_color,   'FaceAlpha', 0.45, 'EdgeColor', 'none');
        p2 = patch(ax, NaN, NaN, sensor_color, 'FaceAlpha', 0.85, 'EdgeColor', 'none');
        hold(ax, 'off');
        if col == 3
            legend(ax, [p1, p2], {'Recon dose (>=10%)', 'Sensor'}, ...
                'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 8);
        end
    end
    drawnow;
end


function plot_dose_panels(original, recon, sensor_mask, density, spacing_mm, titleStr, smooth_sigma)
%PLOT_DOSE_PANELS 2x3 dose comparison: coronal, sagittal, axial.
%  Row 1 = first volume (343 m/s), Row 2 = second volume (1480 m/s).

    if nargin < 7 || isempty(smooth_sigma), smooth_sigma = 0; end

    gridSize = size(original);
    if ~isequal(size(sensor_mask), gridSize)
        sensor_mask = sensor_mask(1:gridSize(1), 1:gridSize(2), 1:gridSize(3));
    end
    have_density = ~isempty(density) && isequal(size(density), gridSize);

    [~, max_idx] = max(original(:));
    [cx, cy, cz] = ind2sub(gridSize, max_idx);

    max_dose = max(original(:));
    if max_dose == 0, max_dose = 1; end

    x_ax = (1:gridSize(1)) * spacing_mm(1);
    y_ax = (1:gridSize(2)) * spacing_mm(2);
    z_ax = (1:gridSize(3)) * spacing_mm(3);

    cmap_jet = jet(256);
    wl_center = 1050; wl_width = 350;
    wl_min    = wl_center - wl_width / 2;

    figure('Name', titleStr, 'Color', 'w', 'NumberTitle', 'off', ...
        'Position', [50, 50, 1380, 700]);
    sgtitle(sprintf('%s\nIsocenter (max dose): X=%d  Y=%d  Z=%d voxel  |  Dose clim [0, %.4f] Gy', ...
        titleStr, cx, cy, cz, max_dose), 'FontWeight', 'bold', 'FontSize', 11);

    row_labels = {'Air 343 m/s', 'Air 1480 m/s'};
    doses      = {original, recon};

    for row = 1:2
        d   = doses{row};
        lbl = row_labels{row};

        dose_slices   = { squeeze(d(:, cy, :))',  squeeze(d(cx, :, :))',  squeeze(d(:, :, cz))' };
        sensor_slices = { squeeze(sensor_mask(:, cy, :))', ...
                          squeeze(sensor_mask(cx, :, :))', ...
                          squeeze(sensor_mask(:, :, cz))' };
        xvecs  = { x_ax, y_ax, x_ax };
        yvecs  = { z_ax, z_ax, y_ax };
        xlbls  = { 'X (mm)', 'Y (mm)', 'X (mm)' };
        ylbls  = { 'Z (mm)', 'Z (mm)', 'Y (mm)' };
        tsuffs = { sprintf('Coronal (Y=%d)', cy), ...
                   sprintf('Sagittal (X=%d)', cx), ...
                   sprintf('Axial (Z=%d)', cz) };
        if have_density
            density_slices = { squeeze(density(:, cy, :))', ...
                               squeeze(density(cx, :, :))', ...
                               squeeze(density(:, :, cz))' };
        end

        for col = 1:3
            ax         = subplot(2, 3, (row-1)*3 + col);
            dose_slice = double(dose_slices{col});
            if smooth_sigma > 0
                dose_slice = imgaussfilt(dose_slice, smooth_sigma);
            end
            xv         = xvecs{col};
            yv         = yvecs{col};

            if have_density
                dn     = (density_slices{col} - wl_min) / wl_width;
                dn     = max(0, min(1, dn));
                bg_rgb = repmat(dn, [1, 1, 3]);
            else
                bg_rgb = zeros([size(dose_slice), 3]);
            end
            image(ax, xv, yv, bg_rgb);
            hold(ax, 'on');

            norm_d   = max(0, min(1, dose_slice / max_dose));
            idx      = max(1, min(256, round(norm_d * 255) + 1));
            sz       = size(dose_slice);
            dose_rgb = reshape(cmap_jet(idx(:), :), [sz, 3]);

            above      = dose_slice >= 0.10 * max_dose;
            ramp       = 0.45 + 0.35 * min(1, ...
                (dose_slice - 0.10*max_dose) / max(0.40*max_dose, 1e-12));
            dose_alpha = above .* ramp;

            h_dose           = image(ax, xv, yv, dose_rgb);
            h_dose.AlphaData = dose_alpha;

            s = sensor_slices{col};
            if any(s(:))
                contour(ax, xv, yv, s, [0.5, 0.5], 'r-', 'LineWidth', 2);
            end
            hold(ax, 'off');

            axis(ax, 'xy'); axis(ax, 'image');
            colormap(ax, cmap_jet);
            caxis(ax, [0, max_dose]);
            cb = colorbar(ax); cb.Label.String = 'Dose (Gy)';
            xlabel(ax, xlbls{col}); ylabel(ax, ylbls{col});
            title(ax, sprintf('%s  %s', lbl, tsuffs{col}));
        end
    end
    drawnow;
end


function plot_gamma_and_error_axial(gamma_results, original, recon, sensor_mask, cz)
%PLOT_GAMMA_AND_ERROR_AXIAL 1x4 axial figure:
%  gamma 10/10 | gamma 5/5 | gamma 3/3 | absolute dose difference.

    gridSize = size(original);
    cz = max(1, min(gridSize(3), cz));
    if ~isequal(size(sensor_mask), gridSize)
        sensor_mask = sensor_mask(1:gridSize(1), 1:gridSize(2), 1:gridSize(3));
    end

    eval_mask  = gamma_results.eval_mask;
    maps       = gamma_results.maps;
    criteria   = gamma_results.criteria;
    pass_rates = gamma_results.pass_rates;
    nCrit      = size(criteria, 1);

    figure('Name', 'Gamma & Absolute Difference  Axial', 'Color', 'w', ...
        'NumberTitle', 'off', 'Position', [50, 300, 1400, 370]);
    sgtitle(sprintf('Axial Plane (Z = %d voxel)  Gamma (343 vs 1480) & |Difference|', cz), ...
        'FontWeight', 'bold', 'FontSize', 11);

    gamma_clim   = [0, 2];
    sensor_slice = squeeze(sensor_mask(:, :, cz))';

    for g = 1:nCrit
        ax = subplot(1, 4, g);
        gmap = maps{g};
        if ~isempty(gmap)
            gmap_disp = double(gmap);
            gmap_disp(~eval_mask) = NaN;
            slice2d = squeeze(gmap_disp(:, :, cz))';
        else
            slice2d = nan(gridSize(2), gridSize(1));
        end
        imagesc(ax, slice2d, gamma_clim);
        axis(ax, 'xy'); axis(ax, 'image');
        colormap(ax, gamma_colormap_gyr());
        cb = colorbar(ax); cb.Label.String = '\gamma';
        caxis(ax, gamma_clim);
        hold(ax, 'on');
        if any(sensor_slice(:))
            contour(ax, sensor_slice, [0.5, 0.5], 'r-', 'LineWidth', 2);
        end
        hold(ax, 'off');
        pr_str = 'FAILED';
        if ~isnan(pass_rates(g))
            pr_str = sprintf('%.1f%%', pass_rates(g));
        end
        xlabel(ax, 'X (voxel)'); ylabel(ax, 'Y (voxel)');
        title(ax, sprintf('%s\nPass rate: %s', criteria{g, 3}, pr_str));
    end

    ax = subplot(1, 4, 4);
    orig_slice  = squeeze(original(:, :, cz))';
    recon_slice = squeeze(recon(:, :, cz))';
    err_slice   = abs(recon_slice - orig_slice);
    max_err     = max(err_slice(:));
    if max_err == 0, max_err = 1; end
    imagesc(ax, err_slice, [0, max_err]);
    axis(ax, 'xy'); axis(ax, 'image');
    colormap(ax, 'hot');
    cb = colorbar(ax); cb.Label.String = '|Diff| (Gy)';
    hold(ax, 'on');
    if any(sensor_slice(:))
        contour(ax, sensor_slice, [0.5, 0.5], 'r-', 'LineWidth', 2);
    end
    hold(ax, 'off');
    xlabel(ax, 'X (voxel)'); ylabel(ax, 'Y (voxel)');
    title(ax, sprintf('|1480 - 343|\nMax: %.4f Gy', max_err));
    drawnow;
end


function cmap = gamma_colormap_gyr()
%GAMMA_COLORMAP_GYR Green-yellow-red colormap for gamma index.
    n    = 256;
    half = round(n / 2);
    rest = n - half;
    r = [linspace(0, 1, half)'; ones(rest, 1)];
    g = [linspace(0.85, 1, half)'; linspace(1, 0, rest)'];
    b = zeros(n, 1);
    cmap = [r, g, b];
end


function tables = define_tissue_tables()
%DEFINE_TISSUE_TABLES Tissue property lookup tables for HU thresholding.

    tables.threshold_1 = struct();
    tables.threshold_1.hu_boundaries = [-1000, -900, -500, -200, -50, 13, 50, 80, 300, 3000, Inf];
    tables.threshold_1.tissue_names  = {'Air','Lung','Fat','Water','Blood','Muscle','SoftTissue','Bone','Metal'};
    tables.threshold_1.density       = [1.2,  400,  920, 1000, 1060, 1050, 1040, 1900, 7800];
    tables.threshold_1.sound_speed   = [343,  600, 1450, 1480, 1575, 1580, 1540, 3200, 5900];
    tables.threshold_1.alpha_coeff   = [0,   0.5, 0.48, 0.0022, 0.2, 0.5, 0.5,  10,   0];
    tables.threshold_1.alpha_power   = [1.0, 1.5, 1.5,  2.0,   1.3, 1.0, 1.1,  1.0,  1.0];
    tables.threshold_1.gruneisen     = [0,   0.5, 0.7,  0.11,  0.15, 0.2, 1.0,  0,    0];

    tables.threshold_2 = struct();
    tables.threshold_2.hu_boundaries = [-1000, -200, -50, 100, Inf];
    tables.threshold_2.tissue_names  = {'Water','Fat','SoftTissue','Bone'};
    tables.threshold_2.density       = [1000,   920, 1040,         1900];
    tables.threshold_2.sound_speed   = [1480,  1450, 1540,         3200];
    tables.threshold_2.alpha_coeff   = [0.0022, 0.48, 0.5,         10];
    tables.threshold_2.alpha_power   = [2.0,    1.5,  1.1,         1.0];
    tables.threshold_2.gruneisen     = [0.11,   0.7,  1.0,         1.0];
    tables.threshold_2.air_hu_threshold = -300;
    tables.threshold_2.air = struct( ...
        'density',     1.2, ...
        'sound_speed', 343, ...
        'alpha_coeff', 0, ...
        'alpha_power', 1.0, ...
        'gruneisen',   0);
end
