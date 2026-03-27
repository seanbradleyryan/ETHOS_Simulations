%% =========================================================================
%  RUN_STANDALONE_SIMULATION.m
%  Standalone k-Wave Photoacoustic Forward + Time-Reversal Simulation
%  =========================================================================

clear; clc; close all;

%% ========================= CONFIGURATION ================================

CONFIG.working_dir    = '/mnt/weka/home/80030361/ETHOS_Simulations';
CONFIG.patient_id     = '1194203';
CONFIG.session        = 'Session_1';

CONFIG.dose_filename = 'total_rs_dose.mat';
CONFIG.sct_filename  = 'sct_resampled.mat';

CONFIG.dose_file_override = 'field1dose.mat';
CONFIG.sct_file_override  = '';

CONFIG.sensor_placement_method = 'full_plane_lateral';
CONFIG.sensor_x_index = 20;
CONFIG.sensor_y_index = 40;

CONFIG.gruneisen_method = 'uniform';

CONFIG.force_uniform_density     = false;
CONFIG.force_uniform_sound_speed = false;
CONFIG.force_uniform_attenuation = false;
CONFIG.force_uniform_gruneisen   = false;

CONFIG.uniform_density      = 1000;
CONFIG.uniform_sound_speed  = 1540;
CONFIG.uniform_alpha_coeff  = 0;
CONFIG.uniform_alpha_power  = 1.1;
CONFIG.uniform_gruneisen    = 1.0;

CONFIG.dose_per_pulse_cGy     = 0.16;
CONFIG.meterset               = 140;
CONFIG.pml_size               = 10;
CONFIG.cfl_number             = 0.3;
CONFIG.use_gpu                = true;
CONFIG.correction_factor      = 1.9;

CONFIG.num_time_reversal_iter = 30;
CONFIG.convergence_tol        = 1e-3;

CONFIG.use_psf_correction      = false;
CONFIG.regularization_lambda   = 0.05;

CONFIG.downscale_factor = 1;
CONFIG.use_grid_padding = true;

CONFIG.save_results = true;
CONFIG.output_file  = 'standalone_recon_results.mat';
CONFIG.plot_results = true;

%% ========================= RESOLVE FILE PATHS ============================

if ~isempty(CONFIG.dose_file_override)
    dose_filepath = CONFIG.dose_file_override;
else
    processed_dir = fullfile(CONFIG.working_dir, 'RayStationFiles', ...
        CONFIG.patient_id, CONFIG.session, 'processed');
    dose_filepath = fullfile(processed_dir, CONFIG.dose_filename);
end

if ~isempty(CONFIG.sct_file_override)
    sct_filepath = CONFIG.sct_file_override;
else
    if ~exist('processed_dir', 'var')
        processed_dir = fullfile(CONFIG.working_dir, 'RayStationFiles', ...
            CONFIG.patient_id, CONFIG.session, 'processed');
    end
    sct_filepath = fullfile(processed_dir, CONFIG.sct_filename);
end

%% ========================= PRINT CONFIGURATION ===========================

fprintf('=========================================================\n');
fprintf('  Standalone k-Wave Photoacoustic Simulation  (v4.1)\n');
fprintf('=========================================================\n');
fprintf('  Patient:         %s / %s\n', CONFIG.patient_id, CONFIG.session);
fprintf('  Dose file:       %s\n', dose_filepath);
fprintf('  SCT file:        %s\n', sct_filepath);
fprintf('  Sensor:          %s\n', CONFIG.sensor_placement_method);
fprintf('  Tissue model:    %s\n', CONFIG.gruneisen_method);
fprintf('  TR iterations:   %d (tol: %.1e)\n', CONFIG.num_time_reversal_iter, CONFIG.convergence_tol);
fprintf('  GPU:             %s\n', mat2str(CONFIG.use_gpu));
if CONFIG.downscale_factor ~= 1
    fprintf('  Downscale factor: %g\n', CONFIG.downscale_factor);
end
fprintf('=========================================================\n\n');

%% ========================= LOAD DATA ====================================

fprintf('[1/7] Loading dose data...\n');
if ~isfile(dose_filepath)
    error('Dose file not found: %s', dose_filepath);
end
dose_data = load(dose_filepath);

dose_fields = fieldnames(dose_data);
if isfield(dose_data, 'total_rs_dose')
    doseGrid = dose_data.total_rs_dose;
    fprintf('       Loaded variable: total_rs_dose\n');
elseif isfield(dose_data, 'dose_Gy')
    doseGrid = dose_data.dose_Gy;
    fprintf('       Loaded variable: dose_Gy\n');
elseif length(dose_fields) == 1
    doseGrid = dose_data.(dose_fields{1});
    fprintf('       Loaded variable: %s\n', dose_fields{1});
else
    error('Cannot auto-detect dose variable.');
end

if ~isnumeric(doseGrid) || ndims(doseGrid) ~= 3
    error('Dose data must be a 3D numeric array.');
end

gridSize = size(doseGrid);
Nx = gridSize(1); Ny = gridSize(2); Nz = gridSize(3);
fprintf('       Grid size: [%d x %d x %d]\n', Nx, Ny, Nz);
fprintf('       Dose range: [%.6f, %.4f] Gy\n', min(doseGrid(:)), max(doseGrid(:)));

fprintf('[2/7] Loading SCT data...\n');
if ~isfile(sct_filepath)
    error('SCT file not found: %s', sct_filepath);
end
sct_data = load(sct_filepath);
if isfield(sct_data, 'sct_resampled')
    sct = sct_data.sct_resampled;
else
    error('sct_resampled variable not found in %s', sct_filepath);
end

required_sct_fields = {'cubeHU', 'spacing'};
for i = 1:length(required_sct_fields)
    if ~isfield(sct, required_sct_fields{i})
        error('sct_resampled missing required field: %s', required_sct_fields{i});
    end
end

spacing_mm = sct.spacing(:)';
dx = spacing_mm(1) / 1000;
dy = spacing_mm(2) / 1000;
dz = spacing_mm(3) / 1000;

sctSize = size(sct.cubeHU);
if ~isequal(gridSize, sctSize)
    error('Dose grid [%d %d %d] does not match SCT grid [%d %d %d].', ...
        Nx, Ny, Nz, sctSize(1), sctSize(2), sctSize(3));
end

fprintf('       Spacing: [%.2f, %.2f, %.2f] mm\n', spacing_mm);
fprintf('       HU range: [%.0f, %.0f]\n', min(sct.cubeHU(:)), max(sct.cubeHU(:)));

% Mask to body
if isfield(sct, 'bodyMask')
    doseGrid = doseGrid .* double(sct.bodyMask);
end

%% ========================= GRID DOWNSCALING ==============================

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
    dx = spacing_mm(1) / 1000;
    dy = spacing_mm(2) / 1000;
    dz = spacing_mm(3) / 1000;
    sct.spacing = spacing_mm;

    Nx = new_Nx; Ny = new_Ny; Nz = new_Nz;
    gridSize = [Nx, Ny, Nz];
end

%% ========================= CREATE ACOUSTIC MEDIUM ========================

fprintf('[3/7] Creating acoustic medium (method: %s)...\n', CONFIG.gruneisen_method);
medium = create_medium(sct, CONFIG);

fprintf('       Density range:     [%.0f, %.0f] kg/m^3\n', min(medium.density(:)), max(medium.density(:)));
fprintf('       Sound speed range: [%.0f, %.0f] m/s\n',    min(medium.sound_speed(:)), max(medium.sound_speed(:)));
fprintf('       Gruneisen range:   [%.4f, %.4f]\n',        min(medium.gruneisen(:)), max(medium.gruneisen(:)));

%% ========================= INITIAL PRESSURE p0 ==========================

% Apply body mask to dose before p0 calculation
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
    warning('No significant dose or zero initial pressure. Aborting.');
    return;
end

%% ========================= OPTIMAL GRID PADDING ==========================

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

%% ========================= SENSOR PLACEMENT ==============================

fprintf('[5/7] Placing sensor (method: %s)...\n', CONFIG.sensor_placement_method);

sensor      = struct();
sensor.mask = zeros(Nx, Ny, Nz);

switch CONFIG.sensor_placement_method
    case 'full_plane_anterior'
        sensor.mask(CONFIG.sensor_x_index, :, :) = 1;
        fprintf('       Sensor: YZ plane at x = %d\n', CONFIG.sensor_x_index);
    case 'full_plane_lateral'
        sensor.mask(:, CONFIG.sensor_y_index, :) = 1;
        fprintf('       Sensor: XZ plane at y = %d\n', CONFIG.sensor_y_index);
    case 'spherical'
        sph_radius  = floor(min([Nx, Ny, Nz]) / 2) - CONFIG.pml_size;
        sensor.mask = makeSphere(Nx, Ny, Nz, sph_radius);
        fprintf('       Sensor: spherical, radius %d voxels\n', sph_radius);
    otherwise
        error('Unknown sensor_placement_method: "%s"', CONFIG.sensor_placement_method);
end

numSensorPts = sum(sensor.mask(:));
fprintf('       Sensor: %d active points\n', numSensorPts);

if numSensorPts == 0
    warning('Sensor mask is empty. Aborting.');
    return;
end

%% ========================= 3D SENSOR vs DOSE VISUALIZATION ==============
%  Displayed before simulation starts — rotate to inspect geometry.

if CONFIG.plot_results
    sensor_vis = logical(sensor.mask(1:Nx_orig, 1:Ny_orig, 1:Nz_orig));
    plot_3d_sensor_dose(double(doseGrid), sensor_vis, spacing_mm, CONFIG);
    fprintf('       [3D sensor visualization displayed — rotate to inspect]\n');
    drawnow;
end

%% ========================= k-WAVE GRID & MEDIUM SETUP ===================

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

fprintf('       dt = %.2e s, Nt = %d, T_sim = %.2e s\n', dt, Nt, simTime);

kmedium             = struct();
kmedium.density     = medium.density;
kmedium.sound_speed = medium.sound_speed;
kmedium.alpha_coeff = 0 * medium.alpha_coeff;
kmedium.alpha_power = 0 * medium.alpha_power;

if CONFIG.use_gpu
    try
        gpuDevice;
        dataCast = 'gpuArray-single';
        fprintf('       Compute: GPU\n');
    catch
        dataCast = 'single';
        fprintf('       Compute: CPU (GPU unavailable)\n');
    end
else
    dataCast = 'single';
    fprintf('       Compute: CPU\n');
end

inputArgs = {'Smooth', false, 'PMLInside', false, 'PMLSize', CONFIG.pml_size, ...
             'DataCast', dataCast, 'PlotSim', false};

%% ========================= PSF CORRECTION FILTER =========================

psf_filter = [];
if CONFIG.use_psf_correction
    fprintf('\n[PSF] Computing PSF correction filter...\n');
    try
        psf_filter = get_psf(doseGrid, sct, medium_orig, CONFIG);
        fprintf('[PSF] Filter ready.\n');
    catch ME
        warning('get_psf failed: %s. Proceeding without PSF correction.', ME.message);
        psf_filter = [];
    end
end

%% ========================= FORWARD SIMULATION ============================

fprintf('\n[7/7] Running k-Wave forward simulation...\n');

source_fwd    = struct();
source_fwd.p0 = p0;

try
    fwd_tic    = tic;
    sensorData = kspaceFirstOrder3D(kgrid, kmedium, source_fwd, sensor, inputArgs{:});
    fwd_time   = toc(fwd_tic);
    fprintf('       Forward complete (%.1f s). Sensor data: [%d x %d]\n', ...
        fwd_time, size(sensorData, 1), size(sensorData, 2));
catch ME
    fprintf('[ERROR] Forward simulation failed: %s\n', ME.message);
    return;
end

sensorData_measured = sensorData;

%% ========================= ITERATIVE TIME-REVERSAL RECONSTRUCTION ========

fprintf('       Running iterative time reversal (%d iterations, tol=%.1e)...\n', ...
    CONFIG.num_time_reversal_iter, CONFIG.convergence_tol);

reconPressure      = zeros(gridSize);
reconPressure_prev = zeros(gridSize);

% Convergence tracking
conv_max_pressure = zeros(CONFIG.num_time_reversal_iter, 1);
conv_rel_change   = nan(CONFIG.num_time_reversal_iter, 1);
num_iters_done    = 0;

% ---- Set up live TR figure ----
% Axial slice through the max-dose voxel (in original, pre-pad coords).
[~, dose_max_idx] = max(doseGrid(:));
[cx_live, cy_live, cz_live] = ind2sub([Nx_orig, Ny_orig, Nz_orig], dose_max_idx);

if CONFIG.plot_results
    fig_live = figure('Name', 'Live TR Reconstruction', 'Color', 'w', ...
        'NumberTitle', 'off', 'Position', [100, 100, 1060, 440]);

    % Panel 1 — initial p0 (axial slice, fixed reference)
    ax_p0 = subplot(1, 3, 1);
    p0_orig_slice = squeeze(p0(1:Nx_orig, 1:Ny_orig, cz_live))';
    imagesc(ax_p0, p0_orig_slice);
    axis(ax_p0, 'xy'); axis(ax_p0, 'image');
    colormap(ax_p0, 'hot'); colorbar(ax_p0);
    clim_p0 = [0, max(p0_orig_slice(:)) + eps];
    caxis(ax_p0, clim_p0);
    xlabel(ax_p0, 'X (voxel)'); ylabel(ax_p0, 'Y (voxel)');
    title(ax_p0, sprintf('Initial p_0   (Z=%d)', cz_live), 'FontWeight', 'bold');

    % Panel 2 — current reconstructed p0 (updates each iteration)
    ax_recon = subplot(1, 3, 2);
    hImg_recon = imagesc(ax_recon, zeros(Ny_orig, Nx_orig));
    axis(ax_recon, 'xy'); axis(ax_recon, 'image');
    colormap(ax_recon, 'hot'); colorbar(ax_recon);
    xlabel(ax_recon, 'X (voxel)'); ylabel(ax_recon, 'Y (voxel)');
    title(ax_recon, 'Reconstructed p_0   (iter 0)', 'FontWeight', 'bold');

    % Panel 3 — live max-pressure convergence
    ax_conv = subplot(1, 3, 3);
    hLine_max = plot(ax_conv, NaN, NaN, 'b-o', 'LineWidth', 1.6, ...
        'MarkerSize', 4, 'MarkerFaceColor', [0.2, 0.4, 1.0]);
    xlabel(ax_conv, 'TR Iteration');
    ylabel(ax_conv, 'Max Reconstructed p_0 (Pa)');
    title(ax_conv, 'Convergence (live)', 'FontWeight', 'bold');
    grid(ax_conv, 'on');
    xlim(ax_conv, [0.5, CONFIG.num_time_reversal_iter + 0.5]);

    sgtitle(fig_live, sprintf( ...
        'Live TR Reconstruction — Axial Z=%d   |   Patient %s', ...
        cz_live, CONFIG.patient_id), 'FontWeight', 'bold', 'FontSize', 11);
    drawnow;
end

% ---- TR iteration loop ----
try
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

        % Record convergence metrics
        conv_max_pressure(tr_iter) = max(reconPressure(:));
        num_iters_done = tr_iter;

        fprintf('       Max pressure: %.4e Pa\n', conv_max_pressure(tr_iter));

        % Convergence check (from iteration 2 onward)
        converged = false;
        if tr_iter > 1
            norm_prev = norm(reconPressure_prev(:));
            if norm_prev > 0
                rel_change = norm(reconPressure(:) - reconPressure_prev(:)) / norm_prev;
            else
                rel_change = Inf;
            end
            conv_rel_change(tr_iter) = rel_change;
            fprintf('       Rel change: %.4e\n', rel_change);
            if rel_change < CONFIG.convergence_tol
                fprintf('       *** Converged at iteration %d ***\n', tr_iter);
                converged = true;
            end
        end

        reconPressure_prev = reconPressure;

        % ---- Update live figure ----
        if CONFIG.plot_results && ishandle(fig_live)
            recon_slice_crop = squeeze( ...
                reconPressure(1:Nx_orig, 1:Ny_orig, cz_live))';
            set(hImg_recon, 'CData', recon_slice_crop);
            caxis(ax_recon, [0, max(recon_slice_crop(:)) + eps]);
            if converged
                title(ax_recon, ...
                    sprintf('Reconstructed p_0   (iter %d — CONVERGED)', tr_iter), ...
                    'FontWeight', 'bold', 'Color', [0, 0.55, 0]);
            else
                title(ax_recon, ...
                    sprintf('Reconstructed p_0   (iter %d / %d)', ...
                    tr_iter, CONFIG.num_time_reversal_iter), 'FontWeight', 'bold');
            end
            set(hLine_max, 'XData', 1:tr_iter, ...
                'YData', conv_max_pressure(1:tr_iter));
            drawnow;
        end

        if converged
            break;
        end

        % Residual correction for next iteration
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
    fprintf('       Reconstructed pressure: [%.2e, %.2e] Pa\n', ...
        min(reconPressure(:)), max(reconPressure(:)));

catch ME
    fprintf('[ERROR] Time reversal failed: %s\n', ME.message);
    return;
end

%% ========================= CROP TO ORIGINAL SIZE =========================

if did_pad
    fprintf('\n[CROP] Restoring dimensions: [%d %d %d] -> [%d %d %d]\n', ...
        Nx, Ny, Nz, Nx_orig, Ny_orig, Nz_orig);
    reconPressure = reconPressure(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
    Nx = Nx_orig; Ny = Ny_orig; Nz = Nz_orig;
    gridSize    = gridSize_orig;
    medium      = medium_orig;
    sensor.mask = sensor.mask(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
end

%% ========================= PSF CORRECTION ================================

psf_applied = false;
if ~isempty(psf_filter) && isstruct(psf_filter) && isfield(psf_filter, 'F') ...
        && ~isempty(psf_filter.F)
    fprintf('       Applying PSF correction...\n');
    P_field       = fftn(reconPressure);
    corrected     = real(ifftn(P_field .* psf_filter.F));
    reconPressure = max(corrected, 0);
    psf_applied   = true;
    fprintf('       Corrected pressure range: [%.2e, %.2e] Pa\n', ...
        min(reconPressure(:)), max(reconPressure(:)));
end

%% ========================= PRESSURE -> DOSE ==============================

fprintf('\n[Post] Converting pressure to dose...\n');

conversionFactor = medium.gruneisen .* medium.density;
conversionFactor(conversionFactor == 0) = 1;

reconDosePerPulse = reconPressure ./ conversionFactor;

body_mask_plot = ones(gridSize);
if isfield(sct, 'bodyMask') && isequal(size(sct.bodyMask), gridSize)
    body_mask_plot = double(sct.bodyMask);
end
recon_dose = reconDosePerPulse * num_pulses .* double(doseMask) .* body_mask_plot;

fprintf('       Reconstructed dose: [%.4f, %.4f] Gy\n', min(recon_dose(:)), max(recon_dose(:)));

%% Crop p0 to original size (if padded)
if did_pad
    p0 = p0(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
end

%% ========================= RESULTS SUMMARY ===============================

fprintf('\n========= RESULTS =========\n');
fprintf('  Original dose:      [%.6f, %.4f] Gy\n', min(doseGrid(:)), max(doseGrid(:)));
fprintf('  Reconstructed dose: [%.6f, %.4f] Gy\n', min(recon_dose(:)), max(recon_dose(:)));
fprintf('  PSF correction:     %s\n', mat2str(psf_applied));

dose_region = doseGrid > doseThreshold;
if any(dose_region(:))
    abs_error = abs(recon_dose(dose_region) - doseGrid(dose_region));
    rel_error = abs_error ./ max(doseGrid(dose_region), 1e-10) * 100;
    fprintf('  Mean abs err: %.6f Gy\n', mean(abs_error));
    fprintf('  Mean rel err: %.2f%%\n',  mean(rel_error));
    fprintf('  Max  rel err: %.2f%%\n',  max(rel_error));
end
fprintf('===========================\n');

%% ========================= GAMMA ANALYSIS ================================

gamma_results = struct();

if exist('CalcGamma', 'file') == 2

    fprintf('\n[Gamma] Running gamma analysis...\n');

    ref_struct.start = [0, 0, 0];
    ref_struct.width = spacing_mm;
    ref_struct.data  = double(doseGrid);

    tgt_struct.start = [0, 0, 0];
    tgt_struct.width = spacing_mm;
    tgt_struct.data  = double(recon_dose);

    low_dose_cutoff = 0.10 * max(doseGrid(:));
    gamma_eval_mask = doseGrid >= low_dose_cutoff;

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
            eval_vals       = gmap(gamma_eval_mask);
            pass_rate       = 100 * mean(eval_vals <= 1);
            pass_rates(gc)  = pass_rate;
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

    fprintf('\n  ------ Gamma Pass Rates (10%% low-dose cutoff) ------\n');
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
    results.recon_dose    = recon_dose;
    results.original_dose = doseGrid;
    results.p0            = p0;
    results.reconPressure = reconPressure;
    results.sensor_mask   = sensor.mask;
    results.config        = CONFIG;
    results.spacing_mm    = spacing_mm;
    results.grid_size     = gridSize;
    results.fwd_time_sec  = fwd_time;
    results.tr_time_sec   = tr_time;
    results.psf_applied   = psf_applied;
    if ~isempty(gamma_results)
        results.gamma = gamma_results;
    end
    save(CONFIG.output_file, '-struct', 'results', '-v7.3');
    fprintf('\nResults saved to: %s\n', CONFIG.output_file);
end

%% ========================= POST-SIMULATION VISUALIZATION ================

if CONFIG.plot_results

    % Figure 1 — 2x3 dose comparison (original top, recon bottom)
    %            Coronal | Sagittal | Axial
    %            Isocenter = max-dose voxel; sensor contour in red
    plot_dose_panels(doseGrid, recon_dose, sensor.mask, spacing_mm, ...
        'Dose Comparison: Original vs Reconstructed');

    % Figure 2 — p0 convergence (max pressure + relative change)
    plot_convergence_history(conv_max_pressure, conv_rel_change, ...
        num_iters_done, CONFIG.convergence_tol);

    % Figure 3 — Axial gamma (3 criteria) + absolute error
    if ~isempty(gamma_results) && isfield(gamma_results, 'maps')
        [~, max_dose_idx] = max(doseGrid(:));
        [~, ~, cz_gamma]  = ind2sub(gridSize, max_dose_idx);
        plot_gamma_and_error_axial(gamma_results, doseGrid, recon_dose, ...
            sensor.mask, cz_gamma);
    end

end

fprintf('\nStandalone simulation complete.\n');


%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

function plot_3d_sensor_dose(dose_3d, sensor_mask_vis, spacing_mm, config)
%PLOT_3D_SENSOR_DOSE Interactive 3D view of sensor placement vs dose volume.
%  Dose shown as semi-transparent isosurface at 30% of max.
%  Sensor shown as a red translucent plane (or scatter for spherical).
%  Yellow star marks the dose maximum.  Drag to rotate.

    figure('Name', '3D Sensor vs Dose Volume', 'Color', 'w', ...
        'NumberTitle', 'off', 'Position', [80, 80, 800, 660]);
    ax = axes;
    hold(ax, 'on');

    [Nx3, Ny3, Nz3] = size(dose_3d);
    max_d = max(dose_3d(:));

    legend_handles = [];
    legend_entries = {};

    % Dose isosurface at 30% of max
    iso_level = 0.30 * max_d;
    if iso_level > 0
        [f, v] = isosurface(dose_3d, iso_level);
        if ~isempty(f)
            v_mm = v .* [spacing_mm(2), spacing_mm(1), spacing_mm(3)];
            h_dose = patch(ax, 'Faces', f, 'Vertices', v_mm);
            h_dose.FaceColor = [0.20, 0.55, 1.0];
            h_dose.EdgeColor = 'none';
            h_dose.FaceAlpha = 0.35;
            legend_handles(end+1) = h_dose;
            legend_entries{end+1} = 'Dose 30% isosurface';
        end
    end

    % Dose max marker
    [~, midx] = max(dose_3d(:));
    [mi, mj, mk] = ind2sub([Nx3, Ny3, Nz3], midx);
    h_max = scatter3(ax, mi * spacing_mm(1), mj * spacing_mm(2), mk * spacing_mm(3), ...
        160, [1, 0.85, 0], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    legend_handles(end+1) = h_max;
    legend_entries{end+1} = 'Dose maximum';

    % Sensor geometry
    h_sensor = [];
    switch lower(config.sensor_placement_method)
        case 'full_plane_anterior'
            xi      = config.sensor_x_index * spacing_mm(1);
            y_range = linspace(0, Ny3 * spacing_mm(2), 30);
            z_range = linspace(0, Nz3 * spacing_mm(3), 30);
            [Ys, Zs] = meshgrid(y_range, z_range);
            Xs = xi * ones(size(Ys));
            h_sensor = surf(ax, Xs, Ys, Zs, 'FaceColor', [0.9, 0.1, 0.1], ...
                'EdgeColor', 'none', 'FaceAlpha', 0.45);

        case 'full_plane_lateral'
            yi      = config.sensor_y_index * spacing_mm(2);
            x_range = linspace(0, Nx3 * spacing_mm(1), 30);
            z_range = linspace(0, Nz3 * spacing_mm(3), 30);
            [Xs, Zs] = meshgrid(x_range, z_range);
            Ys = yi * ones(size(Xs));
            h_sensor = surf(ax, Xs, Ys, Zs, 'FaceColor', [0.9, 0.1, 0.1], ...
                'EdgeColor', 'none', 'FaceAlpha', 0.45);

        case 'spherical'
            [si, sj, sk] = ind2sub(size(sensor_mask_vis), find(sensor_mask_vis));
            h_sensor = scatter3(ax, si * spacing_mm(1), sj * spacing_mm(2), sk * spacing_mm(3), ...
                3, [0.9, 0.1, 0.1], 'filled', 'MarkerFaceAlpha', 0.25);
    end
    if ~isempty(h_sensor)
        legend_handles(end+1) = h_sensor;
        legend_entries{end+1} = sprintf('Sensor (%s)', config.sensor_placement_method);
    end

    camlight(ax, 'headlight');
    lighting(ax, 'gouraud');
    axis(ax, 'equal');
    grid(ax, 'on');
    view(ax, [-38, 24]);
    rotate3d(ax, 'on');

    xlabel(ax, 'X (mm)'); ylabel(ax, 'Y (mm)'); zlabel(ax, 'Z (mm)');
    title(ax, sprintf( ...
        '3D Sensor Placement vs Dose Volume\nSensor: %s   |   Yellow star = dose max   |   Drag to rotate', ...
        config.sensor_placement_method), 'FontSize', 10);

    if ~isempty(legend_handles)
        legend(ax, legend_handles, legend_entries, ...
            'Location', 'northeastoutside', 'FontSize', 8);
    end

    hold(ax, 'off');
end


function plot_dose_panels(original, recon, sensor_mask, spacing_mm, titleStr)
%PLOT_DOSE_PANELS 2x3 dose comparison: coronal, sagittal, axial.
%  Row 1 = original dose,  Row 2 = reconstructed dose.
%  Isocenter at max-dose voxel.  Sensor contour in red on every panel.

    gridSize = size(original);
    if ~isequal(size(sensor_mask), gridSize)
        sensor_mask = sensor_mask(1:gridSize(1), 1:gridSize(2), 1:gridSize(3));
    end

    [~, max_idx] = max(original(:));
    [cx, cy, cz] = ind2sub(gridSize, max_idx);

    max_dose = max(original(:));
    if max_dose == 0, max_dose = 1; end

    x_ax = (1:gridSize(1)) * spacing_mm(1);
    y_ax = (1:gridSize(2)) * spacing_mm(2);
    z_ax = (1:gridSize(3)) * spacing_mm(3);

    figure('Name', titleStr, 'Color', 'w', 'NumberTitle', 'off', ...
        'Position', [50, 50, 1380, 700]);
    sgtitle(sprintf('%s\nIsocenter (max dose): X=%d  Y=%d  Z=%d voxel', ...
        titleStr, cx, cy, cz), 'FontWeight', 'bold', 'FontSize', 11);

    row_labels = {'Original', 'Reconstructed'};
    doses = {original, recon};

    for row = 1:2
        d   = doses{row};
        lbl = row_labels{row};

        % Coronal — XZ plane at y = cy
        ax = subplot(2, 3, (row-1)*3 + 1);
        img = squeeze(d(:, cy, :))';
        imagesc(ax, x_ax, z_ax, img, [0, max_dose]);
        axis(ax, 'xy'); axis(ax, 'image');
        colormap(ax, 'jet');
        cb = colorbar(ax); cb.Label.String = 'Dose (Gy)';
        hold(ax, 'on');
        s = squeeze(sensor_mask(:, cy, :))';
        if any(s(:)), contour(ax, x_ax, z_ax, s, [0.5,0.5], 'r-', 'LineWidth', 2); end
        hold(ax, 'off');
        xlabel(ax, 'X (mm)'); ylabel(ax, 'Z (mm)');
        title(ax, sprintf('%s — Coronal (Y=%d)', lbl, cy));

        % Sagittal — YZ plane at x = cx
        ax = subplot(2, 3, (row-1)*3 + 2);
        img = squeeze(d(cx, :, :))';
        imagesc(ax, y_ax, z_ax, img, [0, max_dose]);
        axis(ax, 'xy'); axis(ax, 'image');
        colormap(ax, 'jet');
        cb = colorbar(ax); cb.Label.String = 'Dose (Gy)';
        hold(ax, 'on');
        s = squeeze(sensor_mask(cx, :, :))';
        if any(s(:)), contour(ax, y_ax, z_ax, s, [0.5,0.5], 'r-', 'LineWidth', 2); end
        hold(ax, 'off');
        xlabel(ax, 'Y (mm)'); ylabel(ax, 'Z (mm)');
        title(ax, sprintf('%s — Sagittal (X=%d)', lbl, cx));

        % Axial — XY plane at z = cz
        ax = subplot(2, 3, (row-1)*3 + 3);
        img = squeeze(d(:, :, cz))';
        imagesc(ax, x_ax, y_ax, img, [0, max_dose]);
        axis(ax, 'xy'); axis(ax, 'image');
        colormap(ax, 'jet');
        cb = colorbar(ax); cb.Label.String = 'Dose (Gy)';
        hold(ax, 'on');
        s = squeeze(sensor_mask(:, :, cz))';
        if any(s(:)), contour(ax, x_ax, y_ax, s, [0.5,0.5], 'r-', 'LineWidth', 2); end
        hold(ax, 'off');
        xlabel(ax, 'X (mm)'); ylabel(ax, 'Y (mm)');
        title(ax, sprintf('%s — Axial (Z=%d)', lbl, cz));
    end
    drawnow;
end


function plot_convergence_history(conv_max_pressure, conv_rel_change, num_iters, tol)
%PLOT_CONVERGENCE_HISTORY p0 convergence over TR iterations.
%  Left axis:  max reconstructed p0 per iteration (blue).
%  Right axis: relative change between iterations (red, log-scale).

    iters   = 1:num_iters;
    p_vals  = conv_max_pressure(iters);
    rc_vals = conv_rel_change(iters);

    figure('Name', 'p0 Convergence', 'Color', 'w', ...
        'NumberTitle', 'off', 'Position', [150, 520, 720, 390]);

    yyaxis left;
    plot(iters, p_vals, 'b-o', 'LineWidth', 1.8, 'MarkerSize', 5, ...
        'MarkerFaceColor', [0.2, 0.4, 1.0]);
    ylabel('Max Reconstructed p_0 (Pa)', 'Color', [0.2, 0.4, 1.0]);
    ylim([0, max(p_vals) * 1.15]);

    yyaxis right;
    valid = ~isnan(rc_vals);
    if any(valid)
        semilogy(iters(valid), rc_vals(valid), 'r-s', 'LineWidth', 1.8, ...
            'MarkerSize', 5, 'MarkerFaceColor', [0.9, 0.1, 0.1]);
        hold on;
        yline(tol, 'k--', sprintf('tol = %.0e', tol), 'LineWidth', 1.2, ...
            'LabelHorizontalAlignment', 'right');
        hold off;
    end
    ylabel('Relative Change ||p_n - p_{n-1}|| / ||p_{n-1}||', ...
        'Color', [0.8, 0.1, 0.1]);

    xlabel('TR Iteration');
    title(sprintf('p_0 Convergence  (%d/%d iterations)', num_iters, numel(conv_max_pressure)), ...
        'FontWeight', 'bold');
    xlim([0.5, num_iters + 0.5]);
    grid on;
    drawnow;
end


function plot_gamma_and_error_axial(gamma_results, original, recon, sensor_mask, cz)
%PLOT_GAMMA_AND_ERROR_AXIAL 1x4 axial figure:
%  gamma 10/10 | gamma 5/5 | gamma 3/3 | absolute dose error.
%  Sensor contour in red.  Slice at cz (max-dose axial index).

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

    figure('Name', 'Gamma & Absolute Error — Axial', 'Color', 'w', ...
        'NumberTitle', 'off', 'Position', [50, 300, 1400, 370]);
    sgtitle(sprintf('Axial Plane (Z = %d voxel) — Gamma Index & Absolute Error', cz), ...
        'FontWeight', 'bold', 'FontSize', 11);

    gamma_clim   = [0, 2];
    sensor_slice = squeeze(sensor_mask(:, :, cz))';

    % Gamma panels
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

    % Absolute error panel
    ax = subplot(1, 4, 4);
    orig_slice  = squeeze(original(:, :, cz))';
    recon_slice = squeeze(recon(:, :, cz))';
    err_slice   = abs(recon_slice - orig_slice);
    max_err     = max(err_slice(:));
    if max_err == 0, max_err = 1; end

    imagesc(ax, err_slice, [0, max_err]);
    axis(ax, 'xy'); axis(ax, 'image');
    colormap(ax, 'hot');
    cb = colorbar(ax); cb.Label.String = '|Error| (Gy)';

    hold(ax, 'on');
    if any(sensor_slice(:))
        contour(ax, sensor_slice, [0.5, 0.5], 'r-', 'LineWidth', 2);
    end
    hold(ax, 'off');

    xlabel(ax, 'X (voxel)'); ylabel(ax, 'Y (voxel)');
    title(ax, sprintf('|Recon - Original|\nMax: %.4f Gy', max_err));

    drawnow;
end


function cmap = gamma_colormap_gyr()
%GAMMA_COLORMAP_GYR Green-yellow-red colormap for gamma index.
%  Green = pass (gamma <= 1),  Red = fail (gamma > 1).
    n    = 256;
    half = round(n / 2);
    rest = n - half;
    r = [linspace(0, 1, half)'; ones(rest, 1)];
    g = [linspace(0.85, 1, half)'; linspace(1, 0, rest)'];
    b = zeros(n, 1);
    cmap = [r, g, b];
end


function medium = create_medium(sct, config)
%CREATE_MEDIUM Build acoustic medium from SCT data and tissue model config.

    HU       = double(sct.cubeHU);
    gridSize = size(HU);
    tables   = define_tissue_tables();

    switch lower(config.gruneisen_method)
        case 'uniform'
            medium.density     = ones(gridSize) * config.uniform_density;
            medium.sound_speed = ones(gridSize) * config.uniform_sound_speed;
            medium.alpha_coeff = ones(gridSize) * config.uniform_alpha_coeff;
            medium.alpha_power = config.uniform_alpha_power;
            medium.gruneisen   = ones(gridSize) * config.uniform_gruneisen;

        case {'threshold_1', 'threshold_2'}
            T          = tables.(config.gruneisen_method);
            nTissues   = length(T.tissue_names);
            boundaries = T.hu_boundaries;

            medium.density     = ones(gridSize) * 1000;
            medium.sound_speed = ones(gridSize) * 1540;
            medium.alpha_coeff = ones(gridSize) * 0.5;
            medium.alpha_power = T.alpha_power(1);
            medium.gruneisen   = ones(gridSize) * 0.11;

            for t = 1:nTissues
                mask = (HU >= boundaries(t)) & (HU < boundaries(t+1));
                medium.density(mask)     = T.density(t);
                medium.sound_speed(mask) = T.sound_speed(t);
                medium.alpha_coeff(mask) = T.alpha_coeff(t);
                medium.gruneisen(mask)   = T.gruneisen(t);
            end

            st_idx = find(contains(lower(T.tissue_names), 'soft'), 1);
            if ~isempty(st_idx)
                medium.alpha_power = T.alpha_power(st_idx);
            end

            fprintf('       Tissue model: %s (%d tissues)\n', config.gruneisen_method, nTissues);
            for t = 1:nTissues
                mask = (HU >= boundaries(t)) & (HU < boundaries(t+1));
                fprintf('         %-12s: %7d voxels (%.1f%%)\n', ...
                    T.tissue_names{t}, sum(mask(:)), 100 * sum(mask(:)) / numel(HU));
            end

        otherwise
            error('Unknown gruneisen_method: %s', config.gruneisen_method);
    end

    if config.force_uniform_density
        medium.density = ones(gridSize) * config.uniform_density;
    end
    if config.force_uniform_sound_speed
        medium.sound_speed = ones(gridSize) * config.uniform_sound_speed;
    end
    if config.force_uniform_attenuation
        medium.alpha_coeff = ones(gridSize) * config.uniform_alpha_coeff;
        medium.alpha_power = config.uniform_alpha_power;
    end
    if config.force_uniform_gruneisen
        medium.gruneisen = ones(gridSize) * config.uniform_gruneisen;
    end

    if isfield(sct, 'bodyMask')
        outside_body = ~logical(sct.bodyMask);
        if isfield(sct, 'couchMask')
            outside_body = outside_body | logical(sct.couchMask);
        end
        medium.density(outside_body)     = config.uniform_density;
        medium.sound_speed(outside_body) = config.uniform_sound_speed;
        medium.alpha_coeff(outside_body) = config.uniform_alpha_coeff;
        medium.gruneisen(outside_body)   = config.uniform_gruneisen;
    end

    medium.grid_size = gridSize;
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
end