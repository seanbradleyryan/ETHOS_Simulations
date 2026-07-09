%% RUN_MEDIUM_COMPARISON.m
%  Photoacoustic reconstruction comparison across 4 tissue media.
%  Media layers stack along the Z axis (k-Wave dim 3).
%  Each medium is simulated with both a spherical and a planar (XY-plane)
%  sensor, reconstructed via iterative TR (<=10 iter, 1% rel-change tol).
%
%  Media:
%    (a) Homogeneous oil
%    (b) Half oil / half water  — planar sensor on oil face (z = 2)
%    (c) Half oil / half water  — planar sensor on water face (z = Nz-1)
%    (d) Third oil / third bone / third water — planar sensor on oil face
%
%  Per-medium figure layout  (3 rows x 3 columns):
%    Row 1 : [True p0]          [Spherical recon]    [Planar recon]
%    Row 2 : [Medium schematic] [Spherical gamma]    [Planar gamma]
%    Row 3 : [Convergence history — spanning all three columns]
%
%  Dependencies: k-Wave toolbox, CalcGamma.m (must be on MATLAB path)
%
%  Author: ETHOS Pipeline Team   |  Date: March 2026

clear; clc; close all;

fprintf('=====================================================\n');
fprintf('  Multi-Medium Photoacoustic Reconstruction Study\n');
fprintf('=====================================================\n\n');

%% ========================= CONFIGURATION ================================

CONFIG.working_dir        = '/mnt/weka/home/80030361/ETHOS_Simulations';
CONFIG.patient_id         = '1194203';
CONFIG.session            = 'Session_1';
CONFIG.dose_filename      = 'total_rs_dose.mat';

CONFIG.downscale_factor   = 2;      % Spatial downscale for speed (1 = off)
CONFIG.meterset           = 140;    % Monitor units
CONFIG.dose_per_pulse_cGy = 0.16;  % cGy per LINAC pulse
CONFIG.pml_size           = 10;    % PML thickness in voxels
CONFIG.cfl_number         = 0.3;   % CFL stability number
CONFIG.Nt_scaling         = 6;     % >0: when air sets minC, divide Nt by this to shorten recording (0 = off)
CONFIG.use_gpu            = true;  % GPU acceleration
CONFIG.max_tr_iter        = 10;    % Max iterative TR iterations
CONFIG.conv_tol           = 0.01;  % 1% relative-change convergence threshold

% Gamma analysis parameters
CONFIG.gamma_pct    = 3;    % Dose difference criterion (%)
CONFIG.gamma_dta    = 3;    % Distance-to-agreement (mm)
CONFIG.gamma_cutoff = 10;   % Dose cutoff (% of max) for pass-rate mask

% Acoustic tissue properties: density (kg/m^3), sound_speed (m/s), Gruneisen, RGB colour
PROPS.oil   = struct('name','Oil',   'density',  900, 'sound_speed', 1460, ...
                     'gruneisen', 0.70, 'color', [1.00, 0.85, 0.00]);
PROPS.water = struct('name','Water', 'density', 1000, 'sound_speed', 1480, ...
                     'gruneisen', 0.11, 'color', [0.20, 0.45, 0.90]);
PROPS.bone  = struct('name','Bone',  'density', 1900, 'sound_speed', 3200, ...
                     'gruneisen', 0.50, 'color', [0.75, 0.75, 0.75]);

%% ========================= LOAD DATA ====================================

fprintf('[Load] Reading processed data...\n');
processed_dir = fullfile(CONFIG.working_dir, 'RayStationFiles', ...
    CONFIG.patient_id, CONFIG.session, 'processed');

dose_loaded = load(fullfile(processed_dir, CONFIG.dose_filename));
if     isfield(dose_loaded, 'total_rs_dose'), doseGrid = double(dose_loaded.total_rs_dose);
elseif isfield(dose_loaded, 'dose_Gy'),       doseGrid = double(dose_loaded.dose_Gy);
else,  f = fieldnames(dose_loaded); doseGrid = double(dose_loaded.(f{1})); end

sct_loaded = load(fullfile(processed_dir, 'sct_resampled.mat'));
spacing_mm = double(sct_loaded.sct_resampled.spacing(:)');

fprintf('[Load] Dose: [%d x %d x %d]  max=%.4f Gy   Spacing: [%.2f %.2f %.2f] mm\n', ...
    size(doseGrid), max(doseGrid(:)), spacing_mm);

%% ========================= OPTIONAL DOWNSCALE ===========================

if CONFIG.downscale_factor ~= 1
    df  = CONFIG.downscale_factor;
    sz0 = size(doseGrid);
    sz1 = max([1 1 1], round(sz0 / df));
    fprintf('[Downscale] [%d %d %d] -> [%d %d %d] (factor %g)\n', sz0, sz1, df);
    doseGrid   = max(imresize3(doseGrid, sz1), 0);
    spacing_mm = spacing_mm .* double(sz0 ./ sz1);
end

[Nx, Ny, Nz] = size(doseGrid);
dx = spacing_mm(1) / 1000;
dy = spacing_mm(2) / 1000;
dz = spacing_mm(3) / 1000;
gridSize = [Nx, Ny, Nz];
fprintf('[Grid] %d x %d x %d   spacing: [%.2f %.2f %.2f] mm\n', ...
    Nx, Ny, Nz, spacing_mm);

%% ========================= DOSE CENTROID ================================

d_thresh = 0.01 * max(doseGrid(:));
d_mask   = doseGrid > d_thresh;
d_sum    = max(sum(doseGrid(d_mask)), 1e-30);
[xi, yi, zi] = ndgrid(1:Nx, 1:Ny, 1:Nz);
cx = round(sum(xi(d_mask) .* doseGrid(d_mask)) / d_sum);
cy = round(sum(yi(d_mask) .* doseGrid(d_mask)) / d_sum);
cz = round(sum(zi(d_mask) .* doseGrid(d_mask)) / d_sum);
cx = max(1, min(Nx, cx));
cy = max(1, min(Ny, cy));
cz = max(1, min(Nz, cz));
fprintf('[Centroid] [%d, %d, %d]\n\n', cx, cy, cz);

%% ========================= GPU / kWave INPUT ARGS =======================

if CONFIG.use_gpu
    try
        gpuDevice;
        dataCast = 'gpuArray-single';
        fprintf('[GPU] Using GPU.\n');
    catch
        dataCast = 'single';
        fprintf('[GPU] Not available; using CPU.\n');
    end
else
    dataCast = 'single';
    fprintf('[GPU] CPU mode.\n');
end

inputArgs = {'Smooth', false, 'PMLInside', false, ...
             'PMLSize', CONFIG.pml_size, 'DataCast', dataCast, 'PlotSim', false};

num_pulses = ceil(CONFIG.meterset / CONFIG.dose_per_pulse_cGy);
dose_pp    = doseGrid / num_pulses;  % per-pulse dose (Gy)

%% ========================= MEDIA SPECIFICATIONS (Z-plane layers) ========
%
%  Each medium varies along Z (k-Wave dimension 3).
%  Layers appear as horizontal bands in the XZ view.
%
%  Columns:
%    1  title string
%    2  Z-region index ranges   (cell of vectors)
%    3  property structs         (cell matching col 2)
%    4  planar sensor Z index
%    5  planar sensor label

split_h  = floor(Nz / 2);         % half-grid Z boundary
split_t1 = floor(Nz / 3);         % first-third Z boundary
split_t2 = floor(2 * Nz / 3);     % second-third Z boundary

media_specs = { ...
    '(a) Homogeneous Oil', ...
        {1:Nz}, ...
        {PROPS.oil}, ...
        2, 'z=2 (oil face)'; ...
    ...
    '(b) Half Oil / Half Water  |  Sensor: oil face (z=2)', ...
        {1:split_h, split_h+1:Nz}, ...
        {PROPS.oil, PROPS.water}, ...
        2, 'z=2 (oil)'; ...
    ...
    '(c) Half Oil / Half Water  |  Sensor: water face', ...
        {1:split_h, split_h+1:Nz}, ...
        {PROPS.oil, PROPS.water}, ...
        Nz-1, sprintf('z=%d (water)', Nz-1); ...
    ...
    '(d) Oil / Bone / Water  |  Sensor: oil face (z=2)', ...
        {1:split_t1, split_t1+1:split_t2, split_t2+1:Nz}, ...
        {PROPS.oil, PROPS.bone, PROPS.water}, ...
        2, 'z=2 (oil)'; ...
};

%% ========================= SPHERICAL SENSOR RADIUS ======================
%
%  The sphere is centered at floor(N/2) in each dimension.
%  Maximum safe radius = distance from centre to PML inner edge - 1.
%  Using the MINIMUM over all three dimensions to guarantee full containment.

sph_r = min([ ...
    floor(Nx/2) - CONFIG.pml_size - 1, ...
    floor(Ny/2) - CONFIG.pml_size - 1, ...
    floor(Nz/2) - CONFIG.pml_size - 1]) - 1;
sph_r = max(sph_r, 5);  % safety floor
fprintf('[Sphere] Radius = %d voxels\n\n', sph_r);

%% ========================= MAIN SIMULATION LOOP =========================

for m = 1:4
    m_title  = media_specs{m, 1};
    m_regs   = media_specs{m, 2};
    m_props  = media_specs{m, 3};
    m_pz     = media_specs{m, 4};
    m_plabel = media_specs{m, 5};

    fprintf('\n+-------------------------------------------------------+\n');
    fprintf('|  MEDIUM %d: %s\n', m, m_title);
    fprintf('+-------------------------------------------------------+\n');

    %% Build spatially-varying medium arrays (layers along Z)
    density_a = zeros(gridSize);
    speed_a   = zeros(gridSize);
    gruen_a   = zeros(gridSize);
    label_a   = zeros(gridSize, 'uint8');

    for r = 1:numel(m_regs)
        zi_r = m_regs{r};
        p    = m_props{r};
        density_a(:, :, zi_r) = p.density;
        speed_a(:, :, zi_r)   = p.sound_speed;
        gruen_a(:, :, zi_r)   = p.gruneisen;
        label_a(:, :, zi_r)   = uint8(r);
    end

    kmedium = struct('density',     density_a, ...
                     'sound_speed', speed_a, ...
                     'alpha_coeff', 0, ...
                     'alpha_power', 1.1);

    %% Initial pressure:  p0 = (dose / num_pulses) * Gruneisen * rho
    p0 = double(dose_pp .* gruen_a .* density_a);
    p0 = double(smooth(p0));  % k-Wave Hann-window smoothing
    p0 = max(p0, 0);
    fprintf('  p0 range: [%.2e, %.2e] Pa\n', min(p0(:)), max(p0(:)));

    %% CFL-stable time step for the fastest sound speed in this medium
    maxC    = max(speed_a(:));
    minC    = min(speed_a(speed_a > 0));
    dt      = CONFIG.cfl_number * min([dx, dy, dz]) / maxC;
    diagLen = sqrt((Nx*dx)^2 + (Ny*dy)^2 + (Nz*dz)^2);
    Nt      = ceil(2.5 * diagLen / minC / dt);

    % Nt_scaling: air (~343 m/s) is the slowest tissue, so it sets minC and
    % inflates Nt. When air drives minC, shorten the recording length by
    % CONFIG.Nt_scaling (0 = disabled).
    if isfield(CONFIG, 'Nt_scaling') && CONFIG.Nt_scaling ~= 0
        air_c = 343;
        if isfield(CONFIG, 'tissue_tables') && isfield(CONFIG.tissue_tables, 'threshold_2') ...
                && isfield(CONFIG.tissue_tables.threshold_2, 'air')
            air_c = CONFIG.tissue_tables.threshold_2.air.sound_speed;
        end
        if minC <= air_c
            Nt = max(1, ceil(Nt / CONFIG.Nt_scaling));
        end
    end

    kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
    kgrid.dt = dt;
    kgrid.Nt = Nt;
    fprintf('  dt = %.2e s,  Nt = %d\n', dt, Nt);

    %% Sensor masks
    % Spherical: centred shell of radius sph_r (fully encloses source region)
    sensor_sph.mask  = makeSphere(Nx, Ny, Nz, sph_r);

    % Planar: full XY plane at z = m_pz
    sensor_plan.mask = false(gridSize);
    sensor_plan.mask(:, :, m_pz) = true;

    fprintf('  Spherical: r=%d, %d voxels  |  Planar: z=%d, %d voxels\n', ...
        sph_r, sum(sensor_sph.mask(:)), m_pz, sum(sensor_plan.mask(:)));

    source_fwd.p0 = p0;

    %% ------ Spherical sensor: forward + iterative TR -------------------
    fprintf('\n  [Spherical] Forward simulation...\n');
    t0 = tic;
    sd_sph = kspaceFirstOrder3D(kgrid, kmedium, source_fwd, sensor_sph, inputArgs{:});
    fprintf('  [Spherical] Forward: %.1f s.  Iterative TR...\n', toc(t0));

    [p0r_sph, conv_sph] = run_iterative_tr(kgrid, kmedium, sensor_sph, sd_sph, ...
        gridSize, CONFIG.max_tr_iter, CONFIG.conv_tol, inputArgs);

    p0r_sph = double(gather(p0r_sph));
    clear sd_sph

    %% ------ Planar sensor: forward + iterative TR ----------------------
    fprintf('\n  [Planar / %s] Forward simulation...\n', m_plabel);
    t0 = tic;
    sd_plan = kspaceFirstOrder3D(kgrid, kmedium, source_fwd, sensor_plan, inputArgs{:});
    fprintf('  [Planar] Forward: %.1f s.  Iterative TR...\n', toc(t0));

    [p0r_plan, conv_plan] = run_iterative_tr(kgrid, kmedium, sensor_plan, sd_plan, ...
        gridSize, CONFIG.max_tr_iter, CONFIG.conv_tol, inputArgs);

    p0r_plan = double(gather(p0r_plan));
    clear sd_plan

    %% ------ Gamma analysis (3D, restricted search) ---------------------
    fprintf('  Computing gamma (%.0f%%/%.0fmm)...\n', ...
        CONFIG.gamma_pct, CONFIG.gamma_dta);

    ref_s.start = [0, 0, 0];
    ref_s.width = spacing_mm;
    ref_s.data  = p0;

    tar_sph.start = [0, 0, 0];
    tar_sph.width = spacing_mm;
    tar_sph.data  = p0r_sph;

    tar_plan.start = [0, 0, 0];
    tar_plan.width = spacing_mm;
    tar_plan.data  = p0r_plan;

    gamma_sph  = CalcGamma(ref_s, tar_sph,  CONFIG.gamma_pct, CONFIG.gamma_dta, ...
        'local', 0, 'restrict', 1, 'res', 20, 'limit', 2, 'cpu', 1);
    gamma_plan = CalcGamma(ref_s, tar_plan, CONFIG.gamma_pct, CONFIG.gamma_dta, ...
        'local', 0, 'restrict', 1, 'res', 20, 'limit', 2, 'cpu', 1);

    % Pass rates (dose-thresholded)
    p0_max   = max(p0(:));
    p0_cutoff = p0_max * CONFIG.gamma_cutoff / 100;
    gam_mask  = p0 >= p0_cutoff;
    pr_sph    = 100 * sum(gamma_sph(gam_mask) <= 1) / max(sum(gam_mask(:)), 1);
    pr_plan   = 100 * sum(gamma_plan(gam_mask) <= 1) / max(sum(gam_mask(:)), 1);
    fprintf('  Gamma pass rates:  Spherical = %.1f%%   Planar = %.1f%%\n', ...
        pr_sph, pr_plan);

    %% ------ Generate figure --------------------------------------------
    make_comparison_figure(m, m_title, ...
        p0, p0r_sph, p0r_plan, ...
        sensor_sph.mask, m_pz, ...
        gamma_sph, gamma_plan, gam_mask, pr_sph, pr_plan, ...
        conv_sph, conv_plan, ...
        label_a, m_props, ...
        cy, cx, cz, gridSize, spacing_mm, CONFIG);

    clear density_a speed_a gruen_a label_a p0 p0r_sph p0r_plan kmedium
    clear gamma_sph gamma_plan gam_mask
    fprintf('\n  Medium %d complete.\n', m);
end

fprintf('\n=====================================================\n');
fprintf('  All simulations complete.\n');
fprintf('=====================================================\n');


%% =========================================================================
%  LOCAL FUNCTION: iterative time-reversal reconstruction
%% =========================================================================

function [p0_out, conv_hist] = run_iterative_tr(kgrid, kmedium, sensor, ...
    sd_meas, gridSize, max_iter, conv_tol, in_args)
%RUN_ITERATIVE_TR  Iterative Dirichlet time-reversal reconstruction.
%
%   Each iteration:
%     1. TR( working_data ) -> p0_est  [Dirichlet BC, positivity constrained]
%     2. Forward( p0_est )  -> sd_recon
%     3. Residual correction: working_data = measured + (measured - sd_recon)
%     4. Convergence: ||p0_n - p0_{n-1}|| / ||p0_{n-1}|| < conv_tol  -> stop

    Nx = gridSize(1);  Ny = gridSize(2);  Nz = gridSize(3);
    sd_work   = sd_meas;
    p0_prev   = zeros(gridSize, 'single');
    conv_hist = [];
    p0_out    = zeros(gridSize);

    for it = 1:max_iter
        % Time-reversal
        src_tr.p_mask = sensor.mask;
        src_tr.p      = fliplr(sd_work);
        src_tr.p_mode = 'dirichlet';

        sen_tr.mask   = ones(Nx, Ny, Nz);
        sen_tr.record = {'p_final'};

        raw = kspaceFirstOrder3D(kgrid, kmedium, src_tr, sen_tr, in_args{:});

        if isstruct(raw) && isfield(raw, 'p_final')
            p0_it = single(gather(reshape(raw.p_final, gridSize)));
        else
            p0_it = single(gather(reshape(raw, gridSize)));
        end
        p0_it  = max(p0_it, 0);  % positivity constraint
        p0_out = double(p0_it);

        % Convergence check (from iter 2 onward)
        if it > 1
            nn = double(norm(p0_prev(:)));
            rc = double(norm(p0_it(:) - p0_prev(:))) / max(nn, 1e-30);
            conv_hist(end+1) = rc; %#ok<AGROW>
            fprintf('    TR iter %2d   rel_change = %.4e\n', it, rc);
            if rc < conv_tol
                fprintf('    Converged at iteration %d.\n', it);
                return;
            end
        else
            fprintf('    TR iter %2d\n', it);
        end
        p0_prev = p0_it;

        % Residual correction for next iteration
        if it < max_iter
            src_res.p0 = p0_it;
            sd_res     = kspaceFirstOrder3D(kgrid, kmedium, src_res, sensor, in_args{:});
            sd_work    = sd_meas + (sd_meas - sd_res);
        end
    end
end


%% =========================================================================
%  LOCAL FUNCTION: per-medium 3x3 comparison figure
%% =========================================================================

function make_comparison_figure(m_idx, m_title, ...
    p0_true, p0_sph, p0_plan, ...
    mask_sph, pz_sensor, ...
    gamma_sph, gamma_plan, gam_mask, pr_sph, pr_plan, ...
    conv_sph, conv_plan, ...
    label_arr, props_list, ...
    cy, cx, cz, gridSize, spacing_mm, CONFIG)
%MAKE_COMPARISON_FIGURE  3x3 figure for one medium.
%
%  All image panels show XZ slices at Y=cy (dose-centroid Y index).
%  X is the horizontal axis, Z is the vertical axis.
%  Medium layers appear as HORIZONTAL bands (varying in Z).

    Nx = gridSize(1);
    Nz = gridSize(3);
    n_tissues = numel(props_list);

    % XZ slices at Y=cy  =>  [Nx x Nz] matrices
    p0t_xz   = squeeze(p0_true(:, cy, :));
    p0s_xz   = squeeze(p0_sph(:,  cy, :));
    p0p_xz   = squeeze(p0_plan(:, cy, :));
    lbl_xz   = squeeze(label_arr(:, cy, :));
    msph_xz  = squeeze(mask_sph(:, cy, :));  % logical [Nx x Nz]
    gs_xz    = squeeze(gamma_sph(:,  cy, :));
    gp_xz    = squeeze(gamma_plan(:, cy, :));
    gmask_xz = squeeze(gam_mask(:, cy, :));

    % Common pressure colour scale (based on true p0)
    cmax = max(p0t_xz(:));
    if cmax < eps, cmax = 1; end
    crange = [0, cmax * 1.02];

    % Spherical sensor points in this XZ slice
    [sph_ix, sph_iz] = find(msph_xz);

    % Colormaps
    hot_cmap = hot(256);
    gam_cmap = gamma_colormap();

    % Create figure
    fig = figure('Name', sprintf('Medium %d', m_idx), 'NumberTitle', 'off', ...
        'Position', [20 + 25*m_idx, 20 + 15*m_idx, 1550, 900], 'Color', 'w');
    sgtitle(sprintf('Medium %d: %s', m_idx, m_title), ...
        'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'none');

    % ---- (1,1)  True p0 ------------------------------------------------
    ax11 = subplot(3, 3, 1);
    imagesc(ax11, 1:Nx, 1:Nz, p0t_xz');
    colormap(ax11, hot_cmap);
    caxis(ax11, crange);  colorbar(ax11);
    set(ax11, 'YDir', 'normal');
    hold(ax11, 'on');
    % Mark dose centroid
    plot(ax11, cx, cz, 'c+', 'MarkerSize', 12, 'LineWidth', 2);
    hold(ax11, 'off');
    xlabel(ax11, 'k-Wave X (voxel)');
    ylabel(ax11, 'k-Wave Z (voxel)');
    title(ax11, 'True  p_0  (Pa)');

    % ---- (1,2)  Spherical reconstruction --------------------------------
    ax12 = subplot(3, 3, 2);
    imagesc(ax12, 1:Nx, 1:Nz, p0s_xz');
    colormap(ax12, hot_cmap);
    caxis(ax12, crange);  colorbar(ax12);
    set(ax12, 'YDir', 'normal');
    hold(ax12, 'on');
    % Draw sphere as scattered red dots in this slice
    if ~isempty(sph_ix)
        scatter(ax12, sph_ix, sph_iz, 4, [1 0.15 0.15], 'filled', ...
            'MarkerFaceAlpha', 0.6);
    end
    hold(ax12, 'off');
    xlabel(ax12, 'k-Wave X (voxel)');
    ylabel(ax12, 'k-Wave Z (voxel)');
    title(ax12, 'Spherical sensor — reconstruction  (Pa)');

    % ---- (1,3)  Planar reconstruction -----------------------------------
    ax13 = subplot(3, 3, 3);
    imagesc(ax13, 1:Nx, 1:Nz, p0p_xz');
    colormap(ax13, hot_cmap);
    caxis(ax13, crange);  colorbar(ax13);
    set(ax13, 'YDir', 'normal');
    hold(ax13, 'on');
    % Horizontal red line at Z = pz_sensor
    plot(ax13, [0.5, Nx + 0.5], [pz_sensor, pz_sensor], 'r-', 'LineWidth', 2.5);
    hold(ax13, 'off');
    xlabel(ax13, 'k-Wave X (voxel)');
    ylabel(ax13, 'k-Wave Z (voxel)');
    title(ax13, sprintf('Planar sensor (z=%d) — reconstruction  (Pa)', pz_sensor));

    % ---- (2,1)  Medium schematic + true p0 contours ---------------------
    ax21 = subplot(3, 3, 4);
    med_rgb = build_medium_rgb(lbl_xz, n_tissues, props_list, Nx, Nz);
    % image() expects [Nz x Nx x 3]
    image(ax21, 1:Nx, 1:Nz, permute(med_rgb, [2, 1, 3]));
    set(ax21, 'YDir', 'normal');
    hold(ax21, 'on');
    % White p0 iso-contours at 10–90% of max
    if cmax > eps
        lvls = linspace(cmax * 0.10, cmax * 0.90, 5);
        contour(ax21, 1:Nx, 1:Nz, p0t_xz', lvls, 'w', 'LineWidth', 0.9);
    end
    % Red horizontal line = planar sensor
    plot(ax21, [0.5, Nx + 0.5], [pz_sensor, pz_sensor], 'r-', 'LineWidth', 2.5);
    % Cyan cross = dose centroid
    plot(ax21, cx, cz, 'c+', 'MarkerSize', 12, 'LineWidth', 2);
    hold(ax21, 'off');
    xlabel(ax21, 'k-Wave X (voxel)');
    ylabel(ax21, 'k-Wave Z (voxel)');
    title(ax21, 'Medium layout  (white = p_0 contours)');
    % Legend
    add_tissue_legend(ax21, props_list, n_tissues, pz_sensor);

    % ---- (2,2)  Spherical gamma map ------------------------------------
    ax22 = subplot(3, 3, 5);
    gs_display = gs_xz;
    gs_display(~gmask_xz) = NaN;
    imagesc(ax22, 1:Nx, 1:Nz, gs_display');
    colormap(ax22, gam_cmap);
    caxis(ax22, [0, 2]);
    cb = colorbar(ax22);
    cb.Label.String = '\gamma';
    set(ax22, 'YDir', 'normal', 'Color', [0.15 0.15 0.15]);
    hold(ax22, 'on');
    if ~isempty(sph_ix)
        scatter(ax22, sph_ix, sph_iz, 4, [1 0.15 0.15], 'filled', ...
            'MarkerFaceAlpha', 0.55);
    end
    hold(ax22, 'off');
    xlabel(ax22, 'k-Wave X (voxel)');
    ylabel(ax22, 'k-Wave Z (voxel)');
    title(ax22, sprintf('Spherical  \\gamma  (pass rate %.1f%%)', pr_sph));

    % ---- (2,3)  Planar gamma map ----------------------------------------
    ax23 = subplot(3, 3, 6);
    gp_display = gp_xz;
    gp_display(~gmask_xz) = NaN;
    imagesc(ax23, 1:Nx, 1:Nz, gp_display');
    colormap(ax23, gam_cmap);
    caxis(ax23, [0, 2]);
    cb = colorbar(ax23);
    cb.Label.String = '\gamma';
    set(ax23, 'YDir', 'normal', 'Color', [0.15 0.15 0.15]);
    hold(ax23, 'on');
    plot(ax23, [0.5, Nx + 0.5], [pz_sensor, pz_sensor], 'r-', 'LineWidth', 2.5);
    hold(ax23, 'off');
    xlabel(ax23, 'k-Wave X (voxel)');
    ylabel(ax23, 'k-Wave Z (voxel)');
    title(ax23, sprintf('Planar  \\gamma  (pass rate %.1f%%)', pr_plan));

    % ---- (3, 1-3)  Convergence — spans full bottom row -----------------
    ax31 = subplot(3, 3, [7 8 9]);
    hold(ax31, 'on');
    if ~isempty(conv_sph)
        semilogy(ax31, 2:numel(conv_sph)+1, conv_sph, 'b-o', ...
            'LineWidth', 1.8, 'MarkerSize', 5, ...
            'DisplayName', sprintf('Spherical (converged=%d)', ...
                find_converged_iter(conv_sph, CONFIG.conv_tol)));
    else
        semilogy(ax31, NaN, NaN, 'b-o', 'LineWidth', 1.8, 'DisplayName', 'Spherical (<2 iter)');
    end
    if ~isempty(conv_plan)
        semilogy(ax31, 2:numel(conv_plan)+1, conv_plan, 'r-s', ...
            'LineWidth', 1.8, 'MarkerSize', 5, ...
            'DisplayName', sprintf('Planar (converged=%d)', ...
                find_converged_iter(conv_plan, CONFIG.conv_tol)));
    else
        semilogy(ax31, NaN, NaN, 'r-s', 'LineWidth', 1.8, 'DisplayName', 'Planar (<2 iter)');
    end
    yline(ax31, CONFIG.conv_tol, 'k--', '1% threshold', ...
        'LineWidth', 1.2, 'LabelHorizontalAlignment', 'right', ...
        'LabelVerticalAlignment', 'top');
    hold(ax31, 'off');
    n_max = max([numel(conv_sph), numel(conv_plan), 1]);
    xlim(ax31, [1.5, n_max + 1.5]);
    xlabel(ax31, 'TR Iteration');
    ylabel(ax31, 'Relative Change  ||p_n - p_{n-1}|| / ||p_{n-1}||');
    title(ax31, 'Iterative TR Convergence History');
    legend(ax31, 'Location', 'best');
    grid(ax31, 'on');

    drawnow;
end


%% =========================================================================
%  HELPER: build RGB tissue image
%% =========================================================================

function med_rgb = build_medium_rgb(label_xz, n_tissues, props_list, Nx, Nz)
    med_rgb = zeros(Nx, Nz, 3);
    for r = 1:n_tissues
        tmask = (label_xz == r);
        c = props_list{r}.color;
        for ch = 1:3
            layer = med_rgb(:, :, ch);
            layer(tmask) = c(ch);
            med_rgb(:, :, ch) = layer;
        end
    end
end


%% =========================================================================
%  HELPER: add tissue legend patches to axes
%% =========================================================================

function add_tissue_legend(ax, props_list, n_tissues, pz_sensor)
    hold(ax, 'on');
    h = gobjects(n_tissues + 1, 1);
    for r = 1:n_tissues
        h(r) = patch(ax, NaN, NaN, props_list{r}.color, 'EdgeColor', 'k', ...
            'DisplayName', props_list{r}.name);
    end
    h(end) = plot(ax, NaN, NaN, 'r-', 'LineWidth', 2, ...
        'DisplayName', sprintf('Planar sensor (z=%d)', pz_sensor));
    hold(ax, 'off');
    legend(ax, h, 'Location', 'southeast', 'FontSize', 8);
end


%% =========================================================================
%  HELPER: find first convergence iteration
%% =========================================================================

function it_conv = find_converged_iter(conv_hist, tol)
    idx = find(conv_hist <= tol, 1);
    if isempty(idx)
        it_conv = numel(conv_hist) + 1;  % did not converge within iterations
    else
        it_conv = idx + 1;  % conv_hist starts at iteration 2
    end
end


%% =========================================================================
%  HELPER: green-yellow-red gamma colormap
%% =========================================================================

function cmap = gamma_colormap()
    n    = 256;
    half = round(n / 2);
    rest = n - half;
    r = [linspace(0,   1, half)'; ones(rest, 1)];
    g = [ones(half, 1) * 0.8;    linspace(0.8, 0, rest)'];
    b = [linspace(0.2, 0, half)'; zeros(rest, 1)];
    cmap = [r, g, b];
end
