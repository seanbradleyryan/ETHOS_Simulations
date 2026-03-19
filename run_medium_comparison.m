%% RUN_MEDIUM_COMPARISON.m
%  Photoacoustic reconstruction comparison across 4 tissue media.
%  Each medium is simulated with both a spherical and planar sensor,
%  reconstructed via iterative time-reversal (up to 10 iter, 1% tol).
%
%  Media:
%    (a) Homogeneous oil
%    (b) Half oil / half water  — planar sensor on oil face (x = 2)
%    (c) Half oil / half water  — planar sensor on water face (x = Nx-1)
%    (d) Third oil / third bone / third water — planar sensor on oil face
%
%  Per-medium figure layout (2 rows x 2 columns):
%    [1,1] Spherical reconstruction (XZ slice) with sensor outline in red
%    [1,2] Planar reconstruction (XZ slice) with sensor line in red
%    [2,1] Convergence history: spherical vs planar
%    [2,2] Medium schematic with dose contours (white) overlaid
%
%  Loaded from: RayStationFiles/<patient>/<session>/processed/
%    total_rs_dose.mat   -- used as dose distribution (p0 source)
%    sct_resampled.mat   -- used only for grid spacing
%
%  Author: ETHOS Pipeline Team
%  Date:   March 2026
%  See also: run_standalone_simulation.m, create_acoustic_medium.m

clear; clc; close all;

fprintf('=====================================================\n');
fprintf('  Multi-Medium Photoacoustic Reconstruction Study\n');
fprintf('=====================================================\n\n');

%% ========================= CONFIGURATION ================================

CONFIG.working_dir        = '/mnt/weka/home/80030361/ETHOS_Simulations';
CONFIG.patient_id         = '1194203';
CONFIG.session            = 'Session_1';
CONFIG.dose_filename      = 'total_rs_dose.mat';

CONFIG.downscale_factor   = 2;       % Spatial downscale for speed (1 = off)
CONFIG.meterset           = 140;     % Monitor units
CONFIG.dose_per_pulse_cGy = 0.16;   % cGy per LINAC pulse
CONFIG.pml_size           = 10;     % PML thickness in voxels
CONFIG.cfl_number         = 0.3;    % CFL stability number
CONFIG.use_gpu            = true;   % GPU acceleration
CONFIG.max_tr_iter        = 10;     % Max iterative TR iterations
CONFIG.conv_tol           = 0.01;   % 1% relative-change convergence threshold

% Acoustic tissue properties: density (kg/m^3), sound_speed (m/s), Gruneisen, RGB colour
% Bone Gruneisen is set to 0.5 (non-zero) so bone regions produce visible signal.
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

sct_loaded  = load(fullfile(processed_dir, 'sct_resampled.mat'));
spacing_mm  = double(sct_loaded.sct_resampled.spacing(:)');

fprintf('[Load] Dose: [%d x %d x %d], max=%.4f Gy   Spacing: [%.2f %.2f %.2f] mm\n', ...
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
dx = spacing_mm(1)/1000;
dy = spacing_mm(2)/1000;
dz = spacing_mm(3)/1000;
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
dose_pp    = doseGrid / num_pulses;   % dose per pulse (Gy)

%% ========================= MEDIA SPECIFICATIONS =========================

split_h  = floor(Nx / 2);         % half-grid boundary
split_t1 = floor(Nx / 3);         % first-third boundary
split_t2 = floor(2 * Nx / 3);     % second-third boundary

%  Cell array rows:
%    col 1  title string
%    col 2  X-region index ranges  (cell of vectors)
%    col 3  props structs          (cell matching col 2)
%    col 4  planar sensor X index
%    col 5  planar sensor label
media_specs = { ...
    '(a) Homogeneous Oil', ...
        {1:Nx}, ...
        {PROPS.oil}, ...
        2, 'x=2 (oil)'; ...
    ...
    '(b) Half Oil / Half Water  |  Sensor: oil face', ...
        {1:split_h, split_h+1:Nx}, ...
        {PROPS.oil, PROPS.water}, ...
        2, 'x=2 (oil)'; ...
    ...
    '(c) Half Oil / Half Water  |  Sensor: water face', ...
        {1:split_h, split_h+1:Nx}, ...
        {PROPS.oil, PROPS.water}, ...
        Nx-1, sprintf('x=%d (water)', Nx-1); ...
    ...
    '(d) Oil / Bone / Water  |  Sensor: oil face', ...
        {1:split_t1, split_t1+1:split_t2, split_t2+1:Nx}, ...
        {PROPS.oil, PROPS.bone, PROPS.water}, ...
        2, 'x=2 (oil)'; ...
};

%% ========================= MAIN SIMULATION LOOP =========================

for m = 1:4
    m_title  = media_specs{m, 1};
    m_regs   = media_specs{m, 2};
    m_props  = media_specs{m, 3};
    m_px     = media_specs{m, 4};
    m_plabel = media_specs{m, 5};

    fprintf('\n+-------------------------------------------------------+\n');
    fprintf('|  MEDIUM %d: %s\n', m, m_title);
    fprintf('+-------------------------------------------------------+\n');

    %% Build spatially-varying medium arrays
    density_a = zeros(gridSize);
    speed_a   = zeros(gridSize);
    gruen_a   = zeros(gridSize);
    label_a   = zeros(gridSize, 'uint8');

    for r = 1:numel(m_regs)
        xi_r = m_regs{r};
        p    = m_props{r};
        density_a(xi_r, :, :) = p.density;
        speed_a(xi_r, :, :)   = p.sound_speed;
        gruen_a(xi_r, :, :)   = p.gruneisen;
        label_a(xi_r, :, :)   = r;
    end

    kmedium = struct('density',     density_a, ...
                     'sound_speed', speed_a, ...
                     'alpha_coeff', 0, ...
                     'alpha_power', 1.1);

    %% Initial pressure:  p0 = (dose/num_pulses) * Gruneisen * rho
    p0 = double(dose_pp .* gruen_a .* density_a);
    p0 = double(smooth(p0));   % k-Wave Hann-window smoothing
    p0 = max(p0, 0);
    fprintf('  p0 range: [%.2e, %.2e] Pa\n', min(p0(:)), max(p0(:)));

    %% CFL-stable time step for this medium's sound speeds
    maxC    = max(speed_a(:));
    minC    = min(speed_a(speed_a > 0));
    dt      = CONFIG.cfl_number * min([dx, dy, dz]) / maxC;
    diagLen = sqrt((Nx*dx)^2 + (Ny*dy)^2 + (Nz*dz)^2);
    Nt      = ceil(2.5 * diagLen / minC / dt);

    kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
    kgrid.dt = dt;
    kgrid.Nt = Nt;
    fprintf('  dt = %.2e s,  Nt = %d\n', dt, Nt);

    %% Build sensor masks
    sph_r = max(1, floor(min(gridSize) / 2) - CONFIG.pml_size - 2);
    sensor_sph.mask = makeSphere(Nx, Ny, Nz, sph_r);

    sensor_plan.mask = false(gridSize);
    sensor_plan.mask(m_px, :, :) = true;

    fprintf('  Spherical: r=%d, %d voxels  |  Planar: x=%d, %d voxels\n', ...
        sph_r, sum(sensor_sph.mask(:)), m_px, sum(sensor_plan.mask(:)));

    source_fwd.p0 = p0;

    % ==================================================================
    % SPHERICAL SENSOR: forward simulation + iterative TR
    % ==================================================================
    fprintf('\n  [Spherical] Forward simulation...\n');
    t0 = tic;
    sd_sph = kspaceFirstOrder3D(kgrid, kmedium, source_fwd, sensor_sph, inputArgs{:});
    fprintf('  [Spherical] Forward: %.1f s.  Starting iterative TR...\n', toc(t0));

    [p0r_sph, conv_sph] = run_iterative_tr(kgrid, kmedium, sensor_sph, sd_sph, ...
        gridSize, CONFIG.max_tr_iter, CONFIG.conv_tol, inputArgs);

    p0r_sph = double(gather(p0r_sph));
    clear sd_sph

    % ==================================================================
    % PLANAR SENSOR: forward simulation + iterative TR
    % ==================================================================
    fprintf('\n  [Planar / %s] Forward simulation...\n', m_plabel);
    t0 = tic;
    sd_plan = kspaceFirstOrder3D(kgrid, kmedium, source_fwd, sensor_plan, inputArgs{:});
    fprintf('  [Planar] Forward: %.1f s.  Starting iterative TR...\n', toc(t0));

    [p0r_plan, conv_plan] = run_iterative_tr(kgrid, kmedium, sensor_plan, sd_plan, ...
        gridSize, CONFIG.max_tr_iter, CONFIG.conv_tol, inputArgs);

    p0r_plan = double(gather(p0r_plan));
    clear sd_plan

    %% Generate per-medium comparison figure
    make_comparison_figure(m, m_title, ...
        double(gather(p0)), p0r_sph, p0r_plan, ...
        sensor_sph.mask, sensor_plan.mask, ...
        conv_sph, conv_plan, ...
        label_a, m_props, ...
        cy, cx, cz, m_px, m_plabel, gridSize, CONFIG.conv_tol);

    clear density_a speed_a gruen_a label_a p0 p0r_sph p0r_plan kmedium
    fprintf('\n  Medium %d complete.\n', m);
end

fprintf('\n=====================================================\n');
fprintf('  All simulations complete.\n');
fprintf('=====================================================\n');


%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

function [p0_out, conv_hist] = run_iterative_tr(kgrid, kmedium, sensor, ...
    sd_meas, gridSize, max_iter, conv_tol, in_args)
%RUN_ITERATIVE_TR  Iterative Dirichlet time-reversal reconstruction.
%
%   Each iteration:
%     1. TR(current_data)  -> p0_estimate  (positivity constrained)
%     2. Forward(p0_estimate) -> sd_recon
%     3. Residual update:  current_data = measured + (measured - sd_recon)
%     4. Convergence check: ||p0_n - p0_{n-1}|| / ||p0_{n-1}|| < conv_tol
%
%   Returns:
%     p0_out    - final reconstruction (double, CPU)
%     conv_hist - relative change values, starting at iteration 2

    Nx = gridSize(1);  Ny = gridSize(2);  Nz = gridSize(3);
    sd_work = sd_meas;
    p0_prev = zeros(gridSize, 'single');
    conv_hist = [];
    p0_out    = zeros(gridSize);

    for it = 1:max_iter

        % --- Dirichlet time-reversal ---
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

        % --- Convergence check (from iteration 2 onward) ---
        if it > 1
            nn = norm(p0_prev(:));
            rc = double(norm(p0_it(:) - p0_prev(:))) / max(double(nn), 1e-30);
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

        % --- Residual correction for next iteration ---
        if it < max_iter
            src_res.p0 = p0_it;
            sd_res     = kspaceFirstOrder3D(kgrid, kmedium, src_res, sensor, in_args{:});
            sd_work    = sd_meas + (sd_meas - sd_res);
        end
    end
end


function make_comparison_figure(m_idx, m_title, p0_true, p0_sph, p0_plan, ...
    mask_sph, mask_plan, conv_sph, conv_plan, ...
    label_arr, props_list, cy, cx, cz, px_sensor, pl_label, gridSize, conv_tol)
%MAKE_COMPARISON_FIGURE  2x2 per-medium comparison figure.
%
%   All views show XZ slices (k-Wave dim-1 vs dim-3) at Y=cy.
%   imagesc is called as imagesc(1:Nx, 1:Nz, data') so:
%     horizontal axis = k-Wave X  (shows medium layering)
%     vertical   axis = k-Wave Z

    Nx = gridSize(1);
    Nz = gridSize(3);
    n_tissues = numel(props_list);

    % Extract XZ slices at Y=cy  =>  [Nx x Nz] matrices
    p0_true_xz  = squeeze(p0_true(:, cy, :));
    p0_sph_xz   = squeeze(p0_sph(:,  cy, :));
    p0_plan_xz  = squeeze(p0_plan(:, cy, :));
    label_xz    = squeeze(label_arr(:, cy, :));
    msph_xz     = squeeze(mask_sph(:, cy, :));  % logical [Nx x Nz]

    % Common colour scale
    cmax = max(p0_true_xz(:));
    if cmax < eps, cmax = 1; end
    crange = [0, cmax * 1.02];

    % Spherical sensor positions visible in this XZ slice
    [sph_ix, sph_iz] = find(msph_xz);

    % Create figure
    figure('Name', sprintf('Medium %d', m_idx), 'NumberTitle', 'off', ...
        'Position', [30 + 30*m_idx, 30 + 20*m_idx, 1500, 760], 'Color', 'w');
    sgtitle(sprintf('Medium %d: %s', m_idx, m_title), ...
        'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'none');

    hot_cmap = hot(256);

    % ------------------------------------------------------------------
    % Panel (1,1):  Spherical reconstruction with sensor outline
    % ------------------------------------------------------------------
    ax1 = subplot(2, 2, 1);
    imagesc(ax1, 1:Nx, 1:Nz, p0_sph_xz');
    colormap(ax1, hot_cmap);
    caxis(ax1, crange);
    colorbar(ax1, 'Location', 'eastoutside');
    set(ax1, 'YDir', 'normal');
    hold(ax1, 'on');
    % Scatter spherical sensor voxels as small red dots
    if ~isempty(sph_ix)
        scatter(ax1, sph_ix, sph_iz, 3, [1 0.15 0.15], 'filled', ...
            'MarkerFaceAlpha', 0.55);
    end
    hold(ax1, 'off');
    xlabel(ax1, 'k-Wave X (voxel)');
    ylabel(ax1, 'k-Wave Z (voxel)');
    title(ax1, 'Spherical sensor — reconstruction');

    % ------------------------------------------------------------------
    % Panel (1,2):  Planar reconstruction with sensor line
    % ------------------------------------------------------------------
    ax2 = subplot(2, 2, 2);
    imagesc(ax2, 1:Nx, 1:Nz, p0_plan_xz');
    colormap(ax2, hot_cmap);
    caxis(ax2, crange);
    colorbar(ax2, 'Location', 'eastoutside');
    set(ax2, 'YDir', 'normal');
    hold(ax2, 'on');
    % Vertical red line at planar sensor X position
    plot(ax2, [px_sensor px_sensor], [0.5, Nz + 0.5], 'r-', 'LineWidth', 2.5);
    hold(ax2, 'off');
    xlabel(ax2, 'k-Wave X (voxel)');
    ylabel(ax2, 'k-Wave Z (voxel)');
    title(ax2, sprintf('Planar sensor (%s) — reconstruction', pl_label));

    % ------------------------------------------------------------------
    % Panel (2,1):  Convergence history
    % ------------------------------------------------------------------
    ax3 = subplot(2, 2, 3);
    hold(ax3, 'on');
    if ~isempty(conv_sph)
        iters_sph = 2 : numel(conv_sph) + 1;
        semilogy(ax3, iters_sph, conv_sph, 'b-o', ...
            'LineWidth', 1.8, 'MarkerSize', 5, 'DisplayName', 'Spherical');
    else
        semilogy(ax3, NaN, NaN, 'b-o', 'LineWidth', 1.8, ...
            'DisplayName', 'Spherical (<2 iter)');
    end
    if ~isempty(conv_plan)
        iters_plan = 2 : numel(conv_plan) + 1;
        semilogy(ax3, iters_plan, conv_plan, 'r-s', ...
            'LineWidth', 1.8, 'MarkerSize', 5, 'DisplayName', 'Planar');
    else
        semilogy(ax3, NaN, NaN, 'r-s', 'LineWidth', 1.8, ...
            'DisplayName', 'Planar (<2 iter)');
    end
    yline(ax3, conv_tol, 'k--', ...
        sprintf('1%%  threshold  (%g)', conv_tol), ...
        'LineWidth', 1.0, 'LabelHorizontalAlignment', 'right');
    hold(ax3, 'off');
    n_max = max([numel(conv_sph), numel(conv_plan), 1]);
    xlim(ax3, [1.5, n_max + 1.5]);
    xlabel(ax3, 'TR Iteration');
    ylabel(ax3, 'Relative Change');
    title(ax3, 'Convergence history');
    legend(ax3, 'Location', 'best');
    grid(ax3, 'on');

    % ------------------------------------------------------------------
    % Panel (2,2):  Medium schematic + true p0 dose contours
    % ------------------------------------------------------------------
    ax4 = subplot(2, 2, 4);

    % Build [Nx x Nz x 3] RGB image from tissue labels
    med_rgb = zeros(Nx, Nz, 3);
    for r = 1:n_tissues
        tissue_mask = (label_xz == r);
        c = props_list{r}.color;
        for ch = 1:3
            layer = med_rgb(:, :, ch);
            layer(tissue_mask) = c(ch);
            med_rgb(:, :, ch) = layer;
        end
    end

    % image() needs [Nz x Nx x 3] when XData=1:Nx, YData=1:Nz
    image(ax4, 1:Nx, 1:Nz, permute(med_rgb, [2, 1, 3]));
    set(ax4, 'YDir', 'normal');
    hold(ax4, 'on');

    % Overlay true p0 as white contours
    if cmax > eps
        lvls = linspace(cmax * 0.10, cmax * 0.90, 5);
        contour(ax4, 1:Nx, 1:Nz, p0_true_xz', lvls, 'w', 'LineWidth', 0.9);
    end

    % Red line for planar sensor position
    plot(ax4, [px_sensor, px_sensor], [0.5, Nz + 0.5], 'r-', 'LineWidth', 2.5);

    % Red cross at dose centroid (projected onto XZ)
    plot(ax4, cx, cz, 'r+', 'MarkerSize', 12, 'LineWidth', 2);

    hold(ax4, 'off');
    xlabel(ax4, 'k-Wave X (voxel)');
    ylabel(ax4, 'k-Wave Z (voxel)');
    title(ax4, 'Medium layout  (white = p_0 contours,  red = sensor + centroid)');

    % Tissue legend
    hold(ax4, 'on');
    h = gobjects(n_tissues, 1);
    for r = 1:n_tissues
        h(r) = patch(ax4, NaN, NaN, props_list{r}.color, 'EdgeColor', 'k', ...
            'DisplayName', props_list{r}.name);
    end
    hold(ax4, 'off');
    legend(ax4, h, 'Location', 'southeast', 'FontSize', 8);

    drawnow;
end
