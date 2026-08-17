%% =========================================================================
%  PLOT_STANDALONE_RESULTS.m
%  Reproduce all post-simulation figures from a saved standalone_results_*.mat
%  file without re-running any k-Wave simulations.
%
%  Usage:
%    plot_standalone_results                  % prompts for file via GUI
%    plot_standalone_results('my_results.mat')
%
%  Required fields in the .mat file (saved by run_standalone_simulation.m):
%    recon_dose, original_dose, p0, reconPressure
%    sensor_mask, spacing_mm, grid_size, config
%
%  Optional fields (gracefully skipped if absent):
%    gamma            – struct with maps/pass_rates/criteria/eval_mask
%    conv_max_pressure, conv_rel_change, num_iters_done  – convergence history
%    density          – medium density volume for CT-style background
%    fwd_time_sec, tr_time_sec
%  =========================================================================

function plot_standalone_results(results_file)

close all;

%% ---- Resolve file path -----------------------------------------------
if nargin < 1 || isempty(results_file)
    [fname, fpath] = uigetfile('*.mat', 'Select standalone results .mat file');
    if isequal(fname, 0)
        fprintf('No file selected. Exiting.\n');
        return;
    end
    results_file = fullfile(fpath, fname);
end

if ~isfile(results_file)
    error('File not found: %s', results_file);
end

fprintf('=========================================================\n');
fprintf('  plot_standalone_results\n');
fprintf('  Loading: %s\n', results_file);
fprintf('=========================================================\n\n');

%% ---- Load results -------------------------------------------------------
R = load(results_file);

% Mandatory fields
required = {'recon_dose','original_dose','p0','sensor_mask','spacing_mm','grid_size'};
for k = 1:numel(required)
    if ~isfield(R, required{k})
        error('Results file is missing required field: %s', required{k});
    end
end

recon_dose    = double(R.recon_dose);
original_dose = double(R.original_dose);
p0            = double(R.p0);
sensor_mask   = logical(R.sensor_mask);
spacing_mm    = R.spacing_mm(:)';
gridSize      = R.grid_size;

% Optional: CONFIG
if isfield(R, 'config')
    CONFIG = R.config;
    patient_id   = CONFIG.patient_id;
    sensor_label = CONFIG.sensor_placement_method;
    tol          = CONFIG.convergence_tol;
    max_iter_cfg = CONFIG.num_time_reversal_iter;
else
    patient_id   = 'unknown';
    sensor_label = 'unknown';
    tol          = 1e-3;
    max_iter_cfg = NaN;
end

% Optional: density for CT background in dose panels
if isfield(R, 'density')
    density = double(R.density);
else
    density = [];          % plot_dose_panels accepts [] -> black background
end

% Optional: convergence history
have_conv = isfield(R, 'conv_max_pressure') && isfield(R, 'conv_rel_change') ...
            && isfield(R, 'num_iters_done');
if have_conv
    conv_max_pressure = R.conv_max_pressure;
    conv_rel_change   = R.conv_rel_change;
    num_iters_done    = R.num_iters_done;
else
    fprintf('[WARN] Convergence history not found in results file.\n');
    fprintf('       Add conv_max_pressure / conv_rel_change / num_iters_done\n');
    fprintf('       to the save block in run_standalone_simulation.m to enable\n');
    fprintf('       the convergence plot.\n\n');
end

% Optional: gamma
have_gamma = isfield(R, 'gamma') && ~isempty(R.gamma) ...
             && isfield(R.gamma, 'maps');
if have_gamma
    gamma_results = R.gamma;
else
    gamma_results = [];
    fprintf('[INFO] No gamma results found in file – gamma figure skipped.\n\n');
end

% Optional: metadata for display
timing_str   = '';
if isfield(R, 'fwd_time_sec') && isfield(R, 'tr_time_sec')
    timing_str = sprintf('  Fwd: %.1f s   TR: %.1f s', R.fwd_time_sec, R.tr_time_sec);
end

%% ---- Print summary -------------------------------------------------------
fprintf('Patient:        %s\n', patient_id);
fprintf('Sensor:         %s\n', sensor_label);
fprintf('Grid size:      [%d x %d x %d]\n', gridSize(1), gridSize(2), gridSize(3));
fprintf('Spacing (mm):   [%.2f x %.2f x %.2f]\n', spacing_mm);
fprintf('Original dose:  [%.6f, %.4f] Gy\n', min(original_dose(:)), max(original_dose(:)));
fprintf('Recon dose:     [%.6f, %.4f] Gy\n', min(recon_dose(:)), max(recon_dose(:)));
if ~isempty(timing_str),   fprintf('%s\n', timing_str); end

doseThreshold = 0.01 * max(original_dose(:));
dose_region   = original_dose > doseThreshold;
if any(dose_region(:))
    abs_error = abs(recon_dose(dose_region) - original_dose(dose_region));
    rel_error = abs_error ./ max(original_dose(dose_region), 1e-10) * 100;
    fprintf('Mean abs err:   %.6f Gy\n', mean(abs_error));
    fprintf('Mean rel err:   %.2f%%\n',  mean(rel_error));
    fprintf('Max  rel err:   %.2f%%\n',  max(rel_error));
end
fprintf('\n');

%% ---- Figure 1: Sensor placement vs dose mask ----------------------------
dose_mask_vis = original_dose >= 0.10 * max(original_dose(:));
% Trim sensor_mask to original grid size in case of saved padded mask
sm_vis = sensor_mask;
if ~isequal(size(sm_vis), gridSize)
    sm_vis = sm_vis(1:gridSize(1), 1:gridSize(2), 1:gridSize(3));
end
plot_sensor_dose_planes(dose_mask_vis, sm_vis, spacing_mm, sensor_label);

%% ---- Figure 2: Dose comparison panels -----------------------------------
plot_dose_panels(original_dose, recon_dose, sm_vis, density, spacing_mm, ...
    'Dose Comparison: Original vs Reconstructed');

%% ---- Figure 3: Convergence history (only if data present) ---------------
if have_conv
    p0_max_for_plot = max(p0(:));
    plot_convergence_history(conv_max_pressure, conv_rel_change, ...
        num_iters_done, tol, p0_max_for_plot);
end

%% ---- Figure 4: Gamma index + absolute error (axial) --------------------
if have_gamma
    [~, max_dose_idx] = max(original_dose(:));
    [~, ~, cz_gamma]  = ind2sub(gridSize, max_dose_idx);
    plot_gamma_and_error_axial(gamma_results, original_dose, recon_dose, ...
        sm_vis, cz_gamma);
end

fprintf('All figures rendered.\n');
end


%% =========================================================================
%  LOCAL PLOT FUNCTIONS  (kept in sync with run_standalone_simulation.m)
%% =========================================================================

function plot_sensor_dose_planes(dose_mask, sensor_mask, spacing_mm, sensor_label)
%PLOT_SENSOR_DOSE_PLANES  1x3 anatomical view of sensor geometry vs dose mask.

    [Nx3, Ny3, Nz3] = size(dose_mask);

    x_ax = (1:Nx3) * spacing_mm(1);
    y_ax = (1:Ny3) * spacing_mm(2);
    z_ax = (1:Nz3) * spacing_mm(3);

    dose_cor = squeeze(any(dose_mask,   2));
    dose_sag = squeeze(any(dose_mask,   1));
    dose_axi = squeeze(any(dose_mask,   3));

    sens_cor = squeeze(any(sensor_mask, 2));
    sens_sag = squeeze(any(sensor_mask, 1));
    sens_axi = squeeze(any(sensor_mask, 3));

    figure('Name', 'Sensor Placement vs Dose Mask', 'Color', 'w', ...
        'NumberTitle', 'off', 'Position', [80, 80, 1300, 420]);
    sgtitle(sprintf('Sensor Placement vs Dose Mask  (\\geq10%% max)   |   Sensor: %s', ...
        sensor_label), 'FontWeight', 'bold', 'FontSize', 11);

    view_data = {
        dose_cor', sens_cor', x_ax, z_ax, 'X (mm)', 'Z (mm)', 'Coronal  (max-proj along Y)';
        dose_sag', sens_sag', y_ax, z_ax, 'Y (mm)', 'Z (mm)', 'Sagittal  (max-proj along X)';
        dose_axi', sens_axi', x_ax, y_ax, 'X (mm)', 'Y (mm)', 'Axial  (max-proj along Z)';
    };

    dose_color   = [0.20, 0.50, 0.90];
    sensor_color = [0.90, 0.10, 0.10];

    for col = 1:3
        ax  = subplot(1, 3, col);
        d2d = double(view_data{col, 1});
        s2d = double(view_data{col, 2});
        xv  = view_data{col, 3};
        yv  = view_data{col, 4};

        imagesc(ax, xv, yv, zeros(size(d2d)));
        colormap(ax, 'gray'); caxis(ax, [0, 1]);
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
        p2 = patch(ax, NaN, NaN, sensor_color,  'FaceAlpha', 0.85, 'EdgeColor', 'none');
        hold(ax, 'off');
        if col == 3
            legend(ax, [p1, p2], {'Dose mask (>=10%)', 'Sensor'}, ...
                'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 8);
        end
    end
    drawnow;
end


function plot_dose_panels(original, recon, sensor_mask, density, spacing_mm, titleStr)
%PLOT_DOSE_PANELS 2x3 dose comparison: coronal, sagittal, axial.
%  Row 1 = original dose,  Row 2 = reconstructed dose.
%  Pass density=[] for a black background.

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

    cmap_jet  = jet(256);
    wl_center = 1050;
    wl_width  = 350;
    wl_min    = wl_center - wl_width / 2;

    figure('Name', titleStr, 'Color', 'w', 'NumberTitle', 'off', ...
        'Position', [50, 50, 1380, 700]);
    sgtitle(sprintf('%s\nIsocenter (max dose): X=%d  Y=%d  Z=%d voxel  |  Dose clim [0, %.4f] Gy', ...
        titleStr, cx, cy, cz, max_dose), 'FontWeight', 'bold', 'FontSize', 11);

    row_labels = {'Original', 'Reconstructed'};
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
        tsuffs = { sprintf('Coronal (Y=%d)',  cy), ...
                   sprintf('Sagittal (X=%d)', cx), ...
                   sprintf('Axial (Z=%d)',    cz) };
        if have_density
            density_slices = { squeeze(density(:, cy, :))', ...
                               squeeze(density(cx, :, :))', ...
                               squeeze(density(:, :, cz))' };
        end

        for col = 1:3
            ax         = subplot(2, 3, (row-1)*3 + col);
            dose_slice = double(dose_slices{col});
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


function plot_convergence_history(conv_max_pressure, conv_rel_change, num_iters, tol, p0_max_orig)
%PLOT_CONVERGENCE_HISTORY p0 convergence over TR iterations.

    iters  = 1:num_iters;
    p_vals = conv_max_pressure(iters);
    rc_vals = conv_rel_change(iters);

    figure('Name', 'p0 Convergence', 'Color', 'w', ...
        'NumberTitle', 'off', 'Position', [150, 520, 720, 390]);

    yyaxis left;
    plot(iters, p_vals, 'b-o', 'LineWidth', 1.8, 'MarkerSize', 5, ...
        'MarkerFaceColor', [0.2, 0.4, 1.0]);
    hold on;
    if nargin >= 5 && ~isempty(p0_max_orig) && p0_max_orig > 0
        yline(p0_max_orig, 'k--', 'LineWidth', 1.8, ...
            'Label', sprintf('p_{0}^{orig} max = %.2e Pa', p0_max_orig), ...
            'LabelHorizontalAlignment', 'left', ...
            'LabelVerticalAlignment', 'bottom', 'FontSize', 8);
    end
    hold off;
    ylabel('Max Reconstructed p_0 (Pa)', 'Color', [0.2, 0.4, 1.0]);
    top_val = max([max(p_vals), p0_max_orig]) * 1.20;
    if top_val <= 0, top_val = 1; end
    ylim([0, top_val]);

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
%PLOT_GAMMA_AND_ERROR_AXIAL 1x4 axial figure: gamma 10/10 | 5/5 | 3/3 | abs error.

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

    figure('Name', 'Gamma & Absolute Error  Axial', 'Color', 'w', ...
        'NumberTitle', 'off', 'Position', [50, 300, 1400, 370]);
    sgtitle(sprintf('Axial Plane (Z = %d voxel)  Gamma Index & Absolute Error', cz), ...
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
    n    = 256;
    half = round(n / 2);
    rest = n - half;
    r = [linspace(0, 1, half)'; ones(rest, 1)];
    g = [linspace(0.85, 1, half)'; linspace(1, 0, rest)'];
    b = zeros(n, 1);
    cmap = [r, g, b];
end
