function plot_dose_panels(original, recon, sensor_mask, density, spacing_mm, titleStr, smooth_sigma, rowLabels)
%PLOT_DOSE_PANELS 2x3 dose comparison: coronal, sagittal, axial.
%  Row 1 = first volume,  Row 2 = second volume (labelled via rowLabels;
%  default {'Original','Reconstructed'}).
%  Dose (jet, semi-transparent) is overlaid on the density map (grayscale).
%  Voxels with dose < 10% of max are fully transparent (masked out).
%  Both rows share an identical dose colour range [0, max(original)]
%  so magnitudes are directly comparable.
%  Isocenter at max-dose voxel.  Sensor contour in red on every panel.
%  smooth_sigma (voxels) Gaussian-smooths each slice for display only, to
%  fill speckle gaps in a spotty recon. Pass 0 to disable.
%
%  Pass density=[] to show dose only with a black background.
%  (Shared copy; formerly a local in the standalone drivers.)

    if nargin < 7 || isempty(smooth_sigma), smooth_sigma = 0; end
    if nargin < 8 || isempty(rowLabels),    rowLabels = {'Original', 'Reconstructed'}; end

    gridSize = size(original);
    if ~isequal(size(sensor_mask), gridSize)
        sensor_mask = sensor_mask(1:gridSize(1), 1:gridSize(2), 1:gridSize(3));
    end

    have_density = ~isempty(density) && isequal(size(density), gridSize);

    [~, max_idx] = max(original(:));
    [cx, cy, cz] = ind2sub(gridSize, max_idx);

    % Shared colour scale anchored to original (reference) dose
    max_dose = max(original(:));
    if max_dose == 0, max_dose = 1; end

    x_ax = (1:gridSize(1)) * spacing_mm(1);
    y_ax = (1:gridSize(2)) * spacing_mm(2);
    z_ax = (1:gridSize(3)) * spacing_mm(3);

    % Jet LUT for manual RGB blending (avoids per-axes colormap conflict)
    cmap_jet = jet(256);

    % Soft-tissue window/level for the density background (kg/m^3).
    wl_center = 1050;          % kg/m^3  (soft tissue ~1000-1080)
    wl_width  = 350;           % kg/m^3
    wl_min    = wl_center - wl_width / 2;   % 875  kg/m^3
    wl_max    = wl_center + wl_width / 2;   % 1225 kg/m^3 %#ok<NASGU>

    figure('Name', titleStr, 'Color', 'w', 'NumberTitle', 'off', ...
        'Position', [50, 50, 1380, 700]);
    sgtitle(sprintf('%s\nIsocenter (max dose): X=%d  Y=%d  Z=%d voxel  |  Dose clim [0, %.4f] Gy', ...
        titleStr, cx, cy, cz, max_dose), 'FontWeight', 'bold', 'FontSize', 11);

    row_labels = rowLabels;
    doses      = {original, recon};

    for row = 1:2
        d   = doses{row};
        lbl = row_labels{row};

        % Slice data for the three views
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

            % --- Background: density as grayscale RGB ---
            if have_density
                dn     = (density_slices{col} - wl_min) / wl_width;
                dn     = max(0, min(1, dn));   % clip to [0,1]
                bg_rgb = repmat(dn, [1, 1, 3]);
            else
                bg_rgb = zeros([size(dose_slice), 3]);   % black
            end
            image(ax, xv, yv, bg_rgb);
            hold(ax, 'on');

            % --- Foreground: dose as jet RGB with alpha mask ---
            %   mask: dose >= 10% of max_dose  (shared threshold)
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

            % --- Sensor contour ---
            s = sensor_slices{col};
            if any(s(:))
                contour(ax, xv, yv, s, [0.5, 0.5], 'r-', 'LineWidth', 2);
            end
            hold(ax, 'off');

            axis(ax, 'xy'); axis(ax, 'image');

            % Colorbar: attach jet LUT with shared dose clim
            colormap(ax, cmap_jet);
            caxis(ax, [0, max_dose]);
            cb = colorbar(ax); cb.Label.String = 'Dose (Gy)';

            xlabel(ax, xlbls{col}); ylabel(ax, ylbls{col});
            title(ax, sprintf('%s  %s', lbl, tsuffs{col}));
        end
    end
    drawnow;
end
