function plot_gamma_and_error_axial(gamma_results, original, recon, sensor_mask, cz, titleStr, parent)
%PLOT_GAMMA_AND_ERROR_AXIAL 1xN axial figure:
%  one gamma panel per criterion + truth dose + recon dose + absolute error.
%  The truth, recon, and |error| panels share the same 'hot' colour setting;
%  truth and recon additionally share a common dose colour range so their
%  magnitudes are directly comparable.
%  Sensor contour in red.  Slice at cz (max-dose axial index).
%  Optional titleStr is prepended to the figure super-title.
%  Optional parent (figure / uitab / uipanel) hosts the plot; when omitted a
%  new figure is created. Passing a uitab lets the caller collect several of
%  these panels as tabs of one window.
%  (Shared copy; formerly a local in the standalone drivers.)

    if nargin < 6, titleStr = ''; end
    if nargin < 7, parent = []; end

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
    nCol       = nCrit + 3;   % gamma panel(s) + truth + recon + |error|

    if isempty(parent)
        parent = figure('Name', 'Gamma, Dose & Absolute Error  Axial', 'Color', 'w', ...
            'NumberTitle', 'off', 'Position', [50, 300, 1400, 370]);
    end
    t = tiledlayout(parent, 1, nCol, 'TileSpacing', 'compact', 'Padding', 'compact');
    if isempty(titleStr)
        supt = sprintf('Axial Plane (Z = %d voxel)  Gamma Index, Dose & Absolute Error', cz);
    else
        supt = sprintf('%s\nAxial Plane (Z = %d voxel)  Gamma Index, Dose & Absolute Error', titleStr, cz);
    end
    title(t, supt, 'FontWeight', 'bold', 'FontSize', 11);

    gamma_clim   = [0, 2];
    sensor_slice = squeeze(sensor_mask(:, :, cz))';

    % Gamma panels
    for g = 1:nCrit
        ax = nexttile(t, g);
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

    % Truth and recon dose slices, shared 'hot' clim so they are comparable.
    orig_slice  = squeeze(original(:, :, cz))';
    recon_slice = squeeze(recon(:, :, cz))';
    dose_max    = max([orig_slice(:); recon_slice(:)]);
    if dose_max == 0, dose_max = 1; end

    % Truth (reference) dose panel
    ax = nexttile(t, nCrit + 1);
    imagesc(ax, orig_slice, [0, dose_max]);
    axis(ax, 'xy'); axis(ax, 'image');
    colormap(ax, 'hot');
    cb = colorbar(ax); cb.Label.String = 'Dose (Gy)';
    hold(ax, 'on');
    if any(sensor_slice(:))
        contour(ax, sensor_slice, [0.5, 0.5], 'r-', 'LineWidth', 2);
    end
    hold(ax, 'off');
    xlabel(ax, 'X (voxel)'); ylabel(ax, 'Y (voxel)');
    title(ax, sprintf('Truth\nMax: %.4f Gy', max(orig_slice(:))));

    % Recon dose panel (same colour setting / clim as truth)
    ax = nexttile(t, nCrit + 2);
    imagesc(ax, recon_slice, [0, dose_max]);
    axis(ax, 'xy'); axis(ax, 'image');
    colormap(ax, 'hot');
    cb = colorbar(ax); cb.Label.String = 'Dose (Gy)';
    hold(ax, 'on');
    if any(sensor_slice(:))
        contour(ax, sensor_slice, [0.5, 0.5], 'r-', 'LineWidth', 2);
    end
    hold(ax, 'off');
    xlabel(ax, 'X (voxel)'); ylabel(ax, 'Y (voxel)');
    title(ax, sprintf('Recon\nMax: %.4f Gy', max(recon_slice(:))));

    % Absolute error panel (same 'hot' colour setting)
    ax = nexttile(t, nCrit + 3);
    err_slice = abs(recon_slice - orig_slice);
    max_err   = max(err_slice(:));
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
    title(ax, sprintf('|Recon - Truth|\nMax: %.4f Gy', max_err));

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
