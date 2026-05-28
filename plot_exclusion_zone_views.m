function plot_exclusion_zone_views(sct, sensor_info, spacing_mm, fig_title)
%PLOT_EXCLUSION_ZONE_VIEWS Show three anatomical views of the beam exclusion zone
%
%   plot_exclusion_zone_views(sct, sensor_info, spacing_mm, fig_title)
%
%   Displays the body silhouette with the exclusion zone (extruded along the
%   anterior-posterior axis to form a rectangular prism) overlaid in red, in
%   three orthogonal slices through the body centroid:
%     * Transverse (axial, dim 3 fixed): X-Y plane
%     * Sagittal   (dim 2 fixed):        anterior-posterior x sup-inf plane
%     * Coronal    (dim 1 fixed):        left-right x sup-inf plane
%
%   The dim labels follow determine_sensor_mask's convention (dim 1 = Y =
%   anterior-posterior of the patient as the function sees the body mask).
%   No permutation is performed; if the caller's body mask has a different
%   physical-axis layout, the slice labels will be misleading but the
%   overlay geometry will still match.
%
%   INPUTS:
%       sct          - struct with .bodyMask (3D logical)
%       sensor_info  - struct from determine_sensor_mask; uses .exclusion_zone
%                      (2D logical over dim 2 x dim 3 of body mask)
%       spacing_mm   - [dx dy dz] mm (optional; used only for axis labels)
%       fig_title    - optional figure title string
%
%   See also: determine_sensor_mask

    if nargin < 3 || isempty(spacing_mm)
        spacing_mm = [1, 1, 1];
    end
    if nargin < 4 || isempty(fig_title)
        fig_title = 'Beam exclusion zone (red) over body (gray)';
    end

    if ~isfield(sct, 'bodyMask') || isempty(sct.bodyMask)
        warning('plot_exclusion_zone_views:NoBody', 'sct.bodyMask missing; skipping plot.');
        return;
    end
    if ~isfield(sensor_info, 'exclusion_zone') || isempty(sensor_info.exclusion_zone)
        warning('plot_exclusion_zone_views:NoZone', 'sensor_info.exclusion_zone missing; skipping plot.');
        return;
    end

    body    = logical(sct.bodyMask);
    excl_2D = logical(sensor_info.exclusion_zone);

    dims = size(body);
    n1 = dims(1);
    n2 = dims(2);
    n3 = dims(3);

    if ~isequal(size(excl_2D), [n2, n3])
        warning('plot_exclusion_zone_views:SizeMismatch', ...
            'exclusion_zone is %s but body dim 2 x dim 3 is [%d %d]; skipping plot.', ...
            mat2str(size(excl_2D)), n2, n3);
        return;
    end

    % Extrude exclusion zone along dim 1 (anterior-posterior in function coords).
    excl_3D = repmat(reshape(excl_2D, [1, n2, n3]), [n1, 1, 1]);

    % Slice indices: body centroid (fall back to grid center if no body).
    body_idx = find(body);
    if isempty(body_idx)
        s1 = round(n1/2); s2 = round(n2/2); s3 = round(n3/2);
    else
        [i1, i2, i3] = ind2sub(dims, body_idx);
        s1 = round(mean(i1));
        s2 = round(mean(i2));
        s3 = round(mean(i3));
    end

    figure('Name', 'Exclusion zone (3 views)', 'Color', 'w', ...
        'Position', [100 100 1400 450]);

    % --- View 1: dim 3 fixed (transverse / axial in DICOM if dim 3 = Z) ---
    subplot(1, 3, 1);
    show_overlay(squeeze(body(:, :, s3)), squeeze(excl_3D(:, :, s3)), ...
                 spacing_mm(2), spacing_mm(1));
    title(sprintf('Transverse (dim 3 = %d)', s3));
    xlabel('dim 2 (voxels)'); ylabel('dim 1 (voxels)');

    % --- View 2: dim 2 fixed (sagittal if dim 2 = X) ---
    subplot(1, 3, 2);
    show_overlay(squeeze(body(:, s2, :)), squeeze(excl_3D(:, s2, :)), ...
                 spacing_mm(3), spacing_mm(1));
    title(sprintf('Sagittal (dim 2 = %d)', s2));
    xlabel('dim 3 (voxels)'); ylabel('dim 1 (voxels)');

    % --- View 3: dim 1 fixed (coronal if dim 1 = Y) ---
    subplot(1, 3, 3);
    show_overlay(squeeze(body(s1, :, :)), squeeze(excl_3D(s1, :, :)), ...
                 spacing_mm(3), spacing_mm(2));
    title(sprintf('Coronal (dim 1 = %d)', s1));
    xlabel('dim 3 (voxels)'); ylabel('dim 2 (voxels)');

    sgtitle(fig_title);
end


function show_overlay(body_slice, excl_slice, asp_x, asp_y)
%SHOW_OVERLAY Render body in gray with exclusion zone overlaid in red.

    body_d = double(body_slice);
    excl_d = double(excl_slice);

    % RGB composite: gray body, red where exclusion
    R = 0.5 * body_d + excl_d;
    G = 0.5 * body_d .* (1 - excl_d);
    B = 0.5 * body_d .* (1 - excl_d);
    R = min(R, 1);
    img = cat(3, R, G, B);

    image(img);
    axis image; axis xy;
    if nargin >= 4 && asp_x > 0 && asp_y > 0
        daspect([1/asp_y, 1/asp_x, 1]);
    end
    set(gca, 'XTick', [], 'YTick', []);
end
