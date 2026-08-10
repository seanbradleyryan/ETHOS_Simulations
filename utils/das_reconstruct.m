function reconPressure = das_reconstruct(sensorData, sensor, sensor_info, medium, ...
                                         Nx, Ny, Nz, dx, dy, dz, dt, opts)
%DAS_RECONSTRUCT  Delay-and-sum photoacoustic reconstruction for the ETHOS
%  k-Wave simulations, ported from the IRAI 2D-array sample code
%  (sample_code/DAS2D.m, sample_code/pcalculation.m, sample_code/DAS_GG.m).
%
%  reconPressure = das_reconstruct(sensorData, sensor, sensor_info, medium, ...
%                                  Nx, Ny, Nz, dx, dy, dz, dt, opts)
%
%  PURPOSE
%    Back-project recorded acoustic signals onto the simulation voxel grid
%    using the sample delay-and-sum algorithm: per-element channels, a
%    depth/r^2 weighting, a narrow per-element directivity (aperture) cone
%    about the array normal, and a final Hilbert envelope.
%
%  -----------------------------------------------------------------------
%  DATA-FORMAT / ASSUMPTION DIFFERENCES  (sample code  vs  these k-Wave sims)
%  -----------------------------------------------------------------------
%    Sample code (pcalculation.m / DAS2D.m)     These k-Wave sims
%    --------------------------------------     --------------------------
%    PAD[time, 32, 32]  (one trace per          sensorData[nCh, Nt]
%      element, time-first, fixed square          (channel-major, arbitrary
%      grid)                                       channel count)
%    element positions in PIXEL units           physical METRES; arbitrary
%      (elmspace = 0.85*1.54/0.35)                 planar / tilted geometry
%    delay = floor(r_pixels) used directly      delay = r_metres / (c*dt)
%      as a sample index (1 px == 1 sample)
%    weight = depth / r^2                        (old local code used 1/r)
%    aperture: z/r > 0.997 about the +z          (no aperture in old code)
%      depth axis
%    output = Hilbert envelope of a pixel        recon stays on the sim
%      cube centred on the array                   voxel grid [Nx Ny Nz]
%
%  HOW THE BRIDGE IS DONE HERE
%    * Channels = array ELEMENTS. sensorData is collapsed to one trace per
%      element with apply_element_averaging(); each element's position is the
%      centroid of its sensor voxels (metres, in the recon-grid frame) -- the
%      direct analogue of the sample's 32x32 element grid. Falls back to
%      per-voxel channels when no element map is available.
%    * Delay is physical: floor(r/(c*dt)) (or linear interp), using a single
%      homogeneous sound speed c = mean(nonzero medium.sound_speed).
%    * "depth" = signed distance from the array plane along the array normal.
%      The normal is found by PCA on the channel positions (works for flat OR
%      tilted arrays) and oriented to point into the medium. (The sample can
%      assume +z because its array is axis-aligned; the sim array is not.)
%    * Envelope detection runs along the grid axis most aligned with that
%      normal, generalising the sample's "Hilbert along dim 1".
%    * Per-element edge skipping (the sample's iii/ii = 3:30) and its 20-sample
%      near-field blanking are hardware-specific and are NOT reproduced; all
%      available elements and samples >= 1 are used.
%
%  INPUTS
%    sensorData  - [nCh x Nt] recorded pressure. Row order MUST match
%                  find(sensor.mask) (the k-Wave sensor ordering).
%    sensor      - struct with .mask (logical [Nx Ny Nz] on the padded grid).
%    sensor_info - struct from determine_sensor_mask. Uses .voxel_element_idx
%                  and .num_elements when present; otherwise per-voxel fallback.
%    medium      - struct with .sound_speed (scalar or [Nx Ny Nz]).
%    Nx,Ny,Nz    - recon / padded grid size (voxels).
%    dx,dy,dz    - voxel size (metres).
%    dt          - sampling period (seconds).
%    opts        - optional struct (defaults faithful to the sample):
%                    .use_elements (true)   element vs per-voxel channels
%                    .envelope     (true)   Hilbert envelope of the result
%                    .use_aperture (true)   apply the directivity cone
%                    .aperture_cos (0.997)  cone threshold (cos of half-angle)
%                    .depth_weight (true)   depth/r^2 vs plain 1/r^2 weighting
%                    .interp       ('nearest')  'nearest' (floor) | 'linear'
%
%  OUTPUT
%    reconPressure - [Nx x Ny x Nz] double, max(.,0)-clipped. The caller
%                    applies any correction_factor scaling.
%
%  See also: apply_element_averaging, determine_sensor_mask,
%            run_standalone_simulation, run_single_field_simulation

    % ---- 0. Options / defaults (faithful to the sample) -----------------
    if nargin < 12 || isempty(opts), opts = struct(); end
    opts = set_default(opts, 'use_elements', true);
    opts = set_default(opts, 'envelope',     true);
    opts = set_default(opts, 'use_aperture', true);
    opts = set_default(opts, 'aperture_cos', 0.997);
    opts = set_default(opts, 'depth_weight', true);
    opts = set_default(opts, 'interp',       'nearest');

    % ---- 1. Channel signals (element-averaged) + positions (metres) -----
    sensor_lin = find(sensor.mask);
    [vx_i, vy_i, vz_i] = ind2sub([Nx, Ny, Nz], sensor_lin);
    vox_pos = [(double(vx_i) - 1) * dx, ...
               (double(vy_i) - 1) * dy, ...
               (double(vz_i) - 1) * dz];

    have_elements = opts.use_elements ...
        && isfield(sensor_info, 'num_elements') && sensor_info.num_elements > 0 ...
        && isfield(sensor_info, 'voxel_element_idx') ...
        && ~isempty(sensor_info.voxel_element_idx);

    if have_elements
        elem_vec = double(sensor_info.voxel_element_idx(:));
        elem_vec = elem_vec(elem_vec > 0);            % keyed by find(sensor.mask)
        nElem    = sensor_info.num_elements;

        if numel(elem_vec) ~= size(vox_pos, 1)
            warning('das_reconstruct:ElemMismatch', ...
                ['voxel_element_idx valid entries (%d) ~= sensor voxels (%d); ', ...
                 'falling back to per-voxel back-projection.'], ...
                numel(elem_vec), size(vox_pos, 1));
            have_elements = false;
        else
            sd = single(gather(apply_element_averaging(sensorData, sensor_info)));  % [nElem x Nt]
            ch_pos = zeros(nElem, 3);
            keep   = false(nElem, 1);
            for e = 1:nElem
                rows = (elem_vec == e);
                if any(rows)
                    ch_pos(e, :) = mean(vox_pos(rows, :), 1);
                    keep(e) = true;
                end
            end
            ch_pos = ch_pos(keep, :);                 % drop empty elements
            sd     = sd(keep, :);
        end
    end

    if ~have_elements
        sd     = single(gather(sensorData));          % per-voxel channels
        ch_pos = vox_pos;
    end

    nCh   = size(ch_pos, 1);
    Nt_sd = size(sd, 2);
    if nCh == 0 || Nt_sd == 0
        warning('das_reconstruct:NoData', 'No channels/samples; returning zeros.');
        reconPressure = zeros(Nx, Ny, Nz);
        return;
    end

    % ---- 2. Homogeneous sound speed -------------------------------------
    c_vals = medium.sound_speed(:);
    c_vals = c_vals(c_vals > 0);
    c_das  = double(gather(mean(c_vals)));
    if ~isfinite(c_das) || c_das <= 0
        error('das_reconstruct:BadSoundSpeed', 'Invalid mean sound speed (%.3g).', c_das);
    end

    % ---- 3. Array normal via PCA (flat or tilted), oriented into medium -
    arrayCentroid = mean(ch_pos, 1);
    n = array_normal_pca(ch_pos);                     % unit 1x3
    gridCentroid = [(Nx - 1) * dx, (Ny - 1) * dy, (Nz - 1) * dz] / 2;
    if dot(gridCentroid - arrayCentroid, n) < 0
        n = -n;                                       % point into the medium
    end

    % ---- 4. Voxel grid (metres) + precomputed projection / depth --------
    [Xg, Yg, Zg] = ndgrid((0:Nx-1) * dx, (0:Ny-1) * dy, (0:Nz-1) * dz);
    Xg = single(Xg);  Yg = single(Yg);  Zg = single(Zg);
    nS      = single(n);
    proj    = Xg * nS(1) + Yg * nS(2) + Zg * nS(3);   % voxel . n
    depth_v = proj - single(dot(arrayCentroid, n));   % distance from array plane
    depth_pos = depth_v > 0;                          % only in front of the array

    % ---- 5. Back-projection over channels (memory-bounded, single) ------
    reconPressure = zeros(Nx, Ny, Nz, 'single');
    dmin       = single(min([dx, dy, dz]));
    inv_cdt    = single(1 / (c_das * dt));
    use_linear = strcmpi(opts.interp, 'linear');
    use_ap     = logical(opts.use_aperture);
    ap_cos     = single(opts.aperture_cos);
    use_dw     = logical(opts.depth_weight);

    if have_elements, ch_kind = 'elements'; else, ch_kind = 'voxels'; end
    fprintf('         DAS: %d channels (%s), c = %.1f m/s, normal = [%.2f %.2f %.2f]\n', ...
        nCh, ch_kind, c_das, n(1), n(2), n(3));

    for ch = 1:nCh
        px = single(ch_pos(ch, 1));
        py = single(ch_pos(ch, 2));
        pz = single(ch_pos(ch, 3));

        r     = sqrt((Xg - px).^2 + (Yg - py).^2 + (Zg - pz).^2);
        rsafe = max(r, dmin);

        % Aperture: cos(angle between element->voxel ray and array normal)
        cosA  = (proj - single(dot(ch_pos(ch, :), n))) ./ rsafe;

        idx_f = r * inv_cdt + 1;                      % fractional sample index
        idx0  = floor(idx_f);

        valid = depth_pos & (idx0 >= 1) & (idx0 < Nt_sd);
        if use_ap
            valid = valid & (cosA > ap_cos);
        end
        if ~any(valid(:)), continue; end

        samp = zeros(Nx, Ny, Nz, 'single');
        i0   = idx0(valid);
        if use_linear
            frac = idx_f(valid) - i0;
            v0   = sd(ch, i0);
            v1   = sd(ch, i0 + 1);
            samp(valid) = (1 - frac) .* v0(:) + frac .* v1(:);
        else
            v0 = sd(ch, i0);
            samp(valid) = v0(:);                      % nearest (floor), as in the sample
        end

        if use_dw
            w = depth_v ./ (rsafe.^2);                % depth / r^2  (sample weighting)
        else
            w = single(1) ./ (rsafe.^2);
        end

        reconPressure = reconPressure + w .* samp;    % samp == 0 outside 'valid'

        if mod(ch, max(1, round(nCh / 10))) == 0
            fprintf('         DAS channel %d/%d\n', ch, nCh);
        end
    end

    % ---- 6. Envelope detection (Hilbert) along the normal-aligned axis --
    if opts.envelope
        [~, ax] = max(abs(n));
        reconPressure = envelope_along_axis(reconPressure, ax);
    end

    % ---- 7. Finalise ----------------------------------------------------
    reconPressure = max(double(real(reconPressure)), 0);
end


% ========================================================================
% Local helpers
% ========================================================================

function s = set_default(s, field, val)
%SET_DEFAULT  Fill s.(field) with val only when it is missing or empty.
    if ~isfield(s, field) || isempty(s.(field))
        s.(field) = val;
    end
end


function n = array_normal_pca(pos)
%ARRAY_NORMAL_PCA  Unit normal of the best-fit plane through points pos (N x 3).
%  The normal is the eigenvector of the smallest eigenvalue of the centred
%  covariance (the least-variance direction of a planar point cloud).
    if size(pos, 1) < 3
        n = [1, 0, 0];                                % degenerate -> default X
        return;
    end
    p0 = pos - mean(pos, 1);
    C  = (p0' * p0) / size(pos, 1);
    [V, D] = eig((C + C') / 2);                       % symmetric -> real eigvecs
    [~, k] = min(diag(D));
    n  = V(:, k)';
    nn = norm(n);
    if nn < eps
        n = [1, 0, 0];
    else
        n = n / nn;
    end
end


function img = envelope_along_axis(img, ax)
%ENVELOPE_ALONG_AXIS  abs(hilbert()) of img taken along dimension ax.
%  hilbert() operates column-wise (dim 1), so move ax -> 1, reshape to a
%  matrix, transform, then restore the original shape.
    img  = double(img);
    nd   = max(ndims(img), ax);
    perm = 1:nd;
    perm([1, ax]) = [ax, 1];
    imgp = permute(img, perm);
    sz   = size(imgp);
    imgp = reshape(imgp, sz(1), []);
    imgp = abs(hilbert(imgp));
    imgp = reshape(imgp, sz);
    img  = single(ipermute(imgp, perm));
end
