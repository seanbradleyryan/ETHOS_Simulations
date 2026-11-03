function [reconPressure, operator] = das_reconstruct(sensorData, sensor, sensor_info, medium, ...
                                                    Nx, Ny, Nz, dx, dy, dz, dt, opts)
%DAS_RECONSTRUCT  Delay-and-sum photoacoustic reconstruction for the ETHOS
%  k-Wave simulations, ported from the IRAI 2D-array sample code
%  (sample_code/DAS2D.m, sample_code/pcalculation.m, sample_code/DAS_GG.m).
%
%  [reconPressure, operator] = das_reconstruct(sensorData, sensor, sensor_info, medium, ...
%                                              Nx, Ny, Nz, dx, dy, dz, dt, opts)
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
%  PRECOMPUTED DAS MATRIX (OPERATOR)
%    Before the envelope, DAS is a fixed linear map from the data to the
%    image: each voxel is a weighted sum of a few samples (its delay) from the
%    few elements whose cone contains it. That map depends only on the channel
%    positions, the grid, dt, Nt, the sound speed and the weighting options --
%    not on the dose -- so it is built once as a sparse matrix A of size
%    (Nx*Ny*Nz) x (nCh*Nt) and reused: image = A * sd(:). The matrix gives the
%    same image as the channel loop up to rounding (the sums run in double).
%    The path is chosen in this order:
%      1. opts.operator given        -> multiply with it (its key must match).
%      2. opts.operator_cache_dir has das_operator_cache_<key>.mat -> load it.
%      3. Build it, but only when it will be reused (a cache dir is set, or the
%         caller asks for the 2nd output): one build costs about one channel
%         loop. With a cache dir it is saved there (temporary name + rename, so
%         no reader sees a half-written file).
%      4. Otherwise (or opts.use_operator = false, or the matrix would need more
%         than opts.max_operator_gb to build) -> the channel loop.
%    The key is an MD5 of everything the matrix depends on, so a different
%    sensor, grid, dt/Nt, sound speed or DAS setting gets its own matrix.
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
%    opts        - optional struct. Missing DAS settings take the das_options
%                  defaults (faithful to the sample):
%                    .use_elements (true)   element vs per-voxel channels
%                    .envelope     (true)   Hilbert envelope of the result
%                    .use_aperture (true)   apply the directivity cone
%                    .aperture_cos (0.997)  cone threshold (cos of half-angle)
%                    .depth_weight (true)   depth/r^2 vs plain 1/r^2 weighting
%                    .interp       ('nearest')  'nearest' (floor) | 'linear'
%                    .use_operator (true)   allow the precomputed DAS matrix
%                    .max_operator_gb (8)   memory limit for building it
%                    .check_operator (false) also run the channel loop and
%                                           print the largest difference
%                    .use_gpu      (false)  multiply the matrix on the GPU
%                                           (falls back to the CPU on error)
%                  and, not read from CONFIG:
%                    .operator     ([])     matrix from an earlier call
%                    .operator_cache_dir ('') folder of the matrix cache
%                    .clip_negative (true)  max(.,0) on the output; false
%                                           keeps the signed result (with
%                                           envelope = false the operator is
%                                           then linear, e.g. for use as an
%                                           approximate adjoint in demo_mbr)
%
%  OUTPUT
%    reconPressure - [Nx x Ny x Nz] double, max(.,0)-clipped unless
%                    opts.clip_negative = false. The caller applies any
%                    correction_factor scaling.
%    operator      - The DAS matrix used (struct: .key, .A), or [] when the
%                    channel loop ran. Pass it back as opts.operator to skip
%                    loading.
%
%  EXAMPLE
%    opts = das_options(CONFIG);
%    opts.operator_cache_dir = get_simulation_dir(CONFIG.working_dir, pid, session);
%    recon = das_reconstruct(sensorData, sensor, sensor_info, medium, ...
%                            Nx, Ny, Nz, dx, dy, dz, dt, opts);
%
%  See also: das_options, apply_element_averaging, determine_sensor_mask,
%            run_standalone_simulation, run_single_field_simulation

    % ---- 0. Options / defaults (faithful to the sample) -----------------
    if nargin < 12 || isempty(opts), opts = struct(); end
    defaults = das_options(struct());
    names    = fieldnames(defaults);
    for k = 1:numel(names)
        opts = set_default(opts, names{k}, defaults.(names{k}));
    end
    opts = set_default(opts, 'clip_negative', true);
    opts = set_default(opts, 'operator', []);
    opts = set_default(opts, 'operator_cache_dir', '');

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
    operator = [];
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

    if have_elements, ch_kind = 'elements'; else, ch_kind = 'voxels'; end
    fprintf('         DAS: %d channels (%s), c = %.1f m/s, normal = [%.2f %.2f %.2f]\n', ...
        nCh, ch_kind, c_das, n(1), n(2), n(3));

    % ---- 4. DAS matrix: given, cached on disk, or built for reuse -------
    if opts.use_operator
        key = operator_key(ch_pos, Nx, Ny, Nz, dx, dy, dz, dt, Nt_sd, c_das, opts);
        cache_file = '';
        if ~isempty(opts.operator_cache_dir)
            cache_file = fullfile(opts.operator_cache_dir, sprintf('das_operator_cache_%s.mat', key));
        end

        if ~isempty(opts.operator)
            if ~isfield(opts.operator, 'key') || ~strcmp(opts.operator.key, key)
                error('das_reconstruct:OperatorMismatch', ...
                    ['opts.operator was built for a different sensor, grid, time step, ' ...
                     'sound speed or DAS setting (this call needs key %s).'], key);
            end
            operator = opts.operator;

        elseif ~isempty(cache_file) && isfile(cache_file)
            try
                load_tic = tic;
                loaded   = load(cache_file, 'das_operator');
                operator = loaded.das_operator;
                if ~strcmp(operator.key, key)
                    error('das_reconstruct:CacheKeyMismatch', 'stored key %s', operator.key);
                end
                fprintf('         DAS matrix loaded from %s (%.1f s)\n', cache_file, toc(load_tic));
            catch ME
                % Unreadable (e.g. a broken save): rebuild below and overwrite it
                warning('das_reconstruct:CacheReadFail', ...
                    'Could not read DAS matrix cache %s (%s); rebuilding it.', cache_file, ME.message);
                operator = [];
            end
        end

        % Build only when it will be reused: a single build costs about one
        % channel loop, so a one-off call is faster with the loop.
        if isempty(operator) && (~isempty(cache_file) || nargout > 1)
            build_tic = tic;
            [~, A] = das_channel_loop(sd, ch_pos, n, arrayCentroid, Nx, Ny, Nz, ...
                dx, dy, dz, c_das, dt, opts, true);
            if ~isempty(A)
                operator = struct('key', key, 'A', A);
                fprintf('         DAS matrix built: %d entries (about %.2f GB) in %.1f s\n', ...
                    nnz(A), nnz(A) * 16 / 1e9, toc(build_tic));
                clear A;
                if ~isempty(cache_file)
                    try
                        if ~isfolder(opts.operator_cache_dir)
                            mkdir(opts.operator_cache_dir);
                        end
                        % Save under a unique temporary name, then rename. The
                        % rename is atomic, so no worker reads a half-written file.
                        [~, tmp_name] = fileparts(tempname);
                        tmp_file = fullfile(opts.operator_cache_dir, [tmp_name '.mat']);
                        das_operator = operator;  %#ok<NASGU>  saved variable name
                        save(tmp_file, 'das_operator', '-v7.3', '-nocompression');
                        movefile(tmp_file, cache_file, 'f');
                        fprintf('         DAS matrix saved to %s\n', cache_file);
                    catch ME
                        warning('das_reconstruct:CacheSaveFail', ...
                            'Could not save DAS matrix cache %s: %s', cache_file, ME.message);
                    end
                end
            end
        end
    end

    % ---- 5. Back-projection: one matrix multiply, or the channel loop ---
    if ~isempty(operator)
        mult_tic = tic;
        x = double(sd(:));                            % column index = ch + (sample-1)*nCh
        on_gpu = false;
        if opts.use_gpu
            try
                y = gather(gpuArray(operator.A) * gpuArray(x));
                on_gpu = true;
            catch ME
                warning('das_reconstruct:GpuMultiplyFail', ...
                    'GPU multiply failed (%s); using the CPU.', ME.message);
            end
        end
        if ~on_gpu
            y = operator.A * x;
        end
        reconPressure = reshape(single(y), Nx, Ny, Nz);
        if on_gpu, where = 'GPU'; else, where = 'CPU'; end
        fprintf('         DAS matrix multiply on the %s (%.2f s)\n', where, toc(mult_tic));

        if opts.check_operator
            recon_loop = das_channel_loop(sd, ch_pos, n, arrayCentroid, Nx, Ny, Nz, ...
                dx, dy, dz, c_das, dt, opts, false);
            rel_diff = max(abs(reconPressure(:) - recon_loop(:))) / max(max(abs(recon_loop(:))), eps('single'));
            fprintf('         DAS matrix check: max |matrix - loop| / max |loop| = %.2e\n', rel_diff);
            clear recon_loop;
        end
    else
        reconPressure = das_channel_loop(sd, ch_pos, n, arrayCentroid, Nx, Ny, Nz, ...
            dx, dy, dz, c_das, dt, opts, false);
    end

    % ---- 6. Envelope detection (Hilbert) along the normal-aligned axis --
    if opts.envelope
        [~, ax] = max(abs(n));
        reconPressure = envelope_along_axis(reconPressure, ax);
    end

    % ---- 7. Finalise ----------------------------------------------------
    reconPressure = double(real(reconPressure));
    if opts.clip_negative
        reconPressure = max(reconPressure, 0);
    end
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


function [reconPressure, A] = das_channel_loop(sd, ch_pos, n, arrayCentroid, Nx, Ny, Nz, ...
                                               dx, dy, dz, c_das, dt, opts, build_matrix)
%DAS_CHANNEL_LOOP  The delay-and-sum loop over channels, in one of two modes.
%  build_matrix = false: sum every channel's delayed, weighted samples into
%    reconPressure ([Nx Ny Nz] single); A = [].
%  build_matrix = true: record each (voxel, sample, weight) instead and return
%    the sparse DAS matrix A ((Nx*Ny*Nz) x (nCh*Nt)); reconPressure = []. The
%    channel nearest the array centre runs first: its entry count x nCh
%    estimates the matrix size, and when building it would need more than
%    opts.max_operator_gb (about 52 bytes per entry while building) A = [] is
%    returned so the caller uses the summing mode.
%  Both modes use the same geometry lines, so A * sd(:) equals the summed image
%  up to rounding.

    % Voxel grid (metres) + precomputed projection / depth
    [Xg, Yg, Zg] = ndgrid((0:Nx-1) * dx, (0:Ny-1) * dy, (0:Nz-1) * dz);
    Xg = single(Xg);  Yg = single(Yg);  Zg = single(Zg);
    nS      = single(n);
    proj    = Xg * nS(1) + Yg * nS(2) + Zg * nS(3);   % voxel . n
    depth_v = proj - single(dot(arrayCentroid, n));   % distance from array plane
    depth_pos = depth_v > 0;                          % only in front of the array

    nCh        = size(ch_pos, 1);
    Nt_sd      = size(sd, 2);
    dmin       = single(min([dx, dy, dz]));
    inv_cdt    = single(1 / (c_das * dt));
    use_linear = strcmpi(opts.interp, 'linear');
    use_ap     = logical(opts.use_aperture);
    ap_cos     = single(opts.aperture_cos);
    use_dw     = logical(opts.depth_weight);

    reconPressure = [];
    A             = [];
    if build_matrix
        [~, ch_mid] = min(sum((ch_pos - arrayCentroid) .^ 2, 2));
        order  = [ch_mid, 1:ch_mid-1, ch_mid+1:nCh];
        rows_c = cell(nCh, 1);   % voxel linear index (uint32)
        cols_c = cell(nCh, 1);   % column of A = ch + (sample-1)*nCh (uint32)
        vals_c = cell(nCh, 1);   % weight (single)
    else
        order = 1:nCh;
        reconPressure = zeros(Nx, Ny, Nz, 'single');
    end

    for k = 1:nCh
        ch = order(k);
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

        if any(valid(:))
            if build_matrix
                lin = find(valid);
                i0  = idx0(valid);
                if use_dw
                    w = depth_v(valid) ./ (rsafe(valid).^2);   % depth / r^2  (sample weighting)
                else
                    w = single(1) ./ (rsafe(valid).^2);
                end
                col0 = ch + (double(i0) - 1) * nCh;
                if use_linear
                    frac = idx_f(valid) - i0;
                    rows_c{ch} = uint32([lin; lin]);
                    cols_c{ch} = uint32([col0; col0 + nCh]);  % samples i0 and i0+1
                    vals_c{ch} = [(1 - frac) .* w; frac .* w];
                else
                    rows_c{ch} = uint32(lin);
                    cols_c{ch} = uint32(col0);
                    vals_c{ch} = w;
                end
            else
                samp = zeros(Nx, Ny, Nz, 'single');
                i0   = idx0(valid);
                if use_linear
                    frac = idx_f(valid) - i0;
                    v0   = sd(ch, i0);
                    v1   = sd(ch, i0 + 1);
                    samp(valid) = (1 - frac) .* v0(:) + frac .* v1(:);
                else
                    v0 = sd(ch, i0);
                    samp(valid) = v0(:);              % nearest (floor), as in the sample
                end

                if use_dw
                    w = depth_v ./ (rsafe.^2);        % depth / r^2  (sample weighting)
                else
                    w = single(1) ./ (rsafe.^2);
                end

                reconPressure = reconPressure + w .* samp;    % samp == 0 outside 'valid'
            end
        end

        % Size check after the centre channel: stop before using too much memory
        if build_matrix && k == 1
            est_entries = numel(vals_c{ch}) * nCh;
            est_gb      = est_entries * 52 / 1e9;
            if est_gb > opts.max_operator_gb
                warning('das_reconstruct:OperatorTooLarge', ...
                    ['DAS matrix would need about %.1f GB to build (limit %.1f GB, ' ...
                     'das_operator_max_gb); using the channel loop.'], est_gb, opts.max_operator_gb);
                return;
            end
        end

        if mod(k, max(1, round(nCh / 10))) == 0
            fprintf('         DAS channel %d/%d\n', k, nCh);
        end
    end

    if build_matrix
        clear Xg Yg Zg proj depth_v depth_pos r rsafe cosA idx_f idx0 valid
        rows = double(vertcat(rows_c{:}));  clear rows_c
        cols = double(vertcat(cols_c{:}));  clear cols_c
        vals = double(vertcat(vals_c{:}));  clear vals_c
        A = sparse(rows, cols, vals, Nx * Ny * Nz, nCh * Nt_sd);
    end
end


function key = operator_key(ch_pos, Nx, Ny, Nz, dx, dy, dz, dt, Nt_sd, c_das, opts)
%OPERATOR_KEY  8-character MD5 key of everything the DAS matrix depends on:
%  channel positions (sensor + element map), grid size and spacing, dt, Nt,
%  sound speed and the DAS settings that enter the weights. The envelope and
%  clipping act after the multiply, so they are left out.
    parts = {ch_pos, [Nx, Ny, Nz, dx, dy, dz, dt, Nt_sd, c_das], ...
             double([opts.use_elements, opts.use_aperture, opts.depth_weight]), ...
             opts.aperture_cos, double(lower(char(opts.interp)))};
    md = java.security.MessageDigest.getInstance('MD5');
    for k = 1:numel(parts)
        md.update(typecast(double(parts{k}(:)), 'uint8'));
    end
    raw = typecast(md.digest(), 'uint8');
    key = sprintf('%02x', raw(1:4));
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
