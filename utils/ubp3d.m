function reconPressure = ubp3d(sensorData, sensor, sensor_info, medium, ...
                              Nx, Ny, Nz, dx, dy, dz, dt, opts)
%UBP3D  Universal back-projection photoacoustic reconstruction (Xu & Wang 2005).
%
%   reconPressure = ubp3d(sensorData, sensor, sensor_info, medium, ...
%                         Nx, Ny, Nz, dx, dy, dz, dt, opts)
%
%   PURPOSE
%     3D universal back-projection for an ARBITRARY sensor mask + sensor data.
%     Each channel's trace is filtered with the paper's kernel
%         b(r0,t) = 2 p(r0,t) - 2 t d/dt p(r0,t)
%     and back-projected onto the spherical shells t = |r - r0|/c, weighted by
%     the 3D solid angle dOmega = dA * cos(theta) / |r - r0|^2 (theta = angle
%     between the detector-surface normal and the ray to the voxel), then
%     normalised by Omega0 = 2*pi (planar aperture). Unlike delay-and-sum this
%     is an exact inverse for the canonical geometries: filled regions come back
%     filled (no bipolar/edge-only artefact) and amplitudes are quantitative.
%
%     Signature and conventions mirror DAS_RECONSTRUCT exactly, so ubp3d is a
%     drop-in reconstruction option for run_single_field_simulation (alongside
%     'tr' and 'das'). Sensor positions come from find(sensor.mask); channels
%     are array elements when an element map is available, else per-voxel.
%
%   NON-IDEAL SENSOR CORRECTIONS (auto when sensor_info has an element map, i.e.
%   determine_sensor_mask was used):
%     * Element averaging  - apply_element_averaging() collapses the sensor
%       voxels to one trace per finite element at its centroid (models the
%       element's finite receive area). Ideal full-plane sensors (no element
%       map) fall back to per-voxel channels.
%     * Bandwidth limitation - the forward chain band-limits the signal (finite
%       transducer response) BEFORE reconstruction, so a broadband derivative
%       would amplify only out-of-band noise. The UBP time-derivative is instead
%       computed with a regularised, band-limited differentiator whose cutoff is
%       estimated from the data spectrum ('auto'), keeping resolution where there
%       is signal and rejecting noise where there is none.
%     * Directivity - the physical cos(theta) obliquity weight (and an optional
%       aperture cone) down-weight off-normal rays, as a finite element does.
%
%   INPUTS  (identical layout to das_reconstruct)
%     sensorData  - [nCh x Nt] recorded pressure. Row order MUST match
%                   find(sensor.mask) (the k-Wave sensor ordering).
%     sensor      - struct with .mask (logical [Nx Ny Nz] on the recon grid).
%     sensor_info - struct from determine_sensor_mask; uses .voxel_element_idx
%                   and .num_elements when present (else per-voxel). Pass [] to
%                   force per-voxel channels (a bare arbitrary mask).
%     medium      - struct with .sound_speed (scalar or [Nx Ny Nz]); a plain
%                   numeric scalar sound speed is also accepted.
%     Nx,Ny,Nz    - recon grid size (voxels).
%     dx,dy,dz    - voxel size (metres).
%     dt          - sampling period (seconds).
%     opts        - optional struct:
%                     .use_elements (true)   element vs per-voxel channels
%                     .band_limit   ('auto') UBP derivative regularisation:
%                                            'auto' | 'none' | cutoff in Hz
%                     .aperture_cos (0)      keep rays with cos(theta) > this
%                                            (0 = front half-space; raise to
%                                            model a narrower element directivity)
%                     .area_weight  (true)   weight elements by voxel count (area)
%                     .interp       ('linear') 'linear' | 'nearest' time sampling
%
%   OUTPUT
%     reconPressure - [Nx x Ny x Nz] double, max(.,0)-clipped. The caller
%                     applies any correction_factor / pressure scaling.
%
%   NOTE (>200 lines): kept as one file so the UBP filter, the non-ideal
%   corrections and the back-projection read together, matching das_reconstruct.
%
%   See also: das_reconstruct, apply_element_averaging, determine_sensor_mask,
%             run_single_field_simulation

    % ---- 0. Options / defaults ------------------------------------------
    if nargin < 12 || isempty(opts), opts = struct(); end
    if nargin < 3, sensor_info = struct(); end
    if isempty(sensor_info), sensor_info = struct(); end
    opts = set_default(opts, 'use_elements', true);
    opts = set_default(opts, 'band_limit',   'auto');
    opts = set_default(opts, 'aperture_cos', 0);
    opts = set_default(opts, 'area_weight',  true);
    opts = set_default(opts, 'interp',       'linear');

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
            warning('ubp3d:ElemMismatch', ...
                ['voxel_element_idx valid entries (%d) ~= sensor voxels (%d); ', ...
                 'falling back to per-voxel back-projection.'], ...
                numel(elem_vec), size(vox_pos, 1));
            have_elements = false;
        else
            sd     = single(gather(apply_element_averaging(sensorData, sensor_info))); % [nElem x Nt]
            ch_pos = zeros(nElem, 3);
            ch_cnt = zeros(nElem, 1);
            keep   = false(nElem, 1);
            for e = 1:nElem
                rows = (elem_vec == e);
                if any(rows)
                    ch_pos(e, :) = mean(vox_pos(rows, :), 1);   % element centroid
                    ch_cnt(e)    = nnz(rows);                   % voxels ~ area
                    keep(e)      = true;
                end
            end
            ch_pos = ch_pos(keep, :);
            ch_cnt = ch_cnt(keep);
            sd     = sd(keep, :);
        end
    end

    if ~have_elements
        sd     = single(gather(sensorData));          % per-voxel channels
        ch_pos = vox_pos;
        ch_cnt = ones(size(vox_pos, 1), 1);
    end

    nCh   = size(ch_pos, 1);
    Nt_sd = size(sd, 2);
    if nCh == 0 || Nt_sd < 2
        warning('ubp3d:NoData', 'No channels/samples; returning zeros.');
        reconPressure = zeros(Nx, Ny, Nz);
        return;
    end

    % Per-channel receive area (relative element area, scaled by a voxel face).
    % Absolute scale is calibrated downstream (correction_factor), so only the
    % relative weighting matters here.
    if opts.area_weight
        ch_area = single(ch_cnt(:) * (dx * dy));
    else
        ch_area = single(ones(nCh, 1) * (dx * dy));
    end

    % ---- 2. Homogeneous sound speed -------------------------------------
    if isstruct(medium)
        c_vals = medium.sound_speed(:);
        c_vals = c_vals(c_vals > 0);
        c_das  = double(gather(mean(c_vals)));
    else
        c_das  = double(medium);                      % plain scalar sound speed
    end
    if ~isfinite(c_das) || c_das <= 0
        error('ubp3d:BadSoundSpeed', 'Invalid mean sound speed (%.3g).', c_das);
    end

    % ---- 3. UBP temporal filter  b = 2*(p - t*dp/dt) --------------------
    %  The derivative is band-limited (see header) so it does not amplify the
    %  out-of-band noise left by the finite transducer response.
    tt   = (0:Nt_sd - 1) * dt;                         % sample times (r/c grid)
    dpdt = ubp_time_derivative(sd, dt, opts.band_limit);
    b    = single(2 * double(sd) - 2 * (tt .* dpdt));  % [nCh x Nt]

    % ---- 4. Array normal via PCA (flat or tilted), oriented into medium -
    arrayCentroid = mean(ch_pos, 1);
    n = array_normal_pca(ch_pos);                      % unit 1x3
    gridCentroid = [(Nx - 1) * dx, (Ny - 1) * dy, (Nz - 1) * dz] / 2;
    if dot(gridCentroid - arrayCentroid, n) < 0
        n = -n;                                        % point into the medium
    end

    % ---- 5. Voxel grid (metres) + projection along the array normal -----
    [Xg, Yg, Zg] = ndgrid((0:Nx-1) * dx, (0:Ny-1) * dy, (0:Nz-1) * dz);
    Xg = single(Xg);  Yg = single(Yg);  Zg = single(Zg);
    nS   = single(n);
    proj = Xg * nS(1) + Yg * nS(2) + Zg * nS(3);       % voxel . n

    % ---- 6. Back-projection over channels (memory-bounded, single) ------
    Omega0     = single(2 * pi);                       % planar-aperture solid angle
    recon      = zeros(Nx, Ny, Nz, 'single');
    dmin       = single(min([dx, dy, dz]));
    inv_cdt    = single(1 / (c_das * dt));
    ap_cos     = single(opts.aperture_cos);
    use_linear = strcmpi(opts.interp, 'linear');

    if have_elements, ch_kind = 'elements'; else, ch_kind = 'voxels'; end
    fprintf('         UBP: %d channels (%s), c = %.1f m/s, normal = [%.2f %.2f %.2f], deriv = %s\n', ...
        nCh, ch_kind, c_das, n(1), n(2), n(3), band_limit_label(opts.band_limit));

    for ch = 1:nCh
        r     = sqrt((Xg - single(ch_pos(ch,1))).^2 + ...
                     (Yg - single(ch_pos(ch,2))).^2 + ...
                     (Zg - single(ch_pos(ch,3))).^2);
        rsafe = max(r, dmin);

        % cos(theta): obliquity between the array normal and the ray to the voxel
        cosT  = (proj - single(dot(ch_pos(ch,:), n))) ./ rsafe;

        idx_f = r * inv_cdt + 1;                       % fractional sample index
        idx0  = floor(idx_f);
        valid = (cosT > ap_cos) & (idx0 >= 1) & (idx0 < Nt_sd);
        if ~any(valid(:)), continue; end

        bsamp = zeros(Nx, Ny, Nz, 'single');
        i0    = idx0(valid);
        if use_linear
            frac = idx_f(valid) - i0;
            v0   = b(ch, i0);
            v1   = b(ch, i0 + 1);
            bsamp(valid) = (1 - frac) .* v0(:) + frac .* v1(:);
        else
            v0 = b(ch, i0);
            bsamp(valid) = v0(:);
        end

        w = ch_area(ch) * cosT ./ (rsafe.^2);          % dA * cos(theta) / r^2
        recon = recon + w .* bsamp;                    % bsamp == 0 outside 'valid'

        if mod(ch, max(1, round(nCh / 10))) == 0
            fprintf('         UBP channel %d/%d\n', ch, nCh);
        end
    end

    % ---- 7. Normalise and finalise --------------------------------------
    recon = recon / Omega0;
    reconPressure = max(double(recon), 0);
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


function dpdt = ubp_time_derivative(sd, dt, mode)
%UBP_TIME_DERIVATIVE  d/dt of each row of sd, optionally band-limited.
%   mode 'none'  -> central finite difference (broadband).
%   mode 'auto'  -> spectral derivative, Gaussian low-pass at the 99%-power
%                   cutoff of the data (regularises the finite transducer band).
%   mode numeric -> spectral derivative, Gaussian low-pass at that cutoff (Hz).
    sd = double(sd);
    Nt = size(sd, 2);

    if ischar(mode) && strcmpi(mode, 'none')
        dpdt              = zeros(size(sd));
        dpdt(:, 2:end-1)  = (sd(:, 3:end) - sd(:, 1:end-2)) / (2 * dt);
        dpdt(:, 1)        = (sd(:, 2)   - sd(:, 1))     / dt;
        dpdt(:, end)      = (sd(:, end) - sd(:, end-1)) / dt;
        return;
    end

    f = [0:floor(Nt/2), -(ceil(Nt/2)-1):-1] / (Nt * dt);   % 1 x Nt signed freqs
    if ischar(mode) && strcmpi(mode, 'auto')
        fc = estimate_bandwidth(sd, dt, 0.99);
    else
        fc = double(mode);
    end
    if ~isfinite(fc) || fc <= 0
        G = ones(1, Nt);
    else
        G = exp(-(f / fc).^2);                             % Gaussian low-pass
    end
    D    = (1i * 2 * pi * f) .* G;                          % band-limited d/dt
    dpdt = real(ifft(fft(sd, [], 2) .* D, [], 2));
end


function fc = estimate_bandwidth(sd, dt, frac)
%ESTIMATE_BANDWIDTH  Frequency below which `frac` of the mean spectral power
%  lies. Used to set the UBP derivative cutoff to the data's real bandwidth.
    Nt   = size(sd, 2);
    P    = mean(abs(fft(sd, [], 2)).^2, 1);
    half = floor(Nt/2) + 1;
    Ppos = P(1:half);
    Ppos(1) = 0;                        % ignore DC: the derivative kills it anyway
    fpos = (0:half-1) / (Nt * dt);
    cP   = cumsum(Ppos) / max(sum(Ppos), eps);
    k    = find(cP >= frac, 1, 'first');
    if isempty(k), fc = fpos(end); else, fc = fpos(k); end
    if fc <= 0, fc = fpos(end); end
end


function lbl = band_limit_label(mode)
%BAND_LIMIT_LABEL  Short printable description of the derivative mode.
    if ischar(mode)
        lbl = mode;
    else
        lbl = sprintf('%.2f MHz', double(mode) / 1e6);
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
