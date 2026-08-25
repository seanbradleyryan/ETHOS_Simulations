function out = misc_UBP(num_sensors)
%MISC_UBP  Universal back-projection (UBP) photoacoustic demo in a fat phantom.
%
%   out = misc_UBP()            % 60 sensors (default)
%   out = misc_UBP(num_sensors) % choose the number of ring sensors
%
%   Set CONFIG.movie = true (below) to instead sweep the sensor count from
%   CONFIG.movie_min_sensors to CONFIG.movie_max_sensors and animate the
%   reconstruction, with a slider to swoop through num_sensors by hand.
%
%   PURPOSE:
%       Stand-alone teaching / sanity demo of the Universal Back-Projection
%       reconstruction of Xu & Wang, "Universal back-projection algorithm for
%       photoacoustic computed tomography", Phys. Rev. E 71, 016706 (2005).
%       It builds a homogeneous fat phantom, drops an arbitrary-but-sensible
%       radiation-beam pressure source into it, records the acoustic field on a
%       ring of point sensors with a k-Wave forward simulation, and reconstructs
%       the pressure in the plane of that ring with UBP. It then reports the
%       spatial resolution of the reconstruction.
%
%   GEOMETRY (all in the 2D sensor plane):
%       - The computational grid IS the plane of the sensor circle.
%       - The radiation beam travels ORTHOGONAL to this plane ("down the
%         phantom"), so what lives in the plane is the beam's CROSS-SECTION: a
%         flat-topped square field with a soft penumbra, centred on the ring.
%       - A cylindrical fat block (a disc, in cross-section) sits in a water
%         background. The sensor ring and the whole beam are INSIDE the fat, so
%         every source->sensor path is homogeneous fat (sound speed c_fat) and
%         the single-speed UBP formula applies exactly.
%
%   WHY UBP AND NOT CT BACK-PROJECTION:
%       A point sensor at r0 that fires at time t constrains the source to the
%       shell |r - r0| = c*t -- a SPHERE (a CIRCLE in this plane), not a line.
%       UBP therefore smears each filtered sensor trace back over circular arcs
%       and sums them, with the paper's filter  b = 2p - 2t dp/dt  and a
%       solid-angle weight. CT/Radon back-projection instead smears along
%       straight lines, which is wrong for a diverging acoustic wavefront.
%
%   MOVIE MODE (CONFIG.movie = true):
%       ONE dense forward simulation is run at CONFIG.movie_max_sensors. Because
%       point sensors are passive receivers, a SUBSET of those traces is identical
%       to having simulated only that subset, so the whole sweep costs one sim.
%       By default frame counts are the divisors of movie_max_sensors (each frame
%       exactly equidistant). Set CONFIG.movie_all_counts = true to show EVERY
%       count from movie_min_sensors to movie_max_sensors (step movie_step); non-
%       divisor counts pick the nearest dense sensors (<= half the dense spacing).
%
%   INPUTS:
%       num_sensors - (optional) number of equally-spaced point sensors on the
%                     ring. Positive integer, default 60. Ignored in movie mode.
%
%   OUTPUTS:
%       out - single mode: struct with .recon/.p0_true/.img_x/.img_y/.sensor_xy/
%             .resolution_mm/.resolution_edge_mm/.config.
%             movie  mode: struct with .frames (cell of recon images),
%             .num_sensors, .resolution_mm (per frame), .img_x/.img_y/.config.
%
%   EXAMPLE:
%       out = misc_UBP(120);   fprintf('res = %.2f mm\n', out.resolution_mm);
%
%   DEPENDENCIES: k-Wave (kWaveGrid, kspaceFirstOrder2D). No project data / CT.
%
%   NOTE (>200 lines): one self-contained file so the phantom, the forward sim,
%   the UBP maths, the resolution measurement, the figure and the sensor-sweep
%   movie stay readable side by side; splitting would only add single-use files.
%
%   See also: run_single_field_simulation, kspaceFirstOrder2D

    if nargin < 1 || isempty(num_sensors), num_sensors = 60; end
    if ~isscalar(num_sensors) || num_sensors < 8 || mod(num_sensors,1) ~= 0
        error('misc_UBP:InvalidInput', 'num_sensors must be an integer >= 8.');
    end

    %% ============================ CONFIG ============================
    % Geometry (metres). Domain is square; everything is centred on (0,0).
    CONFIG.dx              = 0.25e-3;   % grid spacing (m) -> sub-mm imaging
    CONFIG.domain_mm       = 100;       % square domain side length (mm)
    CONFIG.fat_radius_mm   = 35;        % fat cylinder radius (mm)
    CONFIG.ring_radius_mm  = 28;        % sensor-ring radius (mm), < fat radius
    CONFIG.field_size_mm   = 20;        % beam field width at this plane (mm)
    CONFIG.penumbra_mm     = 0.3;       % beam edge softness, 1-sigma (mm). Kept
                                        % SHARP (< clinical ~3 mm) so the edge
                                        % probes the reconstruction, not the beam.

    % Physics (from the project's threshold_2 tissue table, define_tissue_tables).
    CONFIG.fat_c           = 1450;      % fat sound speed (m/s)   -> used by UBP
    CONFIG.fat_rho         = 920;       % fat density (kg/m^3)
    CONFIG.fat_gruneisen   = 0.7;       % fat Grueneisen (dimensionless)
    CONFIG.water_c         = 1480;      % water sound speed (m/s)
    CONFIG.water_rho       = 1000;      % water density (kg/m^3)
    CONFIG.sample_dose_Gy  = 2.0;       % peak dose of the sample beam (Gy)
    CONFIG.use_attenuation = false;     % lossless demo (see note in summary)

    % Simulation.
    CONFIG.num_sensors     = num_sensors;
    CONFIG.cfl             = 0.3;       % k-Wave stability number
    CONFIG.pml_size        = 10;        % PML thickness (voxels)
    CONFIG.use_gpu         = true;      % GPU if available, else CPU
    CONFIG.record_factor   = 1.2;       % record 1.2 x ring-diameter transit time
    CONFIG.plot            = true;

    % Movie: sweep num_sensors and animate the reconstruction.
    CONFIG.movie             = false;   % true -> run the sensor-count sweep
    CONFIG.movie_min_sensors = 10;      % smallest ring in the sweep
    CONFIG.movie_max_sensors = 360;     % largest ring (and the dense-sim ring)
    CONFIG.movie_all_counts  = false;   % false -> only divisors of the max (each
                                        %   frame exactly equidistant, prettiest).
                                        % true  -> EVERY count from min to max;
                                        %   non-divisors use the nearest dense
                                        %   sensors (<= half the dense spacing off).
    CONFIG.movie_step        = 1;       % count increment when movie_all_counts
    CONFIG.movie_fps         = 8;       % animation / export frame rate
    CONFIG.movie_save        = false;   % also write an MP4 of the sweep
    CONFIG.movie_file        = 'misc_UBP_movie.mp4';

    fprintf('\n=== misc_UBP: Universal Back-Projection PA demo ===\n');

    %% ======================= K-WAVE GRID ===========================
    dx = CONFIG.dx;
    N  = round(CONFIG.domain_mm * 1e-3 / dx);   % voxels per side
    kgrid = kWaveGrid(N, dx, N, dx);            % dim1 = x, dim2 = y

    [X, Y] = ndgrid(kgrid.x_vec, kgrid.y_vec);  % physical coords (m), centred at 0
    Rmap   = sqrt(X.^2 + Y.^2);

    fat_R  = CONFIG.fat_radius_mm  * 1e-3;
    ring_R = CONFIG.ring_radius_mm * 1e-3;

    %% =================== MEDIUM: FAT DISC IN WATER ==================
    medium = struct();
    medium.sound_speed = CONFIG.water_c   * ones(N, N);
    medium.density     = CONFIG.water_rho * ones(N, N);
    fat_mask = Rmap <= fat_R;
    medium.sound_speed(fat_mask) = CONFIG.fat_c;
    medium.density(fat_mask)     = CONFIG.fat_rho;
    if CONFIG.use_attenuation
        medium.alpha_power = 1.5;
        medium.alpha_coeff = 0.48 * fat_mask + 0.0022 * ~fat_mask;  % dB/MHz^y/cm
    end

    %% ============ SOURCE: RADIATION-BEAM CROSS-SECTION ==============
    % Flat-topped square field with an erf penumbra, centred on the ring.
    % p0 = Dose * Grueneisen * density   (Pa), the project's PA source relation.
    half_w   = 0.5 * CONFIG.field_size_mm * 1e-3;
    sigma    = max(CONFIG.penumbra_mm * 1e-3, dx);
    edge_x   = 0.5 * (erf((X + half_w) / (sqrt(2)*sigma)) - erf((X - half_w) / (sqrt(2)*sigma)));
    edge_y   = 0.5 * (erf((Y + half_w) / (sqrt(2)*sigma)) - erf((Y - half_w) / (sqrt(2)*sigma)));
    dose_Gy  = CONFIG.sample_dose_Gy * edge_x .* edge_y;   % separable square field
    dose_Gy(~fat_mask) = 0;                                % beam lives in the fat

    p0 = dose_Gy * CONFIG.fat_gruneisen * CONFIG.fat_rho;  % Pa
    source = struct();
    source.p0 = p0;

    fprintf('  Grid: %d x %d @ %.3f mm | fat R=%.0f mm | ring R=%.0f mm | field=%.0f mm\n', ...
        N, N, dx*1e3, CONFIG.fat_radius_mm, CONFIG.ring_radius_mm, CONFIG.field_size_mm);
    fprintf('  Peak p0 = %.3e Pa\n', max(p0(:)));

    %% ===================== TIME STEPPING ===========================
    dt = CONFIG.cfl * dx / max(CONFIG.water_c, CONFIG.fat_c);
    % Record long enough for any in-ring pixel to hear any sensor: up to one ring
    % diameter of travel in fat, times a safety factor.
    t_end  = CONFIG.record_factor * (2 * ring_R) / CONFIG.fat_c;
    Nt     = ceil(t_end / dt);
    kgrid.setTime(Nt, dt);
    t_axis = (0:Nt-1) * dt;
    fprintf('  dt = %.2e s, Nt = %d (T = %.2e s)\n', dt, Nt, t_end);

    if CONFIG.use_gpu
        try, gpuDevice; dataCast = 'gpuArray-single'; catch, dataCast = 'single'; end
    else
        dataCast = 'single';
    end

    %% =============== IMAGE GRID (inside the ring) ===================
    fov    = 0.98 * ring_R;
    ix_sel = find(abs(kgrid.x_vec) <= fov);
    iy_sel = find(abs(kgrid.y_vec) <= fov);
    img_x  = kgrid.x_vec(ix_sel);
    img_y  = kgrid.y_vec(iy_sel);
    p0_true = p0(ix_sel, iy_sel);

    %% ==================== MOVIE MODE (branch) =======================
    if CONFIG.movie
        out = run_sensor_movie(CONFIG, kgrid, medium, source, ring_R, ...
            t_axis, img_x, img_y, p0_true, dataCast);
        return;
    end

    %% ============ SINGLE RECONSTRUCTION (num_sensors) ===============
    fprintf('  Single reconstruction with %d sensors.\n', CONFIG.num_sensors);
    [sensor_xy, sensor_n] = ring_sensors(CONFIG.num_sensors, ring_R);

    fprintf('  Forward simulation (DataCast = %s)...\n', dataCast);
    sensorData = kspaceFirstOrder2D(kgrid, medium, source, ...
        struct('mask', sensor_xy), 'PMLInside', false, 'PMLSize', CONFIG.pml_size, ...
        'DataCast', dataCast, 'PlotSim', false);
    sensorData = double(gather(sensorData));     % [N x Nt], row i = sensor i

    fprintf('  UBP reconstruction on %d x %d image...\n', numel(img_x), numel(img_y));
    recon = universal_backprojection_2d(sensorData, t_axis, sensor_xy, sensor_n, ...
        img_x, img_y, CONFIG.fat_c, ring_R);

    R = measure_edge_resolution(recon, img_x, img_y, CONFIG);
    arc_mm = 2*pi*CONFIG.ring_radius_mm / CONFIG.num_sensors;

    fprintf('\n  --- RESOLUTION REPORT ---\n');
    fprintf('  Reconstruction resolution (edge FWHM, input removed): %.2f mm  (%.1f voxels)\n', ...
        R.res_mm, R.res_mm / (dx*1e3));
    fprintf('  Raw measured edge LSF FWHM                          : %.2f mm\n', R.res_edge_mm);
    fprintf('  Input beam edge LSF FWHM (removed)                  : %.2f mm\n', R.input_lsf_mm);
    fprintf('  Grid/Nyquist floor (2*dx)                          : %.2f mm\n', 2*dx*1e3);
    fprintf('  Sensor spacing along ring (2*pi*R/N)               : %.2f mm\n', arc_mm);
    fprintf('  Sound speed used for UBP (fat)                     : %d m/s\n', CONFIG.fat_c);
    fprintf('  -------------------------\n\n');

    if CONFIG.plot
        draw_summary(kgrid, p0, medium, sensor_xy, fat_R, ring_R, ...
            sensorData, t_axis, recon, p0_true, img_x, img_y, R, CONFIG);
    end

    out = struct('recon', recon, 'p0_true', p0_true, 'img_x', img_x, ...
        'img_y', img_y, 'sensor_xy', sensor_xy, 'resolution_mm', R.res_mm, ...
        'resolution_edge_mm', R.res_edge_mm, 'config', CONFIG);
end

% ====================================================================
function out = run_sensor_movie(CONFIG, kgrid, medium, source, ring_R, ...
        t_axis, img_x, img_y, p0_true, dataCast)
%RUN_SENSOR_MOVIE  Sweep the sensor count, animate the reconstruction.
%   One dense forward sim at movie_max_sensors; each frame reconstructs from a
%   subset of those traces (divisors of the max, or every count when
%   movie_all_counts is set). See the header MOVIE MODE note.

    Nmax = CONFIG.movie_max_sensors;
    if CONFIG.movie_all_counts
        frame_Ns = CONFIG.movie_min_sensors : CONFIG.movie_step : Nmax;   % every count
    else
        d = 1:Nmax;
        frame_Ns = d(mod(Nmax, d) == 0 & d >= CONFIG.movie_min_sensors); % exact divisors
    end
    nF = numel(frame_Ns);
    fprintf('  [MOVIE] %d frames, N = %d ... %d (all_counts = %d)\n', ...
        nF, frame_Ns(1), frame_Ns(end), CONFIG.movie_all_counts);

    % ---- one dense forward simulation ----
    [sxy_max, sn_max] = ring_sensors(Nmax, ring_R);
    fprintf('  [MOVIE] Dense forward simulation, %d sensors (DataCast = %s)...\n', Nmax, dataCast);
    dense = kspaceFirstOrder2D(kgrid, medium, source, ...
        struct('mask', sxy_max), 'PMLInside', false, 'PMLSize', CONFIG.pml_size, ...
        'DataCast', dataCast, 'PlotSim', false);
    dense = double(gather(dense));          % [Nmax x Nt]

    % ---- reconstruct every frame ----
    frames = cell(1, nF);
    res_mm = zeros(1, nF);
    for f = 1:nF
        N = frame_Ns(f);
        % Nearest dense sensor to each of N equidistant target angles. Because
        % Nmax >= N these indices are always distinct; when N divides Nmax they
        % land exactly on an equidistant subset (identical to simulating only N).
        idx = round((0:N-1) * Nmax / N) + 1;
        recon = universal_backprojection_2d(dense(idx, :), t_axis, ...
            sxy_max(:, idx), sn_max(:, idx), img_x, img_y, CONFIG.fat_c, ring_R);
        frames{f} = recon;
        R = measure_edge_resolution(recon, img_x, img_y, CONFIG);
        res_mm(f) = R.res_mm;
        if nF <= 40 || mod(f, ceil(nF/20)) == 0 || f == nF
            fprintf('    [%3d/%3d] N = %3d  ->  resolution %.2f mm\n', f, nF, N, res_mm(f));
        end
    end

    % ---- interactive slider figure ----
    cmax = max(cellfun(@(r) max(r(:)), frames));
    fig  = figure('Name', 'misc_UBP movie', 'Color', 'w', 'Position', [80 80 1120 620]);

    axR  = axes('Parent', fig, 'Position', [0.06 0.24 0.50 0.70]);
    hImg = imagesc(axR, img_y*1e3, img_x*1e3, frames{1});
    axis(axR, 'image'); set(axR, 'YDir', 'normal', 'CLim', [0 cmax]);
    colormap(axR, 'hot'); colorbar(axR);
    xlabel(axR, 'Y (mm)'); ylabel(axR, 'X (mm)'); tR = title(axR, '');

    axC = axes('Parent', fig, 'Position', [0.66 0.24 0.30 0.70]);
    plot(axC, frame_Ns, res_mm, '-o', 'Color', [0.3 0.3 0.3]); hold(axC, 'on');
    hMark = plot(axC, frame_Ns(1), res_mm(1), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 9);
    xlabel(axC, 'number of sensors'); ylabel(axC, 'resolution (mm)');
    title(axC, 'resolution vs sensors'); grid(axC, 'on');

    sld = uicontrol('Parent', fig, 'Style', 'slider', 'Units', 'normalized', ...
        'Position', [0.06 0.08 0.50 0.05], 'Min', 1, 'Max', nF, 'Value', 1, ...
        'SliderStep', [1 1] / max(nF-1, 1));
    txt = uicontrol('Parent', fig, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.06 0.02 0.50 0.04], 'BackgroundColor', 'w', ...
        'String', 'drag the slider to sweep num\_sensors');
    set(sld, 'Callback', @(s,~) show_frame(round(get(s,'Value'))));
    addlistener(sld, 'ContinuousValueChange', @(s,~) show_frame(round(get(s,'Value'))));

    show_frame(1);

    % ---- optionally save a complete MP4, then a short on-screen preview ----
    if CONFIG.movie_save
        vw = VideoWriter(CONFIG.movie_file, 'MPEG-4');
        vw.FrameRate = CONFIG.movie_fps; open(vw);
        for f = 1:nF                       % every frame, for a smooth movie file
            set(sld, 'Value', f); show_frame(f); drawnow; writeVideo(vw, getframe(fig));
        end
        close(vw); fprintf('  [MOVIE] saved %s\n', CONFIG.movie_file);
    end
    stride = max(1, round(nF / 60));       % preview <= ~60 frames so it stays brief
    for f = 1:stride:nF
        set(sld, 'Value', f); show_frame(f); drawnow; pause(1/CONFIG.movie_fps);
    end
    set(sld, 'Value', nF); show_frame(nF); % leave the slider live at the finest ring

    out = struct('frames', {frames}, 'num_sensors', frame_Ns, ...
        'resolution_mm', res_mm, 'img_x', img_x, 'img_y', img_y, ...
        'p0_true', p0_true, 'config', CONFIG);

    function show_frame(k)
        k = min(max(k, 1), nF);
        set(hImg, 'CData', frames{k});
        set(hMark, 'XData', frame_Ns(k), 'YData', res_mm(k));
        set(tR, 'String', sprintf('UBP recon:  N = %d sensors,  resolution = %.2f mm', ...
            frame_Ns(k), res_mm(k)));
        set(txt, 'String', sprintf('frame %d / %d   (N = %d sensors)', k, nF, frame_Ns(k)));
    end
end

% ====================================================================
function [sxy, sn] = ring_sensors(N, ring_R)
%RING_SENSORS  N equally-spaced ring points (2 x N) and their inward normals.
    theta = (0:N-1) * (2*pi / N);
    sxy = [ring_R * cos(theta); ring_R * sin(theta)];
    sn  = [-cos(theta);         -sin(theta)];
end

% ====================================================================
function recon = universal_backprojection_2d(p, t_axis, s_xy, s_n, img_x, img_y, c, ring_R)
%UNIVERSAL_BACKPROJECTION_2D  UBP of Xu & Wang (2005) in the detector plane.
%
%   For every sensor i the filtered trace  b_i = 2*p_i - 2*t*dp_i/dt  is
%   back-projected onto the circular arcs t = |r - r0_i|/c and summed with the
%   2D solid-angle weight (cos(theta)/|r - r0|), theta = angle between the
%   sensor's inward normal and the direction to the image point. The 1/Omega0
%   normalisation (Omega0 = 2*pi) times the per-sensor arc length (2*pi*R/N)
%   collapses to the prefactor R/N. Absolute scale is otherwise not critical.

    N  = size(p, 1);
    dt = t_axis(2) - t_axis(1);

    % Temporal derivative and the UBP filter b = 2*(p - t*dp/dt).
    dpdt              = zeros(size(p));
    dpdt(:, 2:end-1)  = (p(:, 3:end) - p(:, 1:end-2)) / (2*dt);
    dpdt(:, 1)        = (p(:, 2)   - p(:, 1))     / dt;
    dpdt(:, end)      = (p(:, end) - p(:, end-1)) / dt;
    b = 2*p - 2 * t_axis .* dpdt;    % t_axis is 1 x Nt, broadcasts over sensors

    [IX, IY] = ndgrid(img_x, img_y);
    recon    = zeros(size(IX));
    for i = 1:N
        dxr = IX - s_xy(1, i);
        dyr = IY - s_xy(2, i);
        d   = sqrt(dxr.^2 + dyr.^2);
        d   = max(d, eps);
        cos_th = (s_n(1, i) .* dxr + s_n(2, i) .* dyr) ./ d;   % normal . unit(r-r0)
        tq  = d / c;
        bval = interp1(t_axis, b(i, :), tq(:), 'linear', 0);   % 0 outside record
        recon = recon + (cos_th ./ d) .* reshape(bval, size(IX));
    end
    recon = recon * (ring_R / N);            % (1/Omega0) * arc length per sensor

    Rimg = sqrt(IX.^2 + IY.^2);
    recon(Rimg > max(img_x(end), img_y(end))) = 0;   % blank the corners outside FOV
end

% ====================================================================
function R = measure_edge_resolution(recon, img_x, img_y, CONFIG)
%MEASURE_EDGE_RESOLUTION  Reconstruction resolution from a reconstructed edge.
%   FWHM of the line-spread function (|d/dx| of the central row) at the +x field
%   edge, with the known input beam edge removed in quadrature.
    [~, iy0] = min(abs(img_y));
    xmm  = img_x(:) * 1e3;
    prof = recon(:, iy0);
    prof = prof / max(prof);                     % normalise
    lsf  = abs(gradient(prof, xmm));             % edge -> line-spread function

    win = abs(xmm - 0.5*CONFIG.field_size_mm) <= 4;   % ~4 mm around +x edge
    xw  = xmm(win); lw = lsf(win);
    res_edge = local_fwhm(xw, lw);

    input_lsf = 2.3548 * CONFIG.penumbra_mm;     % FWHM of the erf edge's LSF
    res       = sqrt(max(res_edge^2 - input_lsf^2, 0));

    R = struct('res_mm', res, 'res_edge_mm', res_edge, 'input_lsf_mm', input_lsf, ...
        'iy0', iy0, 'xmm', xmm, 'prof', prof, 'xw', xw, 'lw', lw);
end

% ====================================================================
function w = local_fwhm(x, y)
%LOCAL_FWHM  Full width at half maximum of a single positive peak y(x).
    x = x(:); y = y(:);
    [ymax, imax] = max(y);
    if ymax <= 0, w = NaN; return; end
    half = 0.5 * ymax;
    iL = find(y(1:imax)   <= half, 1, 'last');
    iR = find(y(imax:end) <= half, 1, 'first') + imax - 1;
    if isempty(iL) || isempty(iR), w = NaN; return; end
    % Explicit linear interpolation of each half-max crossing (robust to slope
    % direction, unlike interp1 on a decreasing sample vector).
    xL = x(iL)   + (half - y(iL))   * (x(iL+1) - x(iL))   / (y(iL+1) - y(iL));
    xR = x(iR-1) + (half - y(iR-1)) * (x(iR)   - x(iR-1)) / (y(iR)   - y(iR-1));
    w  = xR - xL;
end

% ====================================================================
function draw_summary(kgrid, p0, medium, s_xy, fat_R, ring_R, sensorData, ...
        t_axis, recon, p0_true, img_x, img_y, R, CONFIG)
%DRAW_SUMMARY  Six-panel overview (2 rows x 3 cols) for the single-N run.
    xg = kgrid.x_vec * 1e3;  yg = kgrid.y_vec * 1e3;   % mm
    th = linspace(0, 2*pi, 200);
    figure('Name', 'misc_UBP', 'Color', 'w', 'Position', [80 80 1300 760]);

    % (1) Source p0 with sensor ring.
    subplot(2,3,1);
    imagesc(yg, xg, p0); axis image; set(gca,'YDir','normal'); colormap(gca,'hot');
    hold on; plot(s_xy(2,:)*1e3, s_xy(1,:)*1e3, 'c.', 'MarkerSize', 8);
    title('Beam cross-section p_0'); xlabel('Y (mm)'); ylabel('X (mm)'); colorbar;

    % (2) Medium (sound speed) with fat/ring geometry.
    subplot(2,3,2);
    imagesc(yg, xg, medium.sound_speed); axis image; set(gca,'YDir','normal');
    colormap(gca,'parula'); hold on;
    plot(fat_R*1e3*cos(th), fat_R*1e3*sin(th), 'w--', 'LineWidth', 1);
    plot(s_xy(2,:)*1e3, s_xy(1,:)*1e3, 'r.', 'MarkerSize', 8);
    title('Medium: fat disc in water'); xlabel('Y (mm)'); ylabel('X (mm)'); colorbar;

    % (3) Sensor traces (sinogram).
    subplot(2,3,3);
    imagesc(t_axis*1e6, 1:size(sensorData,1), sensorData);
    colormap(gca,'gray'); title('Sensor traces'); xlabel('t (\mus)'); ylabel('sensor #'); colorbar;

    % (4) UBP reconstruction.
    subplot(2,3,4);
    imagesc(img_y*1e3, img_x*1e3, recon); axis image; set(gca,'YDir','normal');
    colormap(gca,'hot'); title('UBP reconstruction'); xlabel('Y (mm)'); ylabel('X (mm)'); colorbar;

    % (5) Central-row profile: truth vs recon (normalised).
    subplot(2,3,5);
    pt = p0_true(:, R.iy0); pt = pt / max(pt);
    plot(R.xmm, pt, 'k-', 'LineWidth', 1.5); hold on;
    plot(R.xmm, R.prof, 'r-', 'LineWidth', 1.5);
    xlabel('X (mm)'); ylabel('normalised'); legend('truth','recon','Location','south');
    title('Central profile (y \approx 0)'); grid on; xlim([min(R.xmm) max(R.xmm)]);

    % (6) Line-spread function at the +x edge with FWHM.
    subplot(2,3,6);
    plot(R.xw, R.lw, 'b-', 'LineWidth', 1.5); hold on;
    yl = ylim; plot([1 1]*0.5*CONFIG.field_size_mm, yl, 'k:');
    xlabel('X (mm)'); ylabel('|d(recon)/dx|');
    title(sprintf('Edge LSF: meas %.2f mm, recon %.2f mm', R.res_edge_mm, R.res_mm)); grid on;
end
