function out = misc_UBP(num_sensors)
%MISC_UBP  Universal back-projection (UBP) photoacoustic demo in a fat phantom.
%
%   out = misc_UBP()            % 60 sensors (default)
%   out = misc_UBP(num_sensors) % choose the number of ring sensors
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
%   INPUTS:
%       num_sensors - (optional) number of equally-spaced point sensors on the
%                     ring. Positive integer, default 60.
%
%   OUTPUTS:
%       out - struct with:
%           .recon           - reconstructed pressure image (Pa), [nx_img x ny_img]
%           .p0_true         - ground-truth initial pressure on the same grid
%           .img_x, .img_y   - image axis coordinates (m)
%           .sensor_xy       - [2 x N] sensor coordinates (m)
%           .resolution_mm   - measured reconstruction resolution (LSF FWHM, mm)
%           .config          - the CONFIG struct actually used
%
%   ALGORITHM:
%       1. Build fat-disc-in-water medium and the beam cross-section p0.
%       2. Place N point sensors on a circle inside the fat, beam through centre.
%       3. k-Wave forward sim -> pressure trace at every sensor.
%       4. UBP: filter each trace (b = 2p - 2t dp/dt) and back-project onto the
%          circular arcs t = |r - r0|/c, weighted by solid angle.
%       5. Resolution = FWHM of the line-spread function at a reconstructed edge.
%
%   EXAMPLE:
%       out = misc_UBP(120);   fprintf('res = %.2f mm\n', out.resolution_mm);
%
%   DEPENDENCIES: k-Wave (kWaveGrid, kspaceFirstOrder2D). No project data / CT.
%
%   NOTE (>200 lines): one self-contained file so the phantom, the forward sim,
%   the UBP maths, the resolution measurement and the figure stay readable side
%   by side; splitting would only add single-use files.
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

    fprintf('\n=== misc_UBP: Universal Back-Projection PA demo (%d sensors) ===\n', ...
        CONFIG.num_sensors);

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

    %% =================== SENSORS: RING OF POINTS ====================
    % Equally spaced points on a circle of radius ring_R, inside the fat, with
    % inward normals (pointing at the ring centre) for the UBP solid-angle term.
    theta   = (0:CONFIG.num_sensors-1) * (2*pi / CONFIG.num_sensors);
    sensor_xy = [ring_R * cos(theta); ring_R * sin(theta)];   % [2 x N], centred
    sensor_n  = [-cos(theta);         -sin(theta)];           % inward unit normals

    sensor = struct();
    sensor.mask = sensor_xy;    % Cartesian sensor mask (k-Wave interpolates points)

    %% ===================== TIME STEPPING ===========================
    dt = CONFIG.cfl * dx / max(CONFIG.water_c, CONFIG.fat_c);
    % Record long enough for any in-ring pixel to hear any sensor: up to one ring
    % diameter of travel in fat, times a safety factor.
    t_end = CONFIG.record_factor * (2 * ring_R) / CONFIG.fat_c;
    Nt    = ceil(t_end / dt);
    kgrid.setTime(Nt, dt);
    t_axis = (0:Nt-1) * dt;
    fprintf('  dt = %.2e s, Nt = %d (T = %.2e s)\n', dt, Nt, t_end);

    %% ===================== FORWARD SIMULATION ======================
    if CONFIG.use_gpu
        try, gpuDevice; dataCast = 'gpuArray-single'; catch, dataCast = 'single'; end
    else
        dataCast = 'single';
    end
    fprintf('  Forward simulation (DataCast = %s)...\n', dataCast);

    sensorData = kspaceFirstOrder2D(kgrid, medium, source, sensor, ...
        'PMLInside', false, 'PMLSize', CONFIG.pml_size, ...
        'DataCast', dataCast, 'PlotSim', false);
    sensorData = double(gather(sensorData));     % [N x Nt], row i = sensor i

    %% ================ UNIVERSAL BACK-PROJECTION =====================
    % Reconstruct just inside the ring (recon is only valid there).
    fov     = 0.98 * ring_R;
    ix_sel  = find(abs(kgrid.x_vec) <= fov);
    iy_sel  = find(abs(kgrid.y_vec) <= fov);
    img_x   = kgrid.x_vec(ix_sel);
    img_y   = kgrid.y_vec(iy_sel);

    fprintf('  UBP reconstruction on %d x %d image...\n', numel(img_x), numel(img_y));
    recon = universal_backprojection_2d(sensorData, t_axis, sensor_xy, sensor_n, ...
        img_x, img_y, CONFIG.fat_c, ring_R);

    p0_true = p0(ix_sel, iy_sel);   % ground truth on the image grid

    %% ================== RESOLUTION MEASUREMENT ======================
    % Line-spread function at the +x field edge on the central row (y ~ 0).
    [~, iy0]  = min(abs(img_y));
    xmm       = img_x(:) * 1e3;
    prof      = recon(:, iy0);
    prof      = prof / max(prof);                 % normalise
    lsf       = abs(gradient(prof, xmm));         % edge -> line-spread function

    edge_win  = abs(xmm - 0.5*CONFIG.field_size_mm) <= 4;   % ~4 mm around +x edge
    xw        = xmm(edge_win);
    lw        = lsf(edge_win);
    res_edge_mm = local_fwhm(xw, lw);              % measured edge LSF FWHM (mm)

    % Isolate the reconstruction's own resolution by removing the known input
    % edge softness in quadrature (Gaussian LSFs: FWHM_meas^2 = FWHM_in^2 + FWHM_sys^2).
    input_lsf_mm = 2.3548 * CONFIG.penumbra_mm;    % FWHM of the erf edge's LSF
    res_mm       = sqrt(max(res_edge_mm^2 - input_lsf_mm^2, 0));

    arc_mm    = 2*pi*CONFIG.ring_radius_mm / CONFIG.num_sensors;   % sensor spacing on ring

    fprintf('\n  --- RESOLUTION REPORT ---\n');
    fprintf('  Reconstruction resolution (edge FWHM, input removed): %.2f mm  (%.1f voxels)\n', ...
        res_mm, res_mm / (dx*1e3));
    fprintf('  Raw measured edge LSF FWHM                          : %.2f mm\n', res_edge_mm);
    fprintf('  Input beam edge LSF FWHM (removed)                  : %.2f mm\n', input_lsf_mm);
    fprintf('  Grid/Nyquist floor (2*dx)                          : %.2f mm\n', 2*dx*1e3);
    fprintf('  Sensor spacing along ring (2*pi*R/N)               : %.2f mm\n', arc_mm);
    fprintf('  Sound speed used for UBP (fat)                     : %d m/s\n', CONFIG.fat_c);
    fprintf('  -------------------------\n\n');

    %% ========================= FIGURE ==============================
    if CONFIG.plot
        draw_summary(kgrid, p0, medium, sensor_xy, fat_R, ring_R, ...
            sensorData, t_axis, recon, p0_true, img_x, img_y, iy0, ...
            xmm, prof, xw, lw, res_mm, res_edge_mm, CONFIG);
    end

    %% ========================= OUTPUT ==============================
    out = struct('recon', recon, 'p0_true', p0_true, 'img_x', img_x, ...
        'img_y', img_y, 'sensor_xy', sensor_xy, 'resolution_mm', res_mm, ...
        'resolution_edge_mm', res_edge_mm, 'config', CONFIG);
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
        t_axis, recon, p0_true, img_x, img_y, iy0, xmm, prof, xw, lw, res_mm, res_edge_mm, CONFIG)
%DRAW_SUMMARY  Six-panel overview (2 rows x 3 cols).
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
    pt = p0_true(:, iy0); pt = pt / max(pt);
    plot(xmm, pt, 'k-', 'LineWidth', 1.5); hold on;
    plot(xmm, prof, 'r-', 'LineWidth', 1.5);
    xlabel('X (mm)'); ylabel('normalised'); legend('truth','recon','Location','south');
    title('Central profile (y \approx 0)'); grid on; xlim([min(xmm) max(xmm)]);

    % (6) Line-spread function at the +x edge with FWHM.
    subplot(2,3,6);
    plot(xw, lw, 'b-', 'LineWidth', 1.5); hold on;
    yl = ylim; plot([1 1]*0.5*CONFIG.field_size_mm, yl, 'k:');
    xlabel('X (mm)'); ylabel('|d(recon)/dx|');
    title(sprintf('Edge LSF: meas %.2f mm, recon %.2f mm', res_edge_mm, res_mm)); grid on;
end
