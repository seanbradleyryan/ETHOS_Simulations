% diagnostic_blind_fidelity.m
%
% PURPOSE:
%   Quantify — without ever looking at a truth dose (hence "blind") — whether
%   a MAJOR change in the acoustic medium and/or sensor coupling between CT_1
%   (reference) and CT_3 (adapted "measurement") anatomy explains why the
%   blind time-reversal reconstruction collapses on CT_3. It answers the
%   advisor's question ("would a small medium change really break TR?") by
%   turning the 32x32 array into a handful of single, readable images and one
%   self-calibrated scatter, instead of eyeballing 1024 channels or 3 slices.
%
%   Three diagnostics, all on the SAME 32x32 element grid:
%
%   [1] PER-CHANNEL SIGNAL CHANGE (needs the forward operator).
%       Hold p0 and the sensor FIXED (both taken from CT_1) and propagate the
%       identical p0 through the CT_1 medium and the CT_3 medium. Difference
%       the raw sensor traces element by element:
%         - Pearson (zero-lag) correlation map: a dead/blocked element shows as
%           a CONTIGUOUS dark blob (physical occlusion); salt-and-pepper decorr
%           is numerics, not anatomy.
%         - Energy-ratio map: amplitude loss even when the waveform shape holds.
%
%   [2] IMPEDANCE / REFLECTION RAY INTEGRAL (geometry only, no k-Wave).
%       For each element, cast rays from the high-dose (beam) source set to the
%       element through the CT_1 and CT_3 media and accumulate a reflection
%       proxy R = sum_i ((Z_{i+1}-Z_i)/(Z_{i+1}+Z_i))^2 (Z = rho*c). The map is
%       the CT_3 - CT_1 difference. It is SELF-CALIBRATED against [1]: a scatter
%       of dR vs (1 - corr) per element gives the reflection change at which
%       channels actually die — no need to guess an absolute impedance cutoff.
%
%   [4] STANDOFF / COUPLING + SENSOR-EYE PROJECTION (geometry only, no k-Wave).
%       Per element, distance along its inward look-direction (toward iso) to
%       the first body voxel, for CT_1 and CT_3 (32x32 maps + difference). An
%       air gap (large or NaN standoff) is an unambiguous, catastrophic coupling
%       failure. A per-element sensor-eye DRR (integrated HU along the same ray)
%       shows structure that moved into an element's line of sight. Coincidence
%       of these maps with the [1] dead-channel blob localizes the cause
%       (coupling gap vs anatomy in path vs bulk impedance).
%
%       ([3] "reconstruct CT_3 data on the CT_3 medium" is intentionally omitted:
%        the matched case is already known to converge, better than CT_1.)
%
% HOW IT GETS THE FORWARD DATA:
%   Diagnostic [1] needs two raw forward sinograms with the SAME p0 and SAME
%   sensor. That is produced by the one simulation engine, not re-derived here:
%   run_standalone_field is called on the blind path (forward = CBCT1, recon =
%   CBCT3) with config.export_forward_sinograms / forward_sinograms_only set, so
%   run_single_field_simulation places the sensor, builds p0 from CT_1, and runs
%   the forward operator through BOTH media, returning the raw sinogram pair and
%   sensor_info (element positions + voxel->element map). No time reversal runs.
%
% USAGE:
%   Edit the CONFIG block (patient / session / working_dir / dose file), then
%   run (F5). HIPAA: execute on the remote device, not locally.
%
% OUTPUTS:
%   - 3 figures (one per diagnostic) + a saved results .mat (all maps + metrics).
%   - Console summary: #dead channels, largest contiguous blob, max standoff
%     change, air-gap count, and the three coincidence correlations that say
%     whether the dead-channel pattern is coupling-, anatomy-, or impedance-driven.
%
% DEPENDENCIES:
%   run_standalone_field, run_single_field_simulation (export hook),
%   create_acoustic_medium, get_default_config, Image Processing Toolbox
%   (bwconncomp for blob sizing).
%
% See also: run_standalone_field, run_single_field_simulation,
%           determine_sensor_mask, apply_element_averaging
%
% JUSTIFICATION (>200 lines): three independent 32x32 diagnostics, each with its
% own map/metric build and plot, plus geometry/impedance ray casting shared
% between [2] and [4]; splitting into separate files would scatter one coherent
% report across the repo.

clear; close all; clc;

% Script lives one level below the repo root; add utils/ and pipeline/.
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repoRoot, 'utils')));
addpath(genpath(fullfile(repoRoot, 'pipeline')));

%% ============================ CONFIGURATION ============================ %%

CONFIG = get_default_config();

% --- Data selection (EDIT THESE) ---
CONFIG.working_dir   = '/mnt/weka/home/80030361/ETHOS_Simulations';
CONFIG.patient_id    = '1194203';
CONFIG.session       = 'Session_1';
CONFIG.dose_filename = 'dose_1194203_Session_1_reference_CT_1_B15_112.mat';

% Forward geometry = CT_1 (reference); reconstruction geometry = CT_3 (adapted).
% The forward sinogram PAIR is CT_1-medium vs CT_3-medium with this p0/sensor.
CONFIG.cbct_filename       = 'CBCT1_resampled.mat';
CONFIG.recon_cbct_filename = 'CBCT3_resampled.mat';
CONFIG.blind_recon         = true;

% --- Diagnostic engine flags: return the raw forward sinogram pair, skip TR ---
CONFIG.export_forward_sinograms = true;
CONFIG.forward_sinograms_only   = true;   % return right after the two forward runs
CONFIG.convolution_kernel       = 0;      % RAW sinograms (no pulse/noise/deconv)
CONFIG.plot_results             = false;
CONFIG.save_results             = false;

% --- Diagnostic parameters ---
DIAG = struct();
DIAG.src_dose_frac   = 0.50;   % [2] source set = dose >= frac * max(dose)  (beam-shaped)
DIAG.max_src_points  = 64;     % [2] cap on source points per element (thinned isodose)
DIAG.ray_nsamp       = 200;    % samples per ray for [2] reflection and [4] DRR
DIAG.standoff_step_mm = 0.5;   % [4] march step toward iso
DIAG.standoff_max_mm  = 120;   % [4] give up (air gap / no contact) beyond this
DIAG.corr_dead_thresh = 0.50;  % [1] Pearson below this => channel counted "dead"
DIAG.energy_lo        = 0.50;  % [1] energy ratio (CT3/CT1) outside [lo, 1/lo] => flagged
DIAG.energy_floor_frac = 0.05; % [1] skip elements below this fraction of the peak channel energy
DIAG.save_dir         = fullfile(CONFIG.working_dir, 'AnalysisResults', ...
                                 CONFIG.patient_id, CONFIG.session);

bar = repmat('=', 1, 72);
fprintf('\n%s\n  BLIND FIDELITY DIAGNOSTIC (CT_1 vs CT_3)\n%s\n', bar, bar);
fprintf('Patient/Session : %s / %s\n', CONFIG.patient_id, CONFIG.session);
fprintf('Dose file       : %s\n', CONFIG.dose_filename);
fprintf('Forward medium  : %s\n', CONFIG.cbct_filename);
fprintf('Recon  medium   : %s\n', CONFIG.recon_cbct_filename);
fprintf('%s\n\n', bar);

%% =================== STEP 1: FORWARD SINOGRAM PAIR ===================== %%
%  One engine call; the export hook runs the forward operator on both media
%  with the SAME p0 and SAME sensor, then returns (no time reversal).

fprintf('[1/5] Running forward sinogram pair (CT_1 & CT_3 media)...\n');
out = run_standalone_field(CONFIG);
sr  = out.sim_results;

assert(isfield(sr, 'sinogram_fwd') && isfield(sr, 'sinogram_alt'), ...
    'Engine did not return the sinogram pair. Check export_forward_sinograms hook.');
sino1 = sr.sinogram_fwd;      % [Nvox x Nt] CT_1 medium
sino3 = sr.sinogram_alt;      % [Nvox x Nt] CT_3 medium
sinfo = sr.sensor_info;
epos  = sinfo.element_positions_mm;             % [N_total x 3] mm (DICOM x,y,z)
N     = sinfo.elements_per_side;                % elements per side (may be < 32)
Ntot  = size(epos, 1);
velem = double(sinfo.voxel_element_idx(:));     % [Nvox x 1] element index per sensor voxel
fprintf('      Sinograms: [%d voxels x %d samples]; array %dx%d = %d elements\n', ...
    size(sino1,1), size(sino1,2), N, N, Ntot);

%% =================== STEP 2: GEOMETRY & MEDIA ========================= %%
%  Load CT_1 / CT_3 geometry directly (they share the dose grid). Element
%  positions are ABSOLUTE mm, so [2]/[4] map them onto the native CBCT grid.

fprintf('[2/5] Loading CT_1 / CT_3 geometry and building media...\n');
processed_dir = fullfile(CONFIG.working_dir, 'RayStationFiles', ...
    CONFIG.patient_id, CONFIG.session, 'processed');
sct1 = load_cbct_geo(fullfile(processed_dir, CONFIG.cbct_filename));
sct3 = load_cbct_geo(fullfile(processed_dir, CONFIG.recon_cbct_filename));

origin  = sct1.origin(:)';
spacing = sct1.spacing(:)';
assert(isequal(size(sct1.cubeHU), size(sct3.cubeHU)), ...
    'CT_1 and CT_3 grids differ; they must share the dose grid.');

[Z1, ~]  = build_impedance(sct1, CONFIG);   % Z = rho*c on the CT_1 grid
[Z3, ~]  = build_impedance(sct3, CONFIG);   % Z = rho*c on the CT_3 grid
body1 = logical(sct1.bodyMask);
body3 = logical(sct3.bodyMask);
HU1   = sct1.cubeHU;
HU3   = sct3.cubeHU;

doseGrid = out.doseGrid;                    % CT_1-body-masked dose, native grid
iso_mm   = resolve_iso_mm(sinfo, doseGrid, origin, spacing);
fprintf('      Isocenter (look target): [%.1f %.1f %.1f] mm\n', iso_mm);

%% =================== DIAGNOSTIC [1]: PER-CHANNEL ====================== %%
%  Average sensor-voxel traces per element (authoritative voxel->element map),
%  then per element compute Pearson correlation and energy ratio CT_3 vs CT_1.

fprintf('[3/5] Diagnostic [1]: per-channel signal change...\n');
if numel(velem) ~= size(sino1, 1)
    warning('diagnostic_blind_fidelity:VoxelCountMismatch', ...
        'voxel_element_idx (%d) ~= sinogram rows (%d); channel map may be misaligned.', ...
        numel(velem), size(sino1, 1));
end
corrMap = nan(Ntot, 1);      % Pearson(CT1 trace, CT3 trace)
enerMap = nan(Ntot, 1);      % ||CT3|| / ||CT1||
E1      = nan(Ntot, 1);      % CT1 trace energy (for the low-signal floor)
for e = 1:Ntot
    rows = find(velem == e);
    if isempty(rows); continue; end
    t1 = mean(sino1(rows, :), 1);
    t3 = mean(sino3(rows, :), 1);
    E1(e) = norm(t1);
    corrMap(e) = pearson(t1, t3);
    if E1(e) > 0
        enerMap(e) = norm(t3) / E1(e);
    end
end
% Low-signal elements can't be meaningfully "blocked"; excluding them keeps the
% dead-channel map (and the coincidence stats) from being polluted by noise.
if any(~isnan(E1))
    sig_floor = DIAG.energy_floor_frac * max(E1(~isnan(E1)));
    lowsig = E1 < sig_floor;
    corrMap(lowsig) = NaN;
    enerMap(lowsig) = NaN;
end
decorr = 1 - corrMap;                                    % 0 = identical, up to 2
deadMask = (corrMap < DIAG.corr_dead_thresh) | ...
           (enerMap < DIAG.energy_lo) | (enerMap > 1/DIAG.energy_lo);
deadMask(isnan(corrMap)) = false;                        % aliased elements: no signal

corrImg = reshape(corrMap, [N, N]);   % rows = ez, cols = ex
enerImg = reshape(enerMap, [N, N]);
deadImg = reshape(double(deadMask), [N, N]);

%% =================== DIAGNOSTIC [2]: REFLECTION RAY ==================== %%
%  Cumulative reflection proxy along source->element rays, CT_1 and CT_3.

fprintf('[4/5] Diagnostic [2]: impedance / reflection ray integral...\n');
src_mm = beam_source_points(doseGrid, origin, spacing, ...
    DIAG.src_dose_frac, DIAG.max_src_points);
fprintf('      Source points (beam >= %.0f%% max dose): %d\n', ...
    100*DIAG.src_dose_frac, size(src_mm,1));

refl1 = nan(Ntot, 1);
refl3 = nan(Ntot, 1);
for e = 1:Ntot
    pe = epos(e, :);
    r1 = 0; r3 = 0;
    for s = 1:size(src_mm, 1)
        r1 = r1 + reflection_along_ray(Z1, origin, spacing, src_mm(s,:), pe, DIAG.ray_nsamp);
        r3 = r3 + reflection_along_ray(Z3, origin, spacing, src_mm(s,:), pe, DIAG.ray_nsamp);
    end
    refl1(e) = r1 / size(src_mm, 1);
    refl3(e) = r3 / size(src_mm, 1);
end
dRefl = refl3 - refl1;
refl1Img = reshape(refl1, [N, N]);
refl3Img = reshape(refl3, [N, N]);
dReflImg = reshape(dRefl, [N, N]);

%% =================== DIAGNOSTIC [4]: STANDOFF + SENSOR-EYE ============ %%
%  Per element: distance to body along the inward look-direction (toward iso),
%  and integrated HU along that ray (sensor-eye DRR), for CT_1 and CT_3.

fprintf('[5/5] Diagnostic [4]: standoff / coupling + sensor-eye DRR...\n');
stand1 = nan(Ntot, 1);  stand3 = nan(Ntot, 1);
drr1   = nan(Ntot, 1);  drr3   = nan(Ntot, 1);
for e = 1:Ntot
    pe = epos(e, :);
    d  = iso_mm - pe;
    nd = norm(d);
    if nd < eps; continue; end
    look = d / nd;                              % inward: element -> iso
    stand1(e) = standoff_to_body(body1, origin, spacing, pe, look, ...
        DIAG.standoff_step_mm, DIAG.standoff_max_mm);
    stand3(e) = standoff_to_body(body3, origin, spacing, pe, look, ...
        DIAG.standoff_step_mm, DIAG.standoff_max_mm);
    drr1(e) = drr_along_ray(HU1, origin, spacing, pe, iso_mm, DIAG.ray_nsamp);
    drr3(e) = drr_along_ray(HU3, origin, spacing, pe, iso_mm, DIAG.ray_nsamp);
end
dStand = stand3 - stand1;
dDRR   = drr3   - drr1;
airGap = isnan(stand1) | isnan(stand3);         % element never reaches body

stand1Img = reshape(stand1, [N, N]);
stand3Img = reshape(stand3, [N, N]);
dStandImg = reshape(dStand, [N, N]);
drr1Img   = reshape(drr1,   [N, N]);
drr3Img   = reshape(drr3,   [N, N]);
dDRRImg   = reshape(dDRR,   [N, N]);

%% =================== COINCIDENCE / SELF-CALIBRATION =================== %%
%  Do the geometry changes explain the dead channels? Correlate each candidate
%  cause with the per-channel decorrelation from [1], over elements that have a
%  valid signal AND valid geometry.

valid = ~isnan(decorr);
cReflVsDecorr  = pearson_masked(dRefl,  decorr, valid);
cStandVsDecorr = pearson_masked(dStand, decorr, valid);
cDRRVsDecorr   = pearson_masked(dDRR,   decorr, valid);

n_dead = nnz(deadMask);
n_air  = nnz(airGap);
blob   = largest_blob(deadImg > 0.5);
maxStand = max(abs(dStand(valid)));

%% =========================== PLOTS ==================================== %%

% Figure 1 — Diagnostic [1]
figure('Name', '[1] Per-channel signal change', 'Color', 'w');
subplot(1,3,1); show_map(corrImg, 'Pearson corr (CT1 vs CT3)', [0 1]);      colormap(gca, 'parula');
subplot(1,3,2); show_map(enerImg, 'Energy ratio CT3/CT1', [0 2]);          colormap(gca, 'parula');
subplot(1,3,3); show_map(deadImg, sprintf('Dead channels (%d)', n_dead), [0 1]); colormap(gca, 'hot');
sgtitle('Diagnostic [1]: contiguous dark blob = physical occlusion; scatter = numerics');

% Figure 2 — Diagnostic [2] + self-calibration
figure('Name', '[2] Reflection ray integral', 'Color', 'w');
subplot(2,2,1); show_map(refl1Img, 'Reflection proxy CT1', []);
subplot(2,2,2); show_map(refl3Img, 'Reflection proxy CT3', []);
subplot(2,2,3); show_map(dReflImg, '\Delta Reflection (CT3-CT1)', []); colormap(gca, diverging_map());
subplot(2,2,4);
scatter(dRefl(valid), decorr(valid), 14, 'filled'); grid on;
xlabel('\Delta reflection proxy (CT3 - CT1)'); ylabel('1 - Pearson corr');
title(sprintf('Self-calibration: r = %.2f', cReflVsDecorr));
sgtitle('Diagnostic [2]: which reflection change kills channels (read off the scatter)');

% Figure 3 — Diagnostic [4]
figure('Name', '[4] Standoff + sensor-eye', 'Color', 'w');
subplot(2,3,1); show_map(stand1Img, 'Standoff CT1 (mm)', []);
subplot(2,3,2); show_map(stand3Img, 'Standoff CT3 (mm)', []);
subplot(2,3,3); show_map(dStandImg, '\Delta Standoff (mm)', []); colormap(gca, diverging_map());
subplot(2,3,4); show_map(drr1Img, 'Sensor-eye DRR CT1', []);
subplot(2,3,5); show_map(drr3Img, 'Sensor-eye DRR CT3', []);
subplot(2,3,6); show_map(dDRRImg, '\Delta DRR (CT3-CT1)', []); colormap(gca, diverging_map());
sgtitle(sprintf(['Diagnostic [4]: air-gap channels = %d | ' ...
    'coincidence r(\\Deltastandoff,decorr)=%.2f, r(\\DeltaDRR,decorr)=%.2f'], ...
    n_air, cStandVsDecorr, cDRRVsDecorr));

%% =========================== SAVE + SUMMARY ========================== %%

if ~isfolder(DIAG.save_dir); mkdir(DIAG.save_dir); end
[~, dose_tag] = fileparts(CONFIG.dose_filename);
results = struct();
results.config       = CONFIG;
results.diag_params  = DIAG;
results.elements_per_side = N;
results.maps = struct('corr', corrImg, 'energy', enerImg, 'dead', deadImg, ...
    'refl_CT1', refl1Img, 'refl_CT3', refl3Img, 'd_refl', dReflImg, ...
    'standoff_CT1', stand1Img, 'standoff_CT3', stand3Img, 'd_standoff', dStandImg, ...
    'drr_CT1', drr1Img, 'drr_CT3', drr3Img, 'd_drr', dDRRImg);
results.coincidence = struct('refl_vs_decorr', cReflVsDecorr, ...
    'standoff_vs_decorr', cStandVsDecorr, 'drr_vs_decorr', cDRRVsDecorr);
results.summary = struct('n_dead', n_dead, 'n_air_gap', n_air, ...
    'largest_dead_blob', blob, 'max_abs_d_standoff_mm', maxStand);
out_mat = fullfile(DIAG.save_dir, sprintf('blind_fidelity_%s.mat', dose_tag));
save(out_mat, 'results', '-v7.3');

fprintf('\n%s\n  SUMMARY\n%s\n', bar, bar);
fprintf('  Dead channels (corr<%.2f or energy out of [%.2f,%.2f]) : %d / %d\n', ...
    DIAG.corr_dead_thresh, DIAG.energy_lo, 1/DIAG.energy_lo, n_dead, nnz(valid));
fprintf('  Largest contiguous dead blob (elements)               : %d\n', blob);
fprintf('  Air-gap / no-contact elements (NaN standoff)          : %d\n', n_air);
fprintf('  Max |standoff change| CT3 vs CT1 (mm)                 : %.1f\n', maxStand);
fprintf('  --- coincidence with per-channel decorrelation (Pearson r) ---\n');
fprintf('    reflection change  vs decorr : %+.2f   (bulk impedance / occlusion)\n', cReflVsDecorr);
fprintf('    standoff  change   vs decorr : %+.2f   (coupling gap)\n', cStandVsDecorr);
fprintf('    sensor-eye DRR chg vs decorr : %+.2f   (anatomy in line of sight)\n', cDRRVsDecorr);
fprintf('  Reading: a large contiguous blob + a high coincidence r on ONE of the\n');
fprintf('  three causes = a real, localized medium/coupling change, not TR noise.\n');
fprintf('  Saved: %s\n%s\n\n', out_mat, bar);


%% ========================================================================
%  LOCAL FUNCTIONS
%% ========================================================================

function sct = load_cbct_geo(fpath)
%LOAD_CBCT_GEO Load a resampled CBCT struct (CBCT1 or CBCT3) with geometry.
    if ~isfile(fpath)
        error('diagnostic_blind_fidelity:CBCTNotFound', 'CBCT not found: %s', fpath);
    end
    d = load(fpath);
    if isfield(d, 'CBCT1_resampled')
        sct = d.CBCT1_resampled;
    elseif isfield(d, 'CBCT3_resampled')
        sct = d.CBCT3_resampled;
    else
        error('diagnostic_blind_fidelity:NoCBCTStruct', ...
            'CBCT1_resampled / CBCT3_resampled not found in %s', fpath);
    end
    if ~isfield(sct, 'origin')  || isempty(sct.origin);  sct.origin  = [0 0 0]; end
    if ~isfield(sct, 'couchMask') || isempty(sct.couchMask)
        sct.couchMask = false(size(sct.bodyMask));
    end
end


function [Z, medium] = build_impedance(sct, config)
%BUILD_IMPEDANCE Acoustic impedance Z = rho*c on the CBCT grid.
%   create_acoustic_medium (the single medium builder) + the coupling-bath
%   override (outside body / couch -> uniform water), matching the forward sim.
%   NOTE: the ~8-line bath is duplicated from run_standalone_field's private
%   build_medium_with_bath (not externally callable); flagged in the summary.
    medium = create_acoustic_medium(sct, config);
    ud = getf(config, 'uniform_density',     1000);
    uc = getf(config, 'uniform_sound_speed', 1540);
    if isfield(sct, 'bodyMask')
        outside = ~logical(sct.bodyMask);
        if isfield(sct, 'couchMask') && ~isempty(sct.couchMask)
            outside = outside | logical(sct.couchMask);
        end
        medium.density(outside)     = ud;
        medium.sound_speed(outside) = uc;
    end
    Z = medium.density .* medium.sound_speed;
end


function iso = resolve_iso_mm(sinfo, doseGrid, origin, spacing)
%RESOLVE_ISO_MM Look-direction target: isocenter if known, else dose max voxel.
    iso = [];
    if isfield(sinfo, 'aim_target_mm') && numel(sinfo.aim_target_mm) == 3 ...
            && all(isfinite(sinfo.aim_target_mm))
        iso = sinfo.aim_target_mm(:)';
        return;
    end
    [~, k] = max(doseGrid(:));
    [r, c, s] = ind2sub(size(doseGrid), k);
    iso = [origin(1) + (c-1)*spacing(1), ...
           origin(2) + (r-1)*spacing(2), ...
           origin(3) + (s-1)*spacing(3)];
end


function src_mm = beam_source_points(doseGrid, origin, spacing, frac, maxpts)
%BEAM_SOURCE_POINTS Thinned set of high-dose (beam-shaped) source voxels, in mm.
    dmax = max(doseGrid(:));
    idx = find(doseGrid >= frac * dmax);
    if isempty(idx)
        [~, idx] = max(doseGrid(:));
    end
    if numel(idx) > maxpts
        idx = idx(round(linspace(1, numel(idx), maxpts)));
    end
    [r, c, s] = ind2sub(size(doseGrid), idx);
    src_mm = [origin(1) + (c-1)*spacing(1), ...
              origin(2) + (r-1)*spacing(2), ...
              origin(3) + (s-1)*spacing(3)];
end


function R = reflection_along_ray(Zvol, origin, spacing, p_src, p_dst, nsamp)
%REFLECTION_ALONG_RAY Cumulative reflection proxy sum_i ((dZ)/(sumZ))^2 along
%   a straight ray. Out-of-bounds samples fall back to water impedance so the
%   coupling bath does not manufacture spurious interfaces.
    Z0  = 1000 * 1540;
    t   = linspace(0, 1, nsamp)';
    pts = p_src + t .* (p_dst - p_src);          % nsamp x 3 (implicit expansion)
    Zs  = sample_pts(Zvol, origin, spacing, pts, Z0);
    dZ  = diff(Zs);
    sZ  = Zs(1:end-1) + Zs(2:end);
    sZ(sZ == 0) = eps;
    r   = dZ ./ sZ;
    R   = sum(r.^2);
end


function d = standoff_to_body(bodyMask, origin, spacing, p0mm, dir_unit, step, maxd)
%STANDOFF_TO_BODY March from an element toward iso; return mm to first body
%   voxel, or NaN if the body is never reached within maxd (air gap / no contact).
    k   = (0:floor(maxd / step))';
    pts = p0mm + (k * step) .* dir_unit;         % (n+1) x 3
    hit = sample_pts(double(bodyMask), origin, spacing, pts, 0) > 0.5;
    f   = find(hit, 1, 'first');
    if isempty(f); d = NaN; else; d = (f - 1) * step; end
end


function val = drr_along_ray(HUvol, origin, spacing, p_src, p_dst, nsamp)
%DRR_ALONG_RAY Integrated radiographic thickness (HU above air) along a ray:
%   a sensor-eye projection value for one element's line of sight.
    t   = linspace(0, 1, nsamp)';
    pts = p_src + t .* (p_dst - p_src);
    L   = norm(p_dst - p_src);
    ds  = L / max(nsamp - 1, 1);
    hu  = sample_pts(HUvol, origin, spacing, pts, -1000);
    val = sum(max(hu + 1000, 0)) * ds;           % 0 in air, grows through tissue/bone
end


function vals = sample_pts(vol, origin, spacing, pts, oob_val)
%SAMPLE_PTS Nearest-voxel values of vol at physical points pts (M x 3, x,y,z mm).
%   vol is (row,col,slice) = (Y,X,Z); returns oob_val for points outside the grid.
    c = round((pts(:,1) - origin(1)) / spacing(1)) + 1;   % X -> col
    r = round((pts(:,2) - origin(2)) / spacing(2)) + 1;   % Y -> row
    s = round((pts(:,3) - origin(3)) / spacing(3)) + 1;   % Z -> slice
    sz = size(vol);
    inb = r >= 1 & r <= sz(1) & c >= 1 & c <= sz(2) & s >= 1 & s <= sz(3);
    vals = repmat(double(oob_val), size(pts, 1), 1);
    if any(inb)
        lin = sub2ind(sz, r(inb), c(inb), s(inb));
        vals(inb) = double(vol(lin));
    end
end


function r = pearson(a, b)
%PEARSON Zero-lag Pearson correlation of two vectors (0 if either is flat).
    a = a(:) - mean(a(:));
    b = b(:) - mean(b(:));
    na = norm(a); nb = norm(b);
    if na == 0 || nb == 0
        r = 0;
    else
        r = (a' * b) / (na * nb);
    end
end


function r = pearson_masked(x, y, mask)
%PEARSON_MASKED Pearson correlation of x,y over finite, masked entries.
    x = x(mask); y = y(mask);
    ok = isfinite(x) & isfinite(y);
    if nnz(ok) < 3
        r = NaN; return;
    end
    r = pearson(x(ok), y(ok));
end


function n = largest_blob(bw)
%LARGEST_BLOB Size (element count) of the largest 8-connected true region.
    n = 0;
    if ~any(bw(:)); return; end
    try
        cc = bwconncomp(bw, 8);
        if cc.NumObjects > 0
            n = max(cellfun(@numel, cc.PixelIdxList));
        end
    catch
        n = nnz(bw);   % Image Processing Toolbox absent: fall back to total count
    end
end


function show_map(img, ttl, clim)
%SHOW_MAP imagesc a 32x32 element map with a colorbar and equal aspect.
    imagesc(img);
    if nargin >= 3 && ~isempty(clim); caxis(clim); end
    axis image; colorbar; title(ttl, 'Interpreter', 'tex');
    xlabel('element col (ex)'); ylabel('element row (ez)');
end


function cmap = diverging_map()
%DIVERGING_MAP Simple blue-white-red map for signed difference images.
    n = 128;
    t = linspace(0, 1, n)';
    top = [t, t, ones(n,1)];                 % blue -> white
    bot = [ones(n,1), flipud(t), flipud(t)]; % white -> red
    cmap = [top; bot];
end


function v = getf(s, f, d)
%GETF Struct field with default.
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
