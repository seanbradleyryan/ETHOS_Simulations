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
%       identical p0 through the CT_1 medium and the CT_3 medium. Per element:
%         - Zero-lag Pearson correlation map.
%         - LAG-TOLERANT peak cross-correlation map + a LAG map (us). Together
%           these separate a bulk time-of-flight shift (sound-speed change: low
%           zero-lag corr but high peak corr, coherent lag) from TRUE signal loss
%           (low even after lag alignment => "truly lost" channels). A contiguous
%           truly-lost blob = physical occlusion; scatter = numerics.
%         - Energy-ratio map: amplitude loss even when the waveform shape holds.
%
%   [2] IMPEDANCE / REFLECTION + TISSUE RAY INTEGRAL (geometry only, no k-Wave).
%       For each element, cast rays from the high-dose (beam) source set to the
%       element through the CT_1 and CT_3 media and accumulate a reflection
%       proxy R = sum_i ((Z_{i+1}-Z_i)/(Z_{i+1}+Z_i))^2 (Z = rho*c), plus the
%       GAS and BONE fractions along the ray. Maps are CT_3 - CT_1 differences.
%       Reflection-up with gas-up but bone-flat = a gas pocket (near-total
%       acoustic wall, invisible on a radiograph); bone-up raises both. The
%       reflection scatter vs (1 - peak corr) SELF-CALIBRATES the "how much
%       change kills a channel" threshold against [1].
%
%   [4] STANDOFF / COUPLING + SENSOR-EYE PROJECTION (geometry only, no k-Wave).
%       Per element, ABSOLUTE distance along its inward look-direction (toward
%       iso) to the first body voxel, for CT_1 and CT_3 (32x32 maps + diff). An
%       air gap (NaN standoff) is an unambiguous, catastrophic coupling failure.
%       A per-element sensor-eye DRR (integrated HU along the ray) shows total
%       structure in an element's line of sight. Coincidence of these with the
%       [1] loss localizes the cause (coupling gap vs anatomy vs bulk impedance).
%
%   [F] SENSOR PLACEMENT CT_1 vs CT_3 (geometry only, no k-Wave).
%       The maps above froze the sensor to isolate the medium. [F] instead
%       RE-PLACES the array on each geometry (determine_sensor_mask, same dose,
%       only bodyMask differs) and maps the per-element displacement plus the
%       center shift, aim-normal angle, and tilt change -- the direct test of
%       "does the sensor tilt/shift between CT_1 and CT_3".
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
%   - ONE figure window with four tabs ([1] per-channel, [2] reflection & tissue,
%     [4] coupling & sensor-eye, [F] sensor placement) + a saved results .mat
%     (all maps + metrics). Array maps are labeled by the physical direction each
%     element axis points (from element_positions_mm), so image up/down maps to
%     real anatomy despite the array tilt.
%   - Console summary: truly-lost channel count, median bulk lag, air-gap count,
%     max standoff change, sensor re-placement shift, and the coincidence
%     correlations that say whether the loss is gas-, bone-, coupling-, or
%     bulk-impedance-driven.
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
DIAG.corr_dead_thresh = 0.50;  % [1] correlation below this => channel counted "dead"
DIAG.energy_lo        = 0.50;  % [1] energy ratio (CT3/CT1) outside [lo, 1/lo] => flagged
DIAG.energy_floor_frac = 0.05; % [1] skip elements below this fraction of the peak channel energy
DIAG.maxlag_samples   = 600;   % [1] lag search half-width (samples). Wide: an 18 mm
                               %     water-path change alone is ~12 us of lag. FFT
                               %     cross-correlation makes a big window cheap.
DIAG.air_hu           = -300;  % [3] HU below this = gas/air along a ray
DIAG.bone_hu          = 300;   % [3] HU above this = bone along a ray
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
if max(abs(sct1.origin(:) - sct3.origin(:))) > 1e-3
    warning('diagnostic_blind_fidelity:OriginMismatch', ...
        ['CT_1 and CT_3 origins differ; [2]/[4] sample CT_3 with the CT_1 origin, ' ...
         'so the fixed rays may be misregistered on CT_3.']);
end

[Z1, med1]  = build_impedance(sct1, CONFIG);   % Z = rho*c on the CT_1 grid
[Z3, med3]  = build_impedance(sct3, CONFIG);   % Z = rho*c on the CT_3 grid
body1 = logical(sct1.bodyMask);
body3 = logical(sct3.bodyMask);
HU1   = sct1.cubeHU;
HU3   = sct3.cubeHU;

doseGrid = out.doseGrid;                    % CT_1-body-masked dose, native grid
iso_mm   = resolve_iso_mm(sinfo, doseGrid, origin, spacing);
fprintf('      Isocenter (look target): [%.1f %.1f %.1f] mm\n', iso_mm);

% Bulk medium check: is the whole-body average sound speed shifted between CT_1
% and CT_3? A coherent shift (not local) points at CBCT HU / segmentation drift
% rather than moved anatomy, and it is what produces a uniform time-of-flight
% lag and global TR defocus. Compared inside each body and their intersection.
bulk = struct();
bulk.mean_c_CT1  = mean(med1.sound_speed(body1));
bulk.mean_c_CT3  = mean(med3.sound_speed(body3));
both = body1 & body3;
bulk.mean_c_CT1_shared = mean(med1.sound_speed(both));
bulk.mean_c_CT3_shared = mean(med3.sound_speed(both));
bulk.dc_shared   = bulk.mean_c_CT3_shared - bulk.mean_c_CT1_shared;
bulk.mean_rho_CT1_shared = mean(med1.density(both));
bulk.mean_rho_CT3_shared = mean(med3.density(both));
fprintf('      Body mean sound speed: CT1 %.1f, CT3 %.1f m/s (shared voxels: %+.1f m/s, %.2f%%)\n', ...
    bulk.mean_c_CT1, bulk.mean_c_CT3, bulk.dc_shared, ...
    100 * bulk.dc_shared / max(bulk.mean_c_CT1_shared, eps));

%% =================== DIAGNOSTIC [1]: PER-CHANNEL ====================== %%
%  Average sensor-voxel traces per element (authoritative voxel->element map),
%  then per element compute Pearson correlation and energy ratio CT_3 vs CT_1.

fprintf('[3/5] Diagnostic [1]: per-channel signal change...\n');
if numel(velem) ~= size(sino1, 1)
    warning('diagnostic_blind_fidelity:VoxelCountMismatch', ...
        'voxel_element_idx (%d) ~= sinogram rows (%d); channel map may be misaligned.', ...
        numel(velem), size(sino1, 1));
end
corrMap = nan(Ntot, 1);      % zero-lag Pearson(CT1 trace, CT3 trace)
peakMap = nan(Ntot, 1);      % peak normalized cross-correlation over lags
lagMap  = nan(Ntot, 1);      % lag (samples) at the peak; sign = CT3 delayed vs CT1
enerMap = nan(Ntot, 1);      % ||CT3|| / ||CT1||
E1      = nan(Ntot, 1);      % CT1 trace energy (for the low-signal floor)
for e = 1:Ntot
    rows = find(velem == e);
    if isempty(rows); continue; end
    t1 = mean(sino1(rows, :), 1);
    t3 = mean(sino3(rows, :), 1);
    E1(e) = norm(t1);
    corrMap(e) = pearson(t1, t3);
    [peakMap(e), lagMap(e)] = peak_xcorr(t1, t3, DIAG.maxlag_samples);
    if E1(e) > 0
        enerMap(e) = norm(t3) / E1(e);
    end
end
% Low-signal elements can't be meaningfully "blocked"; excluding them keeps the
% dead-channel map (and the coincidence stats) from being polluted by noise.
if any(~isnan(E1))
    sig_floor = DIAG.energy_floor_frac * max(E1(~isnan(E1)));
    lowsig = E1 < sig_floor;
    corrMap(lowsig) = NaN;  peakMap(lowsig) = NaN;
    lagMap(lowsig)  = NaN;  enerMap(lowsig) = NaN;
end
% Lag in microseconds (bulk sound-speed / time-of-flight change).
dt_us  = getf(sr, 'sinogram_dt', NaN) * 1e6;
lagUs  = lagMap * dt_us;
% Two "decorrelation" measures: zero-lag (swamped by any bulk time shift) and
% lag-tolerant (true signal loss after removing a uniform shift).
decorr      = 1 - corrMap;      % zero-lag
decorr_peak = 1 - peakMap;      % lag-tolerant: this is the physical "signal lost"
% Zero-lag dead (may just be shifted) vs truly-lost (low even after alignment).
deadMask = (corrMap < DIAG.corr_dead_thresh) | ...
           (enerMap < DIAG.energy_lo) | (enerMap > 1/DIAG.energy_lo);
deadMask(isnan(corrMap)) = false;
trueDead = peakMap < DIAG.corr_dead_thresh;               % lost even allowing a lag
trueDead(isnan(peakMap)) = false;

corrImg  = reshape(corrMap,  [N, N]);   % rows = ez, cols = ex
peakImg  = reshape(peakMap,  [N, N]);
lagImg   = reshape(lagUs,    [N, N]);
enerImg  = reshape(enerMap,  [N, N]);
deadImg  = reshape(double(deadMask), [N, N]);
trueDeadImg = reshape(double(trueDead), [N, N]);

%% =================== DIAGNOSTIC [2]: REFLECTION RAY ==================== %%
%  Cumulative reflection proxy along source->element rays, CT_1 and CT_3.

fprintf('[4/5] Diagnostic [2]: impedance / reflection ray integral...\n');
src_mm = beam_source_points(doseGrid, origin, spacing, ...
    DIAG.src_dose_frac, DIAG.max_src_points);
fprintf('      Source points (beam >= %.0f%% max dose): %d\n', ...
    100*DIAG.src_dose_frac, size(src_mm,1));

refl1 = nan(Ntot, 1);  refl3 = nan(Ntot, 1);
air1  = nan(Ntot, 1);  air3  = nan(Ntot, 1);   % [3] gas fraction along the ray
bone1 = nan(Ntot, 1);  bone3 = nan(Ntot, 1);   % [3] bone fraction along the ray
nsrc  = size(src_mm, 1);
for e = 1:Ntot
    pe = epos(e, :);
    acc1 = [0 0 0]; acc3 = [0 0 0];   % [reflection, airFrac, boneFrac]
    for s = 1:nsrc
        acc1 = acc1 + ray_metrics(Z1, HU1, origin, spacing, src_mm(s,:), pe, ...
            DIAG.ray_nsamp, DIAG.air_hu, DIAG.bone_hu);
        acc3 = acc3 + ray_metrics(Z3, HU3, origin, spacing, src_mm(s,:), pe, ...
            DIAG.ray_nsamp, DIAG.air_hu, DIAG.bone_hu);
    end
    refl1(e) = acc1(1)/nsrc;  air1(e) = acc1(2)/nsrc;  bone1(e) = acc1(3)/nsrc;
    refl3(e) = acc3(1)/nsrc;  air3(e) = acc3(2)/nsrc;  bone3(e) = acc3(3)/nsrc;
end
dRefl = refl3 - refl1;
dAir  = air3  - air1;      % gas that entered (or left) each element's paths
dBone = bone3 - bone1;     % bone that entered (or left) each element's paths
refl1Img = reshape(refl1, [N, N]);
refl3Img = reshape(refl3, [N, N]);
dReflImg = reshape(dRefl, [N, N]);
dAirImg  = reshape(dAir,  [N, N]);
dBoneImg = reshape(dBone, [N, N]);

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

%% =================== [F] SENSOR PLACEMENT: CT_1 vs CT_3 =============== %%
%  The maps above froze the sensor (CT_1 placement) to isolate the medium. This
%  block instead RE-PLACES the array on each geometry (same dose, only bodyMask
%  differs) to test whether the sensor itself tilts/shifts between CT_1 and CT_3.
%  Placement is computed natively here, so it may differ slightly from the
%  engine's internal feed; the CT_1 -> CT_3 delta is what matters.

fprintf('[F]   Sensor placement comparison (CT_1 vs CT_3)...\n');
placeShiftImg = []; placeStats = struct('ok', false);
try
    beam_md = load_beam_metadata_local(processed_dir);
    fd_place = struct('dose_Gy', doseGrid, 'gantry_angle', 0, ...
        'origin', origin, 'spacing', spacing, 'dimensions', size(doseGrid));
    [~, si1] = determine_sensor_mask(sct1, fd_place, beam_md, CONFIG);
    [~, si3] = determine_sensor_mask(sct3, fd_place, beam_md, CONFIG);
    ep1 = si1.element_positions_mm;
    ep3 = si3.element_positions_mm;
    placeStats.tilt1_deg   = getf(si1, 'tilt_angle_deg', NaN);
    placeStats.tilt3_deg   = getf(si3, 'tilt_angle_deg', NaN);
    placeStats.center1_mm  = getf(si1, 'sensor_center_mm', [NaN NaN NaN]);
    placeStats.center3_mm  = getf(si3, 'sensor_center_mm', [NaN NaN NaN]);
    placeStats.center_shift_mm = norm(placeStats.center3_mm - placeStats.center1_mm);
    placeStats.normal_angle_deg = vec_angle_deg( ...
        getf(si1, 'aim_normal', [0 -1 0]), getf(si3, 'aim_normal', [0 -1 0]));
    if isequal(size(ep1), size(ep3)) && size(ep1,1) == Ntot
        shift = sqrt(sum((ep3 - ep1).^2, 2));    % per-element displacement (mm)
        placeShiftImg = reshape(shift, [N, N]);
        placeStats.median_elem_shift_mm = median(shift);
        placeStats.max_elem_shift_mm    = max(shift);
    else
        warning('diagnostic_blind_fidelity:PlacementSizeMismatch', ...
            'CT_1/CT_3 element counts differ (%d vs %d); per-element shift map skipped.', ...
            size(ep1,1), size(ep3,1));
    end
    placeStats.ok = true;
catch ME
    warning('diagnostic_blind_fidelity:PlacementCompareFail', ...
        'Sensor placement comparison failed: %s', ME.message);
end

%% =================== COINCIDENCE / SELF-CALIBRATION =================== %%
%  Do the geometry changes explain the SIGNAL LOSS? Correlate each candidate
%  cause with the lag-tolerant decorrelation (1 - peak xcorr) from [1] -- i.e.
%  true signal loss after removing any bulk time shift -- over elements with a
%  valid signal. Zero-lag versions kept alongside for comparison.

valid = ~isnan(decorr_peak);
cReflVsLoss  = pearson_masked(dRefl,  decorr_peak, valid);
cAirVsLoss   = pearson_masked(dAir,   decorr_peak, valid);
cBoneVsLoss  = pearson_masked(dBone,  decorr_peak, valid);
cStandVsLoss = pearson_masked(dStand, decorr_peak, valid);
cDRRVsLoss   = pearson_masked(dDRR,   decorr_peak, valid);
% Zero-lag reference (comparable to the earlier run).
cReflVsDecorr  = pearson_masked(dRefl,  decorr, ~isnan(decorr));
cStandVsDecorr = pearson_masked(dStand, decorr, ~isnan(decorr));
cDRRVsDecorr   = pearson_masked(dDRR,   decorr, ~isnan(decorr));

n_dead     = nnz(deadMask);
n_truedead = nnz(trueDead);
n_air      = nnz(airGap);
blob       = largest_blob(trueDeadImg > 0.5);      % blob of TRULY-lost channels
maxStand   = max(abs(dStand(valid)));
med_lag_us = median(lagUs(valid), 'omitnan');
% Lag-window saturation: if many elements pin near +/- maxlag, the window is too
% small and "truly lost" is an artifact (widen DIAG.maxlag_samples and re-run).
rail_frac  = mean(abs(lagMap(valid)) >= 0.9 * DIAG.maxlag_samples);
lag_lo_us  = min(lagUs(valid));
lag_hi_us  = max(lagUs(valid));

%% =========================== PLOTS (tabbed) =========================== %%
%  One window; each diagnostic on its own tab. Array maps are labeled by the
%  physical direction each element axis (ex, ez) points, derived from
%  element_positions_mm, so "up/down in the image" maps to real anatomy.

[xlab, ylab] = array_axis_labels(epos, N);
hFig = figure('Name', 'Blind Fidelity Diagnostic (CT_1 vs CT_3)', 'Color', 'w');
tg = uitabgroup(hFig);

% --- Tab [1]: per-channel ---
tabA = uitab(tg, 'Title', '[1] Per-channel');
tl = tiledlayout(tabA, 2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile(tl); show_map(corrImg, 'Zero-lag corr', [0 1], xlab, ylab);
nexttile(tl); show_map(peakImg, 'Peak x-corr (lag-tolerant)', [0 1], xlab, ylab);
nexttile(tl); show_map(lagImg,  'Lag at peak (\mus)', [], xlab, ylab); colormap(gca, diverging_map());
nexttile(tl); show_map(enerImg, 'Energy ratio CT3/CT1', [0 2], xlab, ylab);
nexttile(tl); show_map(deadImg, sprintf('Zero-lag dead (%d)', n_dead), [0 1], xlab, ylab); colormap(gca, 'hot');
nexttile(tl); show_map(trueDeadImg, sprintf('Truly lost (%d)', n_truedead), [0 1], xlab, ylab); colormap(gca, 'hot');
title(tl, 'Diagnostic [1]: low zero-lag but high peak x-corr = bulk time shift, not loss');

% --- Tab [2]: reflection + tissue ---
tabB = uitab(tg, 'Title', '[2] Reflection & tissue');
tl = tiledlayout(tabB, 2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile(tl); show_map(refl1Img, 'Reflection proxy CT1', [], xlab, ylab);
nexttile(tl); show_map(refl3Img, 'Reflection proxy CT3', [], xlab, ylab);
nexttile(tl); show_map(dReflImg, '\Delta Reflection (CT3-CT1)', [], xlab, ylab); colormap(gca, diverging_map());
nexttile(tl); show_map(dAirImg,  '\Delta gas fraction in path', [], xlab, ylab); colormap(gca, diverging_map());
nexttile(tl); show_map(dBoneImg, '\Delta bone fraction in path', [], xlab, ylab); colormap(gca, diverging_map());
nexttile(tl);
scatter(dRefl(valid), decorr_peak(valid), 14, 'filled'); grid on;
xlabel('\Delta reflection proxy (CT3 - CT1)'); ylabel('1 - peak x-corr');
title(sprintf('Self-calibration: r = %.2f', cReflVsLoss));
title(tl, 'Diagnostic [2]: \Deltareflection up with \Deltagas up & \Deltabone flat = gas moved in');

% --- Tab [4]: coupling + sensor-eye ---
tabC = uitab(tg, 'Title', '[4] Coupling & sensor-eye');
tl = tiledlayout(tabC, 2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile(tl); show_map(stand1Img, 'Standoff CT1 (mm)', [], xlab, ylab);
nexttile(tl); show_map(stand3Img, 'Standoff CT3 (mm)', [], xlab, ylab);
nexttile(tl); show_map(dStandImg, '\Delta Standoff (mm)', [], xlab, ylab); colormap(gca, diverging_map());
nexttile(tl); show_map(drr1Img, 'Sensor-eye DRR CT1', [], xlab, ylab);
nexttile(tl); show_map(drr3Img, 'Sensor-eye DRR CT3', [], xlab, ylab);
nexttile(tl); show_map(dDRRImg, '\Delta DRR (CT3-CT1)', [], xlab, ylab); colormap(gca, diverging_map());
title(tl, sprintf('Diagnostic [4]: air-gap (no-contact) channels = %d', n_air));

% --- Tab [F]: sensor placement CT_1 vs CT_3 ---
tabD = uitab(tg, 'Title', '[F] Sensor placement');
if placeStats.ok && ~isempty(placeShiftImg)
    tl = tiledlayout(tabD, 1, 1, 'Padding', 'compact');
    nexttile(tl); show_map(placeShiftImg, 'Per-element placement shift CT1\rightarrowCT3 (mm)', ...
        [], xlab, ylab); colormap(gca, 'parula');
    title(tl, sprintf(['center shift %.1f mm | aim-normal angle %.1f deg | ' ...
        'tilt %.0f\\rightarrow%.0f deg | median/max elem shift %.1f/%.1f mm'], ...
        placeStats.center_shift_mm, placeStats.normal_angle_deg, ...
        placeStats.tilt1_deg, placeStats.tilt3_deg, ...
        getf(placeStats,'median_elem_shift_mm',NaN), getf(placeStats,'max_elem_shift_mm',NaN)));
else
    ax = axes('Parent', tabD); axis(ax, 'off');
    text(ax, 0.5, 0.5, 'Sensor placement comparison unavailable (see warning).', ...
        'HorizontalAlignment', 'center');
end

%% =========================== SAVE + SUMMARY ========================== %%

if ~isfolder(DIAG.save_dir); mkdir(DIAG.save_dir); end
[~, dose_tag] = fileparts(CONFIG.dose_filename);
results = struct();
results.config       = CONFIG;
results.diag_params  = DIAG;
results.elements_per_side = N;
results.maps = struct('corr', corrImg, 'peak_xcorr', peakImg, 'lag_us', lagImg, ...
    'energy', enerImg, 'dead', deadImg, 'true_dead', trueDeadImg, ...
    'refl_CT1', refl1Img, 'refl_CT3', refl3Img, 'd_refl', dReflImg, ...
    'd_air_frac', dAirImg, 'd_bone_frac', dBoneImg, ...
    'standoff_CT1', stand1Img, 'standoff_CT3', stand3Img, 'd_standoff', dStandImg, ...
    'drr_CT1', drr1Img, 'drr_CT3', drr3Img, 'd_drr', dDRRImg, ...
    'placement_shift', placeShiftImg);
results.coincidence = struct( ...
    'refl_vs_loss', cReflVsLoss, 'air_vs_loss', cAirVsLoss, 'bone_vs_loss', cBoneVsLoss, ...
    'standoff_vs_loss', cStandVsLoss, 'drr_vs_loss', cDRRVsLoss, ...
    'refl_vs_zerolag', cReflVsDecorr, 'standoff_vs_zerolag', cStandVsDecorr, ...
    'drr_vs_zerolag', cDRRVsDecorr);
results.placement = placeStats;
results.medium_bulk = bulk;
results.summary = struct('n_dead_zerolag', n_dead, 'n_truly_lost', n_truedead, ...
    'n_air_gap', n_air, 'largest_lost_blob', blob, ...
    'max_abs_d_standoff_mm', maxStand, 'median_lag_us', med_lag_us, ...
    'lag_rail_fraction', rail_frac, 'lag_range_us', [lag_lo_us, lag_hi_us]);
out_mat = fullfile(DIAG.save_dir, sprintf('blind_fidelity_%s.mat', dose_tag));
save(out_mat, 'results', '-v7.3');

fprintf('\n%s\n  SUMMARY\n%s\n', bar, bar);
fprintf('  Channels dead at zero lag (corr<%.2f/energy)          : %d / %d\n', ...
    DIAG.corr_dead_thresh, n_dead, nnz(valid));
fprintf('  Channels TRULY lost (peak x-corr<%.2f, lag-tolerant)  : %d / %d\n', ...
    DIAG.corr_dead_thresh, n_truedead, nnz(valid));
fprintf('  Median lag (bulk time-of-flight shift, us)            : %+.2f  [range %+.2f..%+.2f]\n', ...
    med_lag_us, lag_lo_us, lag_hi_us);
fprintf('  Lag-window saturation (|lag| >= 0.9*maxlag)           : %.0f%%', 100*rail_frac);
if rail_frac > 0.1
    fprintf('   <-- WIDEN DIAG.maxlag_samples (%d) AND RE-RUN\n', DIAG.maxlag_samples);
else
    fprintf('   (ok)\n');
end
fprintf('  Body mean sound speed CT1->CT3 (shared voxels)        : %.1f -> %.1f m/s (%+.2f%%)\n', ...
    bulk.mean_c_CT1_shared, bulk.mean_c_CT3_shared, ...
    100 * bulk.dc_shared / max(bulk.mean_c_CT1_shared, eps));
fprintf('  Largest contiguous TRULY-lost blob (elements)         : %d\n', blob);
fprintf('  Air-gap / no-contact elements (NaN standoff)          : %d\n', n_air);
fprintf('  Max |standoff change| CT3 vs CT1 (mm)                 : %.1f\n', maxStand);
if placeStats.ok
    fprintf('  Sensor re-placement CT1->CT3: center %.1f mm, normal %.1f deg, ', ...
        placeStats.center_shift_mm, placeStats.normal_angle_deg);
    fprintf('max elem shift %.1f mm\n', getf(placeStats, 'max_elem_shift_mm', NaN));
end
fprintf('  --- coincidence with TRUE signal loss (1 - peak x-corr), Pearson r ---\n');
fprintf('    reflection change  : %+.2f   (bulk impedance / occlusion)\n', cReflVsLoss);
fprintf('    gas fraction chg   : %+.2f   (bowel/stomach gas moved into path)\n', cAirVsLoss);
fprintf('    bone fraction chg  : %+.2f   (bone moved into path)\n', cBoneVsLoss);
fprintf('    standoff change    : %+.2f   (coupling gap)\n', cStandVsLoss);
fprintf('    sensor-eye DRR chg : %+.2f   (total anatomy in line of sight)\n', cDRRVsLoss);
fprintf('  (zero-lag reference: refl %+.2f, standoff %+.2f, drr %+.2f)\n', ...
    cReflVsDecorr, cStandVsDecorr, cDRRVsDecorr);
fprintf('  Reading: whichever cause has the highest r vs TRUE loss is the mechanism;\n');
fprintf('  gas high + bone flat = a gas pocket, not bone. A large median lag with\n');
fprintf('  high peak x-corr means the medium mainly shifted timing (bulk c change).\n');
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


function m = ray_metrics(Zvol, HUvol, origin, spacing, p_src, p_dst, nsamp, air_hu, bone_hu)
%RAY_METRICS Along one straight source->element ray, return
%   m = [reflection_proxy, gas_fraction, bone_fraction]:
%     - reflection proxy = sum_i ((dZ)/(sumZ))^2 (Z = rho*c); out-of-bounds
%       samples fall back to water so the coupling bath adds no false interface.
%     - gas/bone fraction = fraction of IN-BODY-range samples below air_hu /
%       above bone_hu (HU, so gas and bone are told apart: gas spikes reflection
%       but not radiographic thickness, bone raises both).
    Z0  = 1000 * 1540;
    t   = linspace(0, 1, nsamp)';
    pts = p_src + t .* (p_dst - p_src);          % nsamp x 3 (implicit expansion)
    Zs  = sample_pts(Zvol,  origin, spacing, pts, Z0);
    HUs = sample_pts(HUvol, origin, spacing, pts, -1000);
    dZ  = diff(Zs);
    sZ  = Zs(1:end-1) + Zs(2:end);
    sZ(sZ == 0) = eps;
    R   = sum((dZ ./ sZ).^2);
    m   = [R, mean(HUs < air_hu), mean(HUs > bone_hu)];
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


function [rpk, lpk] = peak_xcorr(a, b, maxlag)
%PEAK_XCORR Peak normalized cross-correlation of a,b over lags [-maxlag,maxlag],
%   via FFT so a wide window is cheap. Returns the peak correlation and the lag
%   (samples) at which it occurs (positive => b delayed relative to a). Globally
%   normalized (norm(a)*norm(b)), so a pure time shift scores ~1 even when
%   zero-lag Pearson is low. No Signal Processing Toolbox dependency.
    a = a(:) - mean(a(:));
    b = b(:) - mean(b(:));
    na = norm(a); nb = norm(b);
    if na == 0 || nb == 0
        rpk = 0; lpk = 0; return;
    end
    n  = numel(a);
    nf = 2^nextpow2(2*n - 1);
    cc = real(ifft(fft(a, nf) .* conj(fft(b, nf))));
    % Reorder to lags -(n-1) .. (n-1); cc(t) = sum a(k) b(k-lag).
    cc   = [cc(nf-n+2:nf); cc(1:n)];
    lags = (-(n-1):(n-1))';
    if ~isempty(maxlag) && maxlag > 0
        keep = abs(lags) <= maxlag;
        cc = cc(keep); lags = lags(keep);
    end
    [mx, ix] = max(cc);
    rpk = mx / (na * nb);
    lpk = lags(ix);
end


function [xlab, ylab] = array_axis_labels(epos, N)
%ARRAY_AXIS_LABELS Physical direction each element axis (ex cols, ez rows)
%   points, from element_positions_mm. e = (ex-1)*N + ez, so reshape([N,N,3])
%   gives E(ez, ex, :). Lets the reader map "up/down in the image" to anatomy.
    E   = reshape(epos, [N, N, 3]);              % E(ez, ex, :)
    gex = squeeze(mean(mean(diff(E, 1, 2), 1), 2));   % mm per +1 ex (columns)
    gez = squeeze(mean(mean(diff(E, 1, 1), 1), 2));   % mm per +1 ez (rows)
    xlab = sprintf('ex \\rightarrow %s', anatomy_of(gex));
    ylab = sprintf('ez \\rightarrow %s', anatomy_of(gez));
end


function s = anatomy_of(v)
%ANATOMY_OF Dominant HFS/DICOM anatomical direction of a physical vector
%   v = [x y z] mm: +x LEFT, +y POSTERIOR, +z SUPERIOR.
    names = {'RIGHT', 'LEFT'; 'ANTERIOR', 'POSTERIOR'; 'INFERIOR', 'SUPERIOR'};
    [~, k] = max(abs(v));
    if v(k) >= 0; s = names{k, 2}; else; s = names{k, 1}; end
end


function ang = vec_angle_deg(u, w)
%VEC_ANGLE_DEG Angle (deg) between two vectors; NaN-safe.
    u = u(:); w = w(:);
    du = norm(u) * norm(w);
    if du == 0; ang = NaN; return; end
    ang = acosd(max(-1, min(1, (u' * w) / du)));
end


function bm = load_beam_metadata_local(processed_dir)
%LOAD_BEAM_METADATA_LOCAL Beam metadata for determine_sensor_mask exclusion
%   zones. Reads processed/metadata.mat (metadata.beam_metadata); [] if absent.
    bm = [];
    f = fullfile(processed_dir, 'metadata.mat');
    if isfile(f)
        md = load(f, 'metadata');
        if isfield(md, 'metadata') && isfield(md.metadata, 'beam_metadata')
            bm = md.metadata.beam_metadata;
        end
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


function show_map(img, ttl, clim, xlab, ylab)
%SHOW_MAP imagesc a 32x32 element map with a colorbar and equal aspect.
%   NaN elements (aliased / no signal) render as the axes background.
    if nargin < 4 || isempty(xlab); xlab = 'element col (ex)'; end
    if nargin < 5 || isempty(ylab); ylab = 'element row (ez)'; end
    h = imagesc(img);
    set(h, 'AlphaData', ~isnan(img));            % aliased/no-signal -> transparent
    if nargin >= 3 && ~isempty(clim); caxis(clim); end
    axis image; colorbar; title(ttl, 'Interpreter', 'tex');
    xlabel(xlab, 'Interpreter', 'tex'); ylabel(ylab, 'Interpreter', 'tex');
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
