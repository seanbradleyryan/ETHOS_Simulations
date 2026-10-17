%% =========================================================================
%  DEMO_MBR.m
%  Model-based reconstruction (MBR) of a radiation dose from IRAI sensor
%  signals: a step-by-step, heavily commented walk-through.
%  =========================================================================
%
%  PURPOSE
%    Teach how model-based reconstruction works on the real total plan dose
%    (every beam segment summed, one total per CT) of patient 1194203 /
%    Session_1, and check (gamma analysis) that it lands closer to the true
%    dose than the time-reversal (TR) reconstruction the pipeline uses today.
%    Run it section by section (Ctrl+Enter) and read the comments as you go.
%
%  THE IDEA IN ONE PARAGRAPH
%    TR plays the recorded signals backwards into a model of the patient and
%    looks at where the sound refocuses. It is a one-shot recipe: it never
%    checks whether its image actually explains the measurements. With a
%    sensor on one side of the body only, a band-limited transducer and noise,
%    the refocused image is blurred, has negative lobes and the wrong size.
%    MBR turns the problem around: guess an image, SIMULATE the signals that
%    guess would produce, compare them with the measured signals, and correct
%    the guess where they disagree. Repeat. The answer is the image whose
%    simulated signals best match the measured ones, with a mild preference
%    for images without speckle. In equations, minimize
%        0.5 * || A p - y ||^2  +  lambda * TV(p)      with  p >= 0
%    p = initial pressure image, y = measured signals, A = the "predictor"
%    (k-Wave + linac pulse + transducer), TV = total variation (the prior).
%    The minimizer used here is FISTA (gradient steps with momentum).
%
%  TWO ARMS (both reconstructed on CT1, the only anatomy a live system has)
%    R1: CT1 total dose, signals generated on CT1 anatomy. The model used
%        to reconstruct is exactly the model that made the data (an "inverse
%        crime"), so this is the best case.
%    R3: CT3 total dose, signals generated on CT3 anatomy, reconstructed
%        on CT1. The model is wrong by the anatomy change: the realistic case.
%
%  STEPS (numbering follows the method write-up)
%    1      raw sensor signals (forward k-Wave + pulse + transducer + noise)
%    base   time reversal, as run_single_field_simulation does it
%    2      predictor          A  : image   -> signals
%    3      reverse predictor  A' : signals -> image ("blame map"), either a
%           backwards k-Wave run ('tr') or delay-and-sum ('das')
%    4      check A and A' are exact opposites (adjoint test)
%    5      safe step size (power iteration)
%    6      TV prior and its weight
%    7-9    FISTA loop from the TR image, fixed iteration count
%    10     pressure -> dose
%    11     mismatch vs iteration; seen / unseen split of R1 - R3
%    result gamma analysis and images: truth | TR | MBR
%
%  RUNTIME
%    Every predictor or reverse-predictor call is one full k-Wave run. With
%    the defaults this demo does ~90 runs (setting mbr_tune_tv adds ~190).
%    A GPU is needed in practice. With mbr_adjoint = 'das' about half of
%    those become CPU delay-and-sum passes instead; das_reconstruct loops
%    over every element, so this is not necessarily faster.
%
%  DEPENDENCIES
%    k-Wave (kspaceFirstOrder3D, kWaveGrid, smooth, gaussianFilter), Image
%    Processing Toolbox (padarray), CalcGamma, and in utils/:
%    get_default_config, create_acoustic_medium,
%    determine_sensor_mask_lateral, find_optimal_kwave_size,
%    build_pulse_kernel, least_squares_gain, compute_gamma,
%    plot_sensor_dose_planes, das_reconstruct (for 'das').
%    Data: RayStationFiles/1194203/Session_1/processed/ (Step 1.5 output:
%    total_dose_CT_1/3.mat, CBCT1/3_resampled.mat, metadata.mat, and the
%    dose_*.mat file names for the MU).
%
%  See also: run_single_field_simulation, step205_run_full_field_simulation,
%            noise_ensemble_error_bars, compute_gamma

clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end
addpath(genpath(fullfile(scriptDir, 'utils')));

%% ========================= 0. CONFIGURATION =============================
%  Simulation settings come from get_default_config: the same defaults the
%  pipeline runs with (lateral sensor, threshold_2 tissue, 4 us rectangular
%  linac pulse, noise at 1/8 of the peak signal, single-pass TR). Only the
%  demo and MBR knobs are added here.

CONFIG = get_default_config();
CONFIG.patient_id    = '1194203';
CONFIG.session       = 'Session_1';

% Gamma criteria (same values as pipeline_simulate.m)
CONFIG.gamma_dose_pct        = 3;
CONFIG.gamma_dist_mm         = 3;
CONFIG.gamma_dose_cutoff_pct = 10;

% MBR knobs, each explained in the step that uses it
CONFIG.mbr_adjoint         = 'tr';    % Step 3:  reverse predictor, 'tr' (k-Wave backwards) | 'das' (delay-and-sum)
CONFIG.mbr_num_iter        = 10;      % Step 9:  loop iterations, identical for both arms
CONFIG.mbr_adjoint_pairs   = 3;       % Step 4:  random image/signal pairs
CONFIG.mbr_power_iter      = 10;      % Step 5:  power iterations
CONFIG.mbr_tv_weight       = 0.005;   % Step 6:  TV strength (fraction of max pressure per iteration)
CONFIG.mbr_tv_edge         = 0.05;    % Step 6:  ripple/edge threshold (fraction of max pressure)
CONFIG.mbr_tune_tv         = false;   % Step 6:  true = pick the weight with the null/shift test (slow)
CONFIG.mbr_tv_candidates   = [0.001, 0.005, 0.02];
CONFIG.mbr_tune_shift_vox  = 2;       % Step 6:  known phantom shift (voxels toward the head)
CONFIG.mbr_null_split_iter = 10;      % Step 11: iterations for the seen/unseen split (0 = skip)

processedDir = fullfile(CONFIG.working_dir, 'RayStationFiles', ...
    CONFIG.patient_id, CONFIG.session, 'processed');

% Place the sensor from the summed CT_1 plan dose so the array sits where it
% would for every field of this plan (determine_sensor_mask_lateral reads it).
CONFIG.total_dose_file = fullfile(processedDir, 'total_dose_CT_1.mat');

% Gamma runs at the very end; fail now rather than after an hour of k-Wave.
if exist('CalcGamma', 'file') ~= 2
    error('demo_mbr:MissingDependency', 'CalcGamma must be on the MATLAB path.');
end
if ~any(strcmp(CONFIG.mbr_adjoint, {'tr', 'das'}))
    error('demo_mbr:BadAdjoint', 'CONFIG.mbr_adjoint must be ''tr'' or ''das'' (got ''%s'').', CONFIG.mbr_adjoint);
end

%% ========================= 1. LOAD THE TOTAL DOSES AND ANATOMY ==========
%  The example is the whole plan, not one beam segment: each CT's total
%  dose, the sum of every segment dose RayStation calculated on that CT
%  (Step 1.5 writes total_dose_CT_1.mat and total_dose_CT_3.mat). The CT_1
%  total is the R1 truth and the CT_3 total is the R3 truth. As in
%  step205_run_full_field_simulation, the plan is treated as ONE field
%  delivered over every pulse of that CT's segments, so the pulse count
%  comes from the summed meterset of those segments.

md = load(fullfile(processedDir, 'metadata.mat'), 'metadata');
beamMetadata = md.metadata.beam_metadata;   % all beams: the sensor must avoid every entrance path
beamNumbers  = [beamMetadata.beam_number];
beamMU       = [beamMetadata.meterset];

arm = struct('name',     {'R1', 'R3'}, ...
             'ct_label', {'CT_1', 'CT_3'}, ...
             'cbct_var', {'CBCT1_resampled', 'CBCT3_resampled'});

for a = 1:2
    totalFile = fullfile(processedDir, ['total_dose_', arm(a).ct_label, '.mat']);
    if ~isfile(totalFile)
        error('demo_mbr:NoTotalDose', 'Per-CT total dose not found: %s (run Step 1.5 first).', totalFile);
    end
    T = load(totalFile);
    if isfield(T, 'ct_total_sparse')
        totalDose = reshape(full(T.ct_total_sparse), T.ct_total_dims);   % stored sparse to save disk
    else
        totalDose = T.ct_total;
    end
    S = load(fullfile(processedDir, [arm(a).cbct_var, '.mat']), arm(a).cbct_var);
    arm(a).sct = S.(arm(a).cbct_var);
    if ~isequal(size(totalDose), size(arm(a).sct.cubeHU))
        error('demo_mbr:GridMismatch', 'Total dose and %s grids differ.', arm(a).cbct_var);
    end
    if ~isfield(arm(a).sct, 'couchMask') || isempty(arm(a).sct.couchMask)
        arm(a).sct.couchMask = false(size(arm(a).sct.bodyMask));
    end
    arm(a).dose_truth = double(totalDose) .* double(arm(a).sct.bodyMask);   % dose only exists in the body

    % Summed meterset of this CT's segments. Step 1.5 gives every segment its
    % beam's RTPLAN meterset, so the beam number in each file name is looked
    % up in metadata.mat: no segment dose file is opened.
    segFiles = dir(fullfile(processedDir, ['dose_*_', arm(a).ct_label, '_B*.mat']));
    if isempty(segFiles)
        error('demo_mbr:NoSegments', 'No %s segment doses in %s', arm(a).ct_label, processedDir);
    end
    arm(a).meterset = 0;
    for i = 1:numel(segFiles)
        beamNum = str2double(regexp(segFiles(i).name, '_B(\d+)_\d+\.mat$', 'tokens', 'once'));
        match   = find(beamNumbers == beamNum, 1);
        if isempty(match)
            error('demo_mbr:NoBeamMU', 'Beam of %s not found in metadata.beam_metadata.', segFiles(i).name);
        end
        arm(a).meterset = arm(a).meterset + beamMU(match);
    end
    fprintf('[DEMO] %s truth: %s total, max %.3f Gy, %.1f MU over %d segments\n', ...
        arm(a).name, arm(a).ct_label, max(arm(a).dose_truth(:)), arm(a).meterset, numel(segFiles));
end
clear T S totalDose

doseDims  = size(arm(1).dose_truth);
spacingMm = arm(1).sct.spacing(:)';
if ~isequal(size(arm(2).dose_truth), doseDims)
    error('demo_mbr:GridMismatch', 'CT_1 and CT_3 doses are on different grids.');
end

fprintf('[DEMO] Grid [%d %d %d], spacing [%.2f %.2f %.2f] mm\n', doseDims, spacingMm);

%% ========================= 2. ACOUSTIC MEDIA AND TRUE PRESSURE ==========
%  CT numbers (HU) become the acoustic properties k-Wave needs: density,
%  sound speed, attenuation, and the Gruneisen coefficient (how much
%  pressure one unit of absorbed dose creates). create_acoustic_medium is
%  the pipeline's single source of truth for this; outside the body we force
%  the water coupling bath, exactly as run_standalone_field does.
%
%  Each arm's TRUE initial pressure is built on its OWN anatomy:
%      p0 = (dose per linac pulse) x Gruneisen x density
%  then smoothed slightly (k-Wave's smooth) to avoid grid-scale ringing, as
%  the pipeline does.

for a = 1:2
    arm(a).medium     = build_medium_with_bath(arm(a).sct, CONFIG);
    arm(a).num_pulses = ceil(arm(a).meterset / CONFIG.dose_per_pulse_cGy);
    p0 = arm(a).dose_truth / arm(a).num_pulses .* arm(a).medium.gruneisen .* arm(a).medium.density;
    arm(a).p0_true = smooth(p0);
end
medium1 = arm(1).medium;   % CT1: the reconstruction anatomy for BOTH arms

%% ========================= 3. SENSOR AND k-WAVE GRID ====================
%  The 2D ultrasound array is pressed against the patient's flank by
%  determine_sensor_mask_lateral, on CT1, avoiding every beam's entrance
%  path. If the array does not fit inside the CT grid the function asks for
%  a water border (grid_pad). Every volume is then padded on the far side to
%  a size whose FFT is fast (k-Wave works in Fourier space). All padding is
%  water with zero pressure.

sctForSensor = arm(1).sct;
if ~isfield(sctForSensor, 'origin') || isempty(sctForSensor.origin)
    sctForSensor.origin = [0, 0, 0];
end
sctForSensor.spacing = spacingMm;
doseForSensor = struct('dose_Gy', arm(1).dose_truth, ...   % no gantry angle: the whole plan
    'origin', sctForSensor.origin, 'spacing', spacingMm, 'dimensions', doseDims);
[sensorMask, sensorInfo] = determine_sensor_mask_lateral(sctForSensor, doseForSensor, beamMetadata, CONFIG);
gp = sensorInfo.grid_pad;

offset  = [gp.y_pre, gp.x_pre, gp.z_pre];   % where the dose grid starts inside the k-Wave grid
expDims = doseDims + offset + [gp.y_post, gp.x_post, gp.z_post];
kDims   = zeros(1, 3);
for d = 1:3
    kDims(d) = find_optimal_kwave_size(expDims(d), CONFIG.pml_size);
end
sensorGrid = padarray(double(sensorMask), kDims - size(sensorMask), 0, 'post');

for a = 1:2
    arm(a).kmedium.density     = to_kwave_grid(arm(a).medium.density,     gp, kDims, 1000);
    arm(a).kmedium.sound_speed = to_kwave_grid(arm(a).medium.sound_speed, gp, kDims, 1540);
    arm(a).kmedium.alpha_coeff = to_kwave_grid(arm(a).medium.alpha_coeff, gp, kDims, 0);
    arm(a).kmedium.alpha_power = 1.1;   % same exponent as run_single_field_simulation
    arm(a).p0_grid = to_kwave_grid(arm(a).p0_true, gp, kDims, 0);
end

% Support: where p0 is allowed to be non-zero. Dose only exists inside the
% body and air pockets make no sound (Gruneisen 0), so the unknowns are the
% non-air body voxels of CT1, the anatomy a live system knows.
bodyGrid = to_kwave_grid(double(arm(1).sct.bodyMask), gp, kDims, 0) > 0.5;
support  = double(bodyGrid & arm(1).kmedium.density >= 100);

% Time stepping, exactly as run_single_field_simulation: the step is set by
% the fastest sound speed (stability), the record length by the slowest
% (2.5 grid diagonals), over BOTH anatomies since both are simulated.
dxM   = spacingMm * 1e-3;
allC  = [arm(1).kmedium.sound_speed(:); arm(2).kmedium.sound_speed(:)];
cMin  = min(allC(allC > 0));
dt    = CONFIG.cfl_number * min(dxM) / max(allC);
Nt    = ceil(2.5 * norm(kDims .* dxM) / cMin / dt);
if CONFIG.Nt_scaling ~= 0 && cMin <= CONFIG.tissue_tables.threshold_2.air.sound_speed
    Nt = ceil(Nt / CONFIG.Nt_scaling);   % air sets cMin and would inflate Nt
end
kgrid    = kWaveGrid(kDims(1), dxM(1), kDims(2), dxM(2), kDims(3), dxM(3));
kgrid.dt = dt;
kgrid.Nt = Nt;

dataCast = 'single';
if CONFIG.use_gpu
    try
        gpuDevice;
        dataCast = 'gpuArray-single';
    catch
        fprintf('[DEMO] No GPU found: running on CPU (much slower).\n');
    end
end

% M holds everything the predictor and reverse predictor need. Its medium is
% CT1 because that is what we reconstruct on.
M = struct();
M.kgrid         = kgrid;
M.kmedium       = arm(1).kmedium;
M.sensor        = struct('mask', sensorGrid);
M.inputArgs     = {'Smooth', false, 'PMLInside', false, 'PMLSize', CONFIG.pml_size, ...
                   'DataCast', dataCast, 'PlotSim', false};
M.H             = fft(build_pulse_kernel(CONFIG.pulse_shape, CONFIG.convolution_kernel, dt), Nt).';
M.FS            = 1 / dt;
M.adjoint_scale = 1;   % calibrated in Step 4

% Delay-and-sum reverse predictor (used when CONFIG.mbr_adjoint = 'das').
% Pipeline DAS (das_reconstruct) with four settings changed so it can stand
% in for the reverse wave run:
%   envelope / clip_negative off : both are non-linear; a blame map must keep
%       its sign (push DOWN where the guess predicts too much signal)
%   use_aperture off : the predictor's sensor points hear every direction, so
%       the back-projection must reach every direction too
%   interp linear    : smoother blame than taking the nearest time sample
% Element averaging stays on: back-project from the 32x32 element centres
% instead of every sensor voxel (far fewer channels).
M.adjoint     = CONFIG.mbr_adjoint;
M.sensor_info = sensorInfo;
M.dx          = dxM;
M.das_opts    = struct('use_elements', true, 'envelope', false, 'clip_negative', false, ...
                       'use_aperture', false, 'depth_weight', true, 'interp', 'linear');

% Pressure -> dose conversion data (Step 10): CT1 Gruneisen x density.
G = struct('offset', offset, 'dims', doseDims);
G.gamma_rho = medium1.gruneisen .* medium1.density;
G.gamma_rho(G.gamma_rho == 0) = 1;   % same guard as the pipeline

% Final masks shared by TR and MBR (pipeline convention): own body, own
% truth dose region (>1% of max), no CT1 air.
for a = 1:2
    inDoseRegion = arm(a).dose_truth > 0.01 * max(arm(a).dose_truth(:));
    arm(a).final_mask = double(arm(a).sct.bodyMask & inDoseRegion & medium1.density >= 100);
end

fprintf('[DEMO] k-Wave grid [%d %d %d], %d sensor points, dt = %.2e s, Nt = %d\n', ...
    kDims, nnz(sensorGrid), dt, Nt);

dose1Grid = to_kwave_grid(arm(1).dose_truth, gp, kDims, 0);
plot_sensor_dose_planes(dose1Grid >= 0.1 * max(dose1Grid(:)), sensorGrid > 0, spacingMm, ...
    arm(1).kmedium.density, CONFIG);
clear dose1Grid

%% ========================= STEP 1: RAW SENSOR SIGNALS ===================
%  This is the "experiment". Each arm's true p0 launches sound that travels
%  through that arm's TRUE anatomy (CT1 for R1, CT3 for R3) to the sensor.
%  The recorded traces then go through what the hardware does to them, in
%  acquisition order (same chain as run_single_field_simulation):
%    a. the linac pulse lasts a few us, so every trace is smeared by the
%       pulse shape
%    b. the transducer only hears a band of frequencies around 0.35 MHz
%    c. electronic noise is added
%  a + b are exactly what the predictor (Step 2) does, so the noise-free
%  signals come from calling the predictor with the TRUE anatomy.
%  Everything below starts from these traces, NOT from a TR image.

for a = 1:2
    Mtrue = M;
    Mtrue.kmedium = arm(a).kmedium;   % the anatomy the waves really crossed
    arm(a).clean  = predict_signals(arm(a).p0_grid, Mtrue);

    if ~isempty(CONFIG.noise_amp_Pa) && CONFIG.noise_amp_Pa > 0
        arm(a).noise_amp = CONFIG.noise_amp_Pa;                                 % calibrated hardware noise
    else
        arm(a).noise_amp = CONFIG.conv_noise_level * max(abs(arm(a).clean(:))); % fraction of this dose's peak
    end
    noise = arm(a).noise_amp * randn(size(arm(a).clean));
    arm(a).y = arm(a).clean + noise;
    arm(a).noise_floor = norm(noise(:)) / norm(arm(a).y(:));   % used in Step 11
    fprintf('[MBR step 1] %s: peak signal %.3e Pa, noise %.3e Pa, SNR %.1f\n', arm(a).name, ...
        max(abs(arm(a).clean(:))), arm(a).noise_amp, max(abs(arm(a).clean(:))) / arm(a).noise_amp);
end
clear noise Mtrue

timeUs = (0:Nt - 1) * dt * 1e6;
[~, loudest] = max(max(abs(arm(1).clean), [], 2));
sigLim = max(abs(arm(1).clean(:))) * [-1, 1];
figure('Name', 'Step 1: raw sensor signals', 'Color', 'w', 'Position', [60, 60, 1400, 420]);
subplot(1, 3, 1);
imagesc(timeUs, 1:size(arm(1).clean, 1), arm(1).clean); caxis(sigLim); colormap(gca, diverging_cmap());
xlabel('time (\mus)'); ylabel('sensor point'); title('R1 noise-free signals'); colorbar;
subplot(1, 3, 2);
imagesc(timeUs, 1:size(arm(1).y, 1), arm(1).y); caxis(sigLim); colormap(gca, diverging_cmap());
xlabel('time (\mus)'); ylabel('sensor point'); title('R1 measured (with noise)'); colorbar;
subplot(1, 3, 3);
plot(timeUs, arm(1).y(loudest, :), 'Color', [0.7 0.7 0.7]); hold on;
plot(timeUs, arm(1).clean(loudest, :), 'b', 'LineWidth', 1.4);
plot(timeUs, arm(2).clean(loudest, :), 'r', 'LineWidth', 1.4); hold off;
xlabel('time (\mus)'); ylabel('pressure (Pa)'); grid on;
legend('R1 measured', 'R1 noise-free', 'R3 noise-free (CT3 anatomy)');
title(sprintf('Loudest sensor point (#%d)', loudest));

%% ========================= BASELINE: TIME REVERSAL ======================
%  What the pipeline does today (run_single_field_simulation, method 'tr'),
%  on the CT1 anatomy for both arms:
%    1. undo the pulse smear with a Wiener deconvolution
%    2. play the time-flipped traces back from the sensor points, FORCING
%       the pressure there to the recorded values, and keep the pressure
%       field left at the end
%  The result is both our comparison and the MBR starting point (Step 7).

for a = 1:2
    fprintf('[BASELINE] time reversal, %s\n', arm(a).name);
    arm(a).p_tr = time_reversal(arm(a).y, M, CONFIG);
end

%% ========================= STEPS 2 + 3: PREDICTOR AND REVERSE PREDICTOR ==
%  Both live at the bottom of this file; read them now.
%
%  predict_signals (A): any guessed pressure image -> the signals the sensor
%    would record from it. k-Wave on the CT1 medium, then the pulse smear,
%    then the transducer response. Rule: anything the real system does to a
%    signal belongs in the predictor; whatever is left out, the
%    reconstruction will wrongly blame on the image.
%
%  reverse_predict (A'): a set of signals -> a "blame map". For a signal
%    mismatch it says which voxels are responsible, and by how much. It runs
%    each predictor stage backwards in reverse order. Mathematically A' is
%    the adjoint (transpose) of A; Step 4 checks this. The wave part has two
%    options (CONFIG.mbr_adjoint):
%    'tr'  a k-Wave run that looks like TR but differs in two ways:
%      * the flipped signals are ADDED at the sensor points (TR overwrites
%        the pressure there)
%      * absorption is not flipped: sound still loses energy going backwards
%    'das' delay-and-sum: each voxel collects every element's signal at the
%        travel time from that voxel to the element. Straight rays, one
%        average sound speed, no absorption, no anatomy: an APPROXIMATE
%        opposite of the predictor (a "mismatched" pair). The loop can still
%        converge, but to a slightly different image than with 'tr'.

%% ========================= STEP 4: ARE THEY TRUE OPPOSITES? =============
%  For a random image x and random signals r, an exact pair satisfies
%      sum( A(x) .* r )  ==  sum( x .* A'(r) )
%  (signal side vs image side: the same number computed two ways). k-Wave
%  rescales "additive" sources internally by a constant that depends on dt,
%  dx and the sound speed at the sensor, so we expect the two sides to differ
%  by a CONSTANT factor, which we fold into A'. If the factor changes from one
%  random pair to the next, A' is not A's opposite (wrong time alignment, a
%  missing stage, ...) and the loop below cannot be trusted.
%  With 'das' a warning is EXPECTED: DAS is only an approximate opposite, and
%  the spread printed here measures how approximate. What matters then is
%  that the ratio keeps the same sign (a wrong sign would push the image
%  uphill) and that the Step 11 mismatch curve still falls.

ratios = zeros(CONFIG.mbr_adjoint_pairs, 1);
for k = 1:CONFIG.mbr_adjoint_pairs
    xRand = randn(kDims) .* support;
    rRand = randn(size(arm(1).y));
    signalSide = sum(predict_signals(xRand, M) .* rRand, 'all');
    imageSide  = sum(xRand .* reverse_predict(rRand, M), 'all');
    ratios(k)  = signalSide / imageSide;
    fprintf('[MBR step 4] pair %d: signal side %.4e, image side %.4e, ratio %.6g\n', ...
        k, signalSide, imageSide, ratios(k));
end
M.adjoint_scale = mean(ratios);
spread = std(ratios) / abs(mean(ratios));
fprintf('[MBR step 4] reverse predictor scaled by %.6g (spread between pairs %.2f%%)\n', ...
    M.adjoint_scale, 100 * spread);
if spread > 0.05
    warning('demo_mbr:AdjointMismatch', ...
        'Adjoint ratio varies by %.1f%% between random pairs: A'' is not A''s opposite.', 100 * spread);
end
if any(sign(ratios) ~= sign(M.adjoint_scale))
    warning('demo_mbr:AdjointSign', ...
        'Adjoint ratio changes sign between pairs: the reverse predictor does not reliably point downhill.');
end
clear xRand rRand

%% ========================= STEP 5: SAFE STEP SIZE =======================
%  One loop iteration moves the image by step x A'(A(image) - y). If the step
%  is too big, an error gets amplified on every round trip and the loop blows
%  up. The most a round trip A'(A(x)) can amplify any image is a single
%  number L. Find it by power iteration: push a random image through A then
%  A' over and over, rescaling to length 1 each time. Whatever the image
%  starts as, it turns into the most-amplified pattern, and the growth per
%  round trip settles to L. Safe steps are below about 1/L. (With 'das' the
%  round trip is not perfectly symmetric, so L may wobble rather than settle;
%  the last value is used.)

x = randn(kDims) .* support;
x = x / norm(x(:));
L_history = zeros(CONFIG.mbr_power_iter, 1);
for k = 1:CONFIG.mbr_power_iter
    x = reverse_predict(predict_signals(x, M), M) .* support;
    L_history(k) = norm(x(:));   % growth of a length-1 image in one round trip
    x = x / L_history(k);
    fprintf('[MBR step 5] power iteration %d/%d: L = %.6g\n', k, CONFIG.mbr_power_iter, L_history(k));
end
L = L_history(end);
clear x
% Steps 2-5 depend only on CT1 and the sensor. In production, adjoint_scale
% and L would be computed once per patient and reused for every field and
% both arms.

figure('Name', 'Steps 4-5: adjoint check and step size', 'Color', 'w', 'Position', [80, 80, 1000, 380]);
subplot(1, 2, 1);
bar(ratios); hold on; yline(M.adjoint_scale, 'r--', 'mean'); hold off;
xlabel('random pair'); ylabel('signal side / image side');
title(sprintf('Step 4: adjoint ratio (spread %.2f%%)', 100 * spread));
subplot(1, 2, 2);
plot(1:CONFIG.mbr_power_iter, L_history, '-o', 'LineWidth', 1.5); grid on;
xlabel('power iteration'); ylabel('round-trip amplification');
title(sprintf('Step 5: L settles to %.3g  ->  step < 1/L', L));

%% ========================= STEP 6: THE PRIOR (TV) AND ITS WEIGHT ========
%  Noise and the one-sided view leave speckle and streaks the data cannot
%  rule out. The prior is a mild preference that breaks those ties. We use
%  total variation (TV): the sum over voxels of how much the image jumps
%  between neighbours. A noisy ripple jumps up and down many times -> large
%  TV -> penalized. A sharp field edge is ONE jump -> left alone. (A plain
%  blur would smear the edge too.) Two dimensionless knobs:
%    mbr_tv_weight : how hard TV pushes. One iteration of TV alone moves a
%                    voxel by about tv_weight x the max pressure.
%    mbr_tv_edge   : neighbour-to-neighbour jumps smaller than this fraction
%                    of the max pressure count as ripples (smoothed); larger
%                    jumps count as edges (kept).
%
%  Choosing the weight (CONFIG.mbr_tune_tv = true; costly):
%    * null test: the SAME R1 dose measured twice with different noise. Any
%      difference between the two reconstructions is "no change".
%    * known change: the R1 dose shifted by mbr_tune_shift_vox voxels.
%    For each candidate, separation = size of the real change / size of the
%    no-change difference. Too little TV: noise swamps the change. Too much:
%    TV flattens the change itself. Keep the weight that separates best.

if CONFIG.mbr_tune_tv
    nullB   = arm(1).clean + arm(1).noise_amp * randn(size(arm(1).clean));
    shifted = predict_signals(circshift(arm(1).p0_grid, CONFIG.mbr_tune_shift_vox, 3), M);   % +z = toward the head
    shifted = shifted + arm(1).noise_amp * randn(size(shifted));
    tuneData  = {arm(1).y, nullB, shifted};
    tuneStart = cell(1, 3);
    for k = 1:3
        tuneStart{k} = scale_to_signals(time_reversal(tuneData{k}, M, CONFIG), tuneData{k}, M, support);
    end

    separation = zeros(size(CONFIG.mbr_tv_candidates));
    for w = 1:numel(CONFIG.mbr_tv_candidates)
        opt   = mbr_options(CONFIG.mbr_tv_candidates(w), CONFIG, L, max(tuneStart{1}(:)));
        recon = cell(1, 3);
        for k = 1:3
            recon{k} = run_mbr(tuneData{k}, tuneStart{k}, M, support, opt, ...
                sprintf('tune %g, set %d', CONFIG.mbr_tv_candidates(w), k));
        end
        noChange   = recon{2} - recon{1};
        realChange = recon{3} - recon{1};
        separation(w) = norm(realChange(:)) / norm(noChange(:));
        fprintf('[MBR step 6] weight %g: separation %.3f\n', CONFIG.mbr_tv_candidates(w), separation(w));
    end
    [~, best] = max(separation);
    CONFIG.mbr_tv_weight = CONFIG.mbr_tv_candidates(best);
    clear nullB shifted tuneData tuneStart recon noChange realChange

    figure('Name', 'Step 6: TV weight', 'Color', 'w');
    semilogx(CONFIG.mbr_tv_candidates, separation, '-o', 'LineWidth', 1.5); grid on;
    xlabel('TV weight'); ylabel('real change / no-change difference');
    title(sprintf('Best weight: %g', CONFIG.mbr_tv_weight));
end
fprintf('[MBR step 6] TV weight %g, edge threshold %g x max pressure\n', ...
    CONFIG.mbr_tv_weight, CONFIG.mbr_tv_edge);

%% ========================= STEPS 7-9: THE RECONSTRUCTION LOOP ===========
%  7. Start from the TR image with negatives set to zero, scaled so its
%     predicted signals best match the measurement (TR gets the shape roughly
%     right but not the overall size). Starting close saves iterations.
%  8. Repeat (see run_mbr at the bottom of this file):
%       a. Predict  - run the predictor on the current guess
%       b. Compare  - subtract the measured signals: the mismatch
%       c. Blame    - reverse predictor on the mismatch: which way each
%                     voxel should move to shrink it
%       d. Prior    - add TV's push (flatten ripples, keep edges)
%       e. Step     - move the guess against the combined push
%       f. Clamp    - negative pressure is impossible, and there is no dose
%                     outside the body
%       g. Coast    - add part of the previous move, like a ball rolling
%                     downhill that keeps some speed (this is what makes
%                     FISTA faster than plain stepping)
%  9. Stop after a FIXED number of iterations, the same for R1 and R3. A
%     "stop when the mismatch is small" rule would let the arms stop at
%     different points, and the iteration count itself changes how smooth
%     the image is, which would show up as a fake R1-R3 difference.
%  Both arms use the SAME objective (same lambda, same edge threshold), with
%  the pressure scale taken from R1's starting guess.

for a = 1:2
    arm(a).x0 = scale_to_signals(arm(a).p_tr, arm(a).y, M, support);
end
opt = mbr_options(CONFIG.mbr_tv_weight, CONFIG, L, max(arm(1).x0(:)));
fprintf('[MBR steps 7-9] step %.4g (= %.2f / L), lambda %.4g, %d iterations per arm\n', ...
    opt.step, opt.step * L, opt.lambda, opt.num_iter);
for a = 1:2
    [arm(a).p_mbr, arm(a).mismatch] = run_mbr(arm(a).y, arm(a).x0, M, support, opt, arm(a).name);
end

%% ========================= STEP 10: PRESSURE -> DOSE ====================
%  p0 = dose x Gruneisen x density, so dose = p0 / (Gruneisen x density)
%  with the CT1 maps (the anatomy we reconstructed on), times the number of
%  pulses. TR and MBR get identical treatment, including the pipeline's final
%  masks. TR also gets the pipeline's empirical correction_factor; it does
%  not affect the gamma comparison because both reconstructions are
%  least-squares scaled to the truth before gamma (pipeline convention).

for a = 1:2
    arm(a).dose_tr  = pressure_to_dose(arm(a).p_tr * CONFIG.correction_factor, G, arm(a).num_pulses) ...
                      .* arm(a).final_mask;
    arm(a).dose_mbr = pressure_to_dose(arm(a).p_mbr, G, arm(a).num_pulses) .* arm(a).final_mask;
end

%% ========================= STEP 11a: MISMATCH VS ITERATION ==============
%  The relative mismatch ||A p - y|| / ||y|| should fall and then flatten.
%  The dashed lines are the size of the noise we added: no image can beat
%  that. R1 (exact model) should flatten near its noise line. R3 flattens
%  higher: the extra gap is the model-error floor from reconstructing CT3
%  data with the CT1 anatomy.

figure('Name', 'Step 11: mismatch vs iteration', 'Color', 'w', 'Position', [100, 100, 700, 450]);
colors = lines(2);
hold on;
for a = 1:2
    plot(1:CONFIG.mbr_num_iter, arm(a).mismatch, '-o', 'Color', colors(a, :), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('%s mismatch', arm(a).name));
    yline(arm(a).noise_floor, '--', 'Color', colors(a, :), 'LineWidth', 1.2, ...
        'DisplayName', sprintf('%s noise level', arm(a).name));
end
hold off; grid on; legend('Location', 'northeast');
xlabel('iteration'); ylabel('||predicted - measured|| / ||measured||');
title('Mismatch flattens at the noise (R1) or model-error (R3) floor');

%% ========================= RESULT: GAMMA ANALYSIS =======================
%  Gamma combines a dose-difference test and a distance-to-agreement test:
%  a voxel passes (gamma <= 1) if a reconstructed dose within
%  gamma_dist_mm of it matches the truth to gamma_dose_pct of the max dose.
%  Evaluated where the truth exceeds 10% of its max. Each reconstruction is
%  first scaled by its least-squares gain to its own truth (pipeline
%  convention), so we compare SHAPE fidelity, not the arbitrary TR scale.

crit   = {CONFIG.gamma_dose_pct, CONFIG.gamma_dist_mm, ...
          sprintf('%g%%/%gmm', CONFIG.gamma_dose_pct, CONFIG.gamma_dist_mm)};
cutoff = CONFIG.gamma_dose_cutoff_pct / 100;
for a = 1:2
    truth = arm(a).dose_truth;
    arm(a).dose_tr_scaled  = least_squares_gain(truth, arm(a).dose_tr)  * arm(a).dose_tr;
    arm(a).dose_mbr_scaled = least_squares_gain(truth, arm(a).dose_mbr) * arm(a).dose_mbr;
    fprintf('[GAMMA] %s time reversal\n', arm(a).name);
    arm(a).gamma_tr  = compute_gamma(truth, arm(a).dose_tr_scaled,  spacingMm, 'Criteria', crit, 'Cutoff', cutoff);
    fprintf('[GAMMA] %s MBR\n', arm(a).name);
    arm(a).gamma_mbr = compute_gamma(truth, arm(a).dose_mbr_scaled, spacingMm, 'Criteria', crit, 'Cutoff', cutoff);
end

fprintf('\n=============== GAMMA %s PASS RATE ===============\n', crit{3});
fprintf('  Arm    TR         MBR        MBR - TR\n');
mbrWins = true;
for a = 1:2
    passTr  = arm(a).gamma_tr.pass_rates(1);
    passMbr = arm(a).gamma_mbr.pass_rates(1);
    fprintf('  %s     %6.2f%%    %6.2f%%    %+6.2f\n', arm(a).name, passTr, passMbr, passMbr - passTr);
    mbrWins = mbrWins && passMbr > passTr;
end
if mbrWins
    fprintf('  CONFIRMED: MBR is closer to the truth than TR in both arms.\n');
else
    fprintf('  NOT confirmed: MBR did not beat TR in every arm (see table).\n');
end
fprintf('==========================================================\n\n');

%% ========================= RESULT: IMAGES ===============================
%  Rows: truth, TR, MBR (both scaled to the truth as in the gamma analysis).
%  Columns: transverse, coronal, sagittal through the truth maximum, body
%  outline in cyan, one colour scale per figure.

views    = {'Transverse', 'Coronal', 'Sagittal'};
rowNames = {'Truth', 'Time reversal', 'MBR'};
for a = 1:2
    [~, iMax] = max(arm(a).dose_truth(:));
    [iy, ix, iz] = ind2sub(doseDims, iMax);
    arm(a).ijk = [iy, ix, iz];
    vols     = {arm(a).dose_truth, arm(a).dose_tr_scaled, arm(a).dose_mbr_scaled};
    doseClim = [0, max(arm(a).dose_truth(:))];

    figure('Name', sprintf('%s: truth vs TR vs MBR', arm(a).name), 'Color', 'w', ...
        'Position', [40, 40, 1300, 950]);
    for r = 1:3
        for c = 1:3
            ax = subplot(3, 3, (r - 1) * 3 + c);
            show_view(ax, vols{r}, arm(a).sct.bodyMask, views{c}, arm(a).ijk, spacingMm, doseClim, hot(256));
            title(ax, sprintf('%s - %s', rowNames{r}, views{c}));
        end
        colorbar(ax);
    end
    sgtitle(sprintf('%s (%s total): dose in Gy  |  gamma %s pass: TR %.1f%%, MBR %.1f%%', ...
        arm(a).name, strrep(arm(a).ct_label, '_', '\_'), crit{3}, ...
        arm(a).gamma_tr.pass_rates(1), arm(a).gamma_mbr.pass_rates(1)));
end

% Gamma maps: blue = pass, white = gamma 1, red = fail; blank outside the
% evaluated (>=10%) region.
figure('Name', 'Gamma maps', 'Color', 'w', 'Position', [60, 60, 1000, 850]);
methodNames = {'TR', 'MBR'};
for a = 1:2
    results = {arm(a).gamma_tr, arm(a).gamma_mbr};
    for m = 1:2
        gammaMap = results{m}.maps{1};
        gammaMap(~results{m}.eval_mask) = NaN;
        ax = subplot(2, 2, (m - 1) * 2 + a);
        show_view(ax, gammaMap, arm(a).sct.bodyMask, 'Transverse', arm(a).ijk, spacingMm, [0, 2], diverging_cmap());
        title(ax, sprintf('%s %s: pass %.1f%%', arm(a).name, methodNames{m}, results{m}.pass_rates(1)));
        colorbar(ax);
    end
end
sgtitle(sprintf('Gamma %s, transverse slice through the truth maximum', crit{3}));

%% ========================= STEP 11b: WHAT THE SENSOR CAN AND CANNOT SEE ==
%  Some images make NO signal at the sensor: structures whose sound never
%  reaches a one-sided array, or detail finer than the transducer band.
%  Adding them to an image does not change the predicted signals at all.
%  Following Bhadra et al. (IEEE TMI 2021, "On hallucinations in tomographic
%  image reconstruction"), split the R1 - R3 difference into
%    seen   : the part the sensor responds to (the data pin it down)
%    unseen : the rest. The data cannot have put it there, so any unseen
%             difference came from the prior (TV, positivity) or the TR
%             starting guess, NOT from a change the system measured.
%  The seen part is the smallest image that reproduces the difference's
%  signals. The plain predict/blame loop (no TV, no clamp) started from zero
%  converges to exactly that image, because every update is a
%  reverse-predictor output and those only contain seen content. A few
%  iterations capture only the strongly-seen part, so the unseen fraction
%  printed below is an upper bound. This relies on the reverse predictor
%  being the true opposite of A, so with mbr_adjoint = 'das' the split is
%  only approximate.

if CONFIG.mbr_null_split_iter > 0
    diffP  = arm(1).p_mbr - arm(2).p_mbr;   % R1 - R3 (pressure, k-Wave grid)
    target = predict_signals(diffP, M);     % what the sensor records from the difference
    seen   = zeros(kDims);
    stepLs = 1 / (1.1 * L);
    for k = 1:CONFIG.mbr_null_split_iter
        if k == 1
            resid = -target;                % an all-zero image predicts zero signals: skip that run
        else
            resid = predict_signals(seen, M) - target;
        end
        seen = seen - stepLs * reverse_predict(resid, M) .* support;
    end
    unseen = diffP - seen;
    unseenFraction = norm(unseen(:))^2 / norm(diffP(:))^2;
    fprintf(['[MBR step 11] %.1f%% of the R1 - R3 difference (energy) is in the part the sensor ' ...
             'cannot see: it came from the prior, not the data.\n'], 100 * unseenFraction);

    % Display in dose units (linear conversion, R1 pulse count)
    parts     = {diffP, seen, unseen};
    partNames = {'R1 - R3 (MBR)', 'seen by the sensor', 'unseen (from the prior)'};
    partDose  = cell(1, 3);
    for k = 1:3
        partDose{k} = pressure_to_dose(parts{k}, G, arm(1).num_pulses) .* double(arm(1).sct.bodyMask);
    end
    diffLim = max(abs(partDose{1}(:))) * [-1, 1];
    figure('Name', 'Step 11: seen / unseen split', 'Color', 'w', 'Position', [80, 80, 1400, 420]);
    for k = 1:3
        ax = subplot(1, 3, k);
        show_view(ax, partDose{k}, arm(1).sct.bodyMask, 'Transverse', arm(1).ijk, spacingMm, diffLim, diverging_cmap());
        title(ax, partNames{k});
        colorbar(ax);
    end
    sgtitle(sprintf('R1 - R3 split (Gy): %.1f%% of the difference is unseen', 100 * unseenFraction));
    clear diffP target seen unseen resid parts
end

fprintf('[DEMO] Done.\n');


%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

function signals = predict_signals(p0, M)
%PREDICT_SIGNALS The predictor A (Step 2): pressure image -> recorded signals.
%   1. k-Wave: p0 launches sound through M.kmedium; record the pressure at
%      every sensor point (one row per point, one column per time step)
%   2. smear by the linac pulse (multiply by the pulse spectrum M.H)
%   3. transducer frequency response: 0.35 MHz centre, 100% bandwidth, the
%      same filter as run_single_field_simulation
    source.p0 = p0;
    raw = kspaceFirstOrder3D(M.kgrid, M.kmedium, source, M.sensor, M.inputArgs{:});
    raw = double(gather(raw));
    smeared = real(ifft(fft(raw, [], 2) .* M.H, [], 2));
    signals = gaussianFilter(smeared, M.FS, 0.35e6, 100, false);
end


function img = reverse_predict(signals, M)
%REVERSE_PREDICT The reverse predictor A' (Step 3): signals -> blame map.
%  Runs the stages of predict_signals backwards, each replaced by its mirror
%  image (adjoint):
%   3'. transducer filter: symmetric in frequency, so it is its own mirror
%   2'. pulse smear: correlate with the pulse (conj(M.H)) instead of
%       convolving
%   1'. wave propagation backwards, by M.adjoint:
%       'tr'  k-Wave: inject the time-flipped signals at the sensor points as
%             ADDED sources (TR would overwrite the pressure there), propagate
%             with normal absorption, record the pressure everywhere at the end
%       'das' delay-and-sum on straight rays (das_reconstruct); no k-Wave run
%  M.adjoint_scale is the constant found in Step 4.
    filtered  = gaussianFilter(signals, M.FS, 0.35e6, 100, false);
    unsmeared = real(ifft(fft(filtered, [], 2) .* conj(M.H), [], 2));
    gridSize  = size(M.sensor.mask);
    if strcmp(M.adjoint, 'das')
        % A sensor records the TIME DERIVATIVE of the pressure averaged over a
        % sphere around it, so the matching back-projection sums -d/dt of the
        % signals. Plain DAS of these band-limited signals would be 90 degrees
        % out of phase with the true blame.
        [dSdt, ~] = gradient(unsmeared, 1 / M.FS);
        img = das_reconstruct(-dSdt, M.sensor, M.sensor_info, M.kmedium, ...
            gridSize(1), gridSize(2), gridSize(3), M.dx(1), M.dx(2), M.dx(3), 1 / M.FS, M.das_opts);
    else
        source.p_mask   = M.sensor.mask;
        source.p        = fliplr(unsmeared);
        source.p_mode   = 'additive';
        recorder.mask   = ones(gridSize);
        recorder.record = {'p_final'};
        out = kspaceFirstOrder3D(M.kgrid, M.kmedium, source, recorder, M.inputArgs{:});
        img = reshape(double(gather(out.p_final)), gridSize);
    end
    img = M.adjoint_scale * img;
end


function pTr = time_reversal(y, M, CONFIG)
%TIME_REVERSAL The pipeline's reconstruction ('tr' branch of
%  run_single_field_simulation) on the medium in M:
%   1. Wiener-deconvolve the pulse smear
%   2. play the time-flipped traces back from the sensor, FORCING the
%      pressure there to the recorded values (Dirichlet); keep the final
%      pressure field, negatives set to zero
%   3. only if num_time_reversal_iter > 1: re-simulate the image, add back
%      what it fails to explain, repeat (stops early if the image settles)
    yDeconv  = real(ifft(fft(y, [], 2) .* conj(M.H) ./ (abs(M.H).^2 + CONFIG.conv_deconv_lambda), [], 2));
    gridSize = size(M.sensor.mask);
    playback = yDeconv;
    pPrev    = zeros(gridSize);
    for it = 1:CONFIG.num_time_reversal_iter
        source.p_mask   = M.sensor.mask;
        source.p        = fliplr(playback);
        source.p_mode   = 'dirichlet';
        recorder.mask   = ones(gridSize);
        recorder.record = {'p_final'};
        out = kspaceFirstOrder3D(M.kgrid, M.kmedium, source, recorder, M.inputArgs{:});
        pTr = max(reshape(double(gather(out.p_final)), gridSize), 0);
        if it > 1 && norm(pTr(:) - pPrev(:)) < CONFIG.convergence_tol * norm(pPrev(:))
            break;
        end
        pPrev = pTr;
        if it < CONFIG.num_time_reversal_iter
            resim = kspaceFirstOrder3D(M.kgrid, M.kmedium, struct('p0', pTr), M.sensor, M.inputArgs{:});
            playback = playback + (yDeconv - double(gather(resim)));
        end
    end
end


function x0 = scale_to_signals(pTr, y, M, support)
%SCALE_TO_SIGNALS MBR starting guess (Step 7) from a TR image.
%  Clamp negatives to zero and keep only the support. Then pick the one
%  number g that makes the predicted signals g*A(x0) best match the measured
%  y (least squares). The predictor is linear, so this costs one run.
    x0   = max(pTr, 0) .* support;
    pred = predict_signals(x0, M);
    g    = sum(pred(:) .* y(:)) / sum(pred(:) .^ 2);
    x0   = g * x0;
end


function opt = mbr_options(tvWeight, CONFIG, L, pressureScale)
%MBR_OPTIONS Turn the dimensionless TV knobs (Step 6) into absolute numbers.
%  lambda : tvWeight * L * pressureScale. With a step of 1/L, one iteration
%           of TV alone moves a voxel by about tvWeight x pressureScale.
%  edge   : CONFIG.mbr_tv_edge * pressureScale (ripple/edge threshold, Pa).
%  step   : 1 over the stiffness of the WHOLE objective. The data part has
%           stiffness L (Step 5); smoothed TV adds at most 12*lambda/edge
%           (4 per direction in 3D). The extra 10% covers power iteration
%           approaching L from below.
    opt.lambda   = tvWeight * L * pressureScale;
    opt.edge     = CONFIG.mbr_tv_edge * pressureScale;
    opt.step     = 1 / (1.1 * (L + 12 * opt.lambda / opt.edge));
    opt.num_iter = CONFIG.mbr_num_iter;
end


function [x, mismatch] = run_mbr(y, x0, M, support, opt, label)
%RUN_MBR The reconstruction loop (Step 8 a-g), opt.num_iter iterations.
%  x : current image (always >= 0 and inside the support)
%  z : look-ahead point = x plus part of the last move ("coasting")
%  mismatch(k) is measured at z, the image the predictor actually ran on.
    x = x0;  xPrev = x0;  z = x0;  t = 1;
    mismatch = zeros(opt.num_iter, 1);
    for k = 1:opt.num_iter
        predicted   = predict_signals(z, M);                    % a. predict
        resid       = predicted - y;                            % b. compare
        mismatch(k) = norm(resid(:)) / norm(y(:));
        dataPush    = reverse_predict(resid, M);                % c. blame
        tvPush      = opt.lambda * tv_gradient(z, opt.edge);    % d. prior
        x = z - opt.step * (dataPush + tvPush);                 % e. step
        x = max(x, 0) .* support;                               % f. clamp
        tNext = (1 + sqrt(1 + 4 * t^2)) / 2;                    % g. coast
        z = x + ((t - 1) / tNext) * (x - xPrev);
        xPrev = x;
        t = tNext;
        fprintf('[MBR %s] iteration %d/%d: relative mismatch %.4f\n', label, k, opt.num_iter, mismatch(k));
    end
end


function g = tv_gradient(p, edge)
%TV_GRADIENT Gradient of smoothed total variation,
%  TV(p) = sum over voxels of sqrt(|jump to neighbours|^2 + edge^2).
%  Where jumps are much smaller than `edge` (ripples) this acts like a gentle
%  smoothing: bumps are pushed down, dips up. Where a jump is much larger
%  (an edge) the push has a fixed size however big the edge is, so edges are
%  not smeared the way a blur would smear them.
    jx  = fwd_diff(p, 1);
    jy  = fwd_diff(p, 2);
    jz  = fwd_diff(p, 3);
    mag = sqrt(jx.^2 + jy.^2 + jz.^2 + edge^2);
    % Minus the divergence of the unit "edge direction" field. back_diff is
    % the exact mirror of fwd_diff, so this is the true gradient of TV.
    g = -(back_diff(jx ./ mag, 1) + back_diff(jy ./ mag, 2) + back_diff(jz ./ mag, 3));
end


function d = fwd_diff(v, dim)
%FWD_DIFF v(i+1) - v(i) along dim, 0 at the last voxel (flat beyond the edge).
    sz = size(v);
    sz(dim) = 1;
    d = cat(dim, diff(v, 1, dim), zeros(sz));
end


function d = back_diff(v, dim)
%BACK_DIFF v(i) - v(i-1) along dim, with v(0) = 0. Mirror image of fwd_diff.
    sz = size(v);
    sz(dim) = 1;
    d = diff(cat(dim, zeros(sz), v), 1, dim);
end


function dose = pressure_to_dose(p, G, numPulses)
%PRESSURE_TO_DOSE Step 10. Crop the k-Wave padding, then
%  dose = p / (Gruneisen x density of CT1) x number of pulses
%  (the simulation worked with the dose of ONE pulse).
    p = p(G.offset(1) + (1:G.dims(1)), G.offset(2) + (1:G.dims(2)), G.offset(3) + (1:G.dims(3)));
    dose = p ./ G.gamma_rho * numPulses;
end


function v = to_kwave_grid(v, gp, kDims, fillValue)
%TO_KWAVE_GRID Move a dose-grid volume onto the k-Wave grid: add the water
%  border determine_sensor_mask_lateral asked for (gp; zero width when none),
%  then pad the far side up to the FFT-friendly size kDims. New voxels get
%  fillValue. gp.y_* pads array dim 1, gp.x_* dim 2, gp.z_* dim 3 (same
%  mapping as run_single_field_simulation).
    v = padarray(v, [gp.y_pre,  gp.x_pre,  gp.z_pre],  fillValue, 'pre');
    v = padarray(v, [gp.y_post, gp.x_post, gp.z_post], fillValue, 'post');
    v = padarray(v, kDims - size(v), fillValue, 'post');
end


function medium = build_medium_with_bath(sct, CONFIG)
%BUILD_MEDIUM_WITH_BATH create_acoustic_medium (the pipeline's only medium
%  builder, which handles in-body air) + water coupling bath outside the body
%  and on the couch. Same as the helper in run_standalone_field.
    medium  = create_acoustic_medium(sct, CONFIG);
    outside = ~logical(sct.bodyMask) | logical(sct.couchMask);
    medium.density(outside)     = CONFIG.uniform_density;
    medium.sound_speed(outside) = CONFIG.uniform_sound_speed;
    if numel(medium.alpha_coeff) > 1
        medium.alpha_coeff(outside) = CONFIG.uniform_alpha_coeff;
    end
    medium.gruneisen(outside)   = CONFIG.uniform_gruneisen;
end


function show_view(ax, vol, body, viewName, ijk, spacingMm, clim, cmap)
%SHOW_VIEW One orthogonal slice through voxel ijk = [row col slice] with the
%  body outline in cyan. Arrays are (Y, X, Z): transverse shows anterior at
%  the top; coronal and sagittal show the head at the top. NaN = blank.
    x = (0:size(vol, 2) - 1) * spacingMm(1);
    y = (0:size(vol, 1) - 1) * spacingMm(2);
    z = (0:size(vol, 3) - 1) * spacingMm(3);
    switch viewName
        case 'Transverse'
            img = vol(:, :, ijk(3));            outline = body(:, :, ijk(3));            h = x; v = y;
        case 'Coronal'
            img = squeeze(vol(ijk(1), :, :))';  outline = squeeze(body(ijk(1), :, :))';  h = x; v = z;
        case 'Sagittal'
            img = squeeze(vol(:, ijk(2), :))';  outline = squeeze(body(:, ijk(2), :))';  h = y; v = z;
    end
    im = imagesc(ax, h, v, img);
    set(im, 'AlphaData', ~isnan(img));
    axis(ax, 'image');
    if ~strcmp(viewName, 'Transverse')
        set(ax, 'YDir', 'normal');
    end
    caxis(ax, clim);
    colormap(ax, cmap);
    if any(outline(:)) && ~all(outline(:))
        hold(ax, 'on');
        contour(ax, h, v, double(outline), [0.5, 0.5], 'c', 'LineWidth', 0.8);
        hold(ax, 'off');
    end
    xlabel(ax, 'mm');
    ylabel(ax, 'mm');
end


function cmap = diverging_cmap()
%DIVERGING_CMAP Blue (low / negative) - white (middle) - red (high / positive).
    ramp = linspace(0, 1, 128)';
    cmap = [ramp, ramp, ones(128, 1); ones(128, 1), flipud(ramp), flipud(ramp)];
end
