function result = calibrate_noise_amp(patient_id, session, config)
%CALIBRATE_NOISE_AMP  Fixed electronic-noise amplitude (Pa) for a target mean SNR.
%
%   result = calibrate_noise_amp(patient_id, session, config)
%
%   PURPOSE:
%   Real electronic noise is a fixed amplitude set by the detection hardware, so
%   the SNR should FLOAT from beam to beam as the acoustic signal strength
%   changes -- weak (low-MU) segments seeing a lower SNR than strong ones. The
%   pipeline models this with one blanket noise amplitude (config.noise_amp_Pa)
%   shared by every field. This util picks that amplitude: it runs a few random
%   beam/segment forward simulations, measures each field's pre-noise sensor
%   signal peak, and returns the noise amplitude whose MEAN SNR across the
%   sampled fields equals config.target_snr (default 8):
%
%       noise_amp_Pa = mean(signal_peak) / target_snr
%
%   With SNR_i = signal_peak_i / noise_amp_Pa, the arithmetic mean of the
%   per-field SNRs is then exactly target_snr, while individual segments float
%   above and below it.
%
%   Only the forward simulation runs (config.forward_peak_only = true), so the
%   expensive time-reversal reconstruction is skipped -- calibration costs one
%   forward sim per sampled field, not a full recon.
%
%   INPUTS:
%       patient_id - Char/string patient identifier (e.g. '1194203').
%       session    - Char/string session name (e.g. 'Session_1').
%       config     - Simulation CONFIG (see get_default_config). Fields used:
%           .working_dir                 - path resolution (required)
%           .target_snr            [8]   - mean SNR to hit
%           .num_calibration_fields [5]  - number of random fields to sample
%           .calibration_seed      [ ]   - optional rng seed for a reproducible
%                                          random sample (unset = shuffle)
%           plus every engine knob (sensor, medium, pulse, ...), which must match
%           the settings the pipeline will run with so the peaks are comparable.
%
%   OUTPUT (struct):
%       .noise_amp_Pa       - Calibrated fixed noise amplitude (Pa). Assign this
%                             to CONFIG.noise_amp_Pa before running the pipeline.
%       .target_snr         - Target mean SNR used.
%       .mean_signal_peak   - Mean pre-noise sensor peak over the sampled fields.
%       .num_fields_sampled - Number of fields actually simulated.
%       .sampled_files      - Cell array of the sampled dose filenames.
%       .signal_peaks       - Per-sample pre-noise sensor peak (Pa).
%       .snr_per_field      - Per-sample SNR at noise_amp_Pa (signal_peak/amp).
%
%   EXAMPLE:
%       CONFIG = get_default_config();
%       CONFIG.working_dir = '/mnt/weka/home/80030361/ETHOS_Simulations';
%       cal = calibrate_noise_amp('1194203', 'Session_1', CONFIG);
%       CONFIG.noise_amp_Pa = cal.noise_amp_Pa;   % pipeline now uses fixed noise
%
%   NOTE: the sampled fields are simulated through run_standalone_field, which
%   wraps the outside-body/couch region as a uniform coupling bath. The pipeline
%   builds its medium straight from create_acoustic_medium (no explicit bath).
%   Outside-body voxels are water in both, so the peaks are directly comparable;
%   the calibration is a mean estimate over a random sample, not an exact match.
%
%   DEPENDENCIES:
%       list_processed_field_doses, load_field_dose_file, run_standalone_field.
%
%   See also: run_standalone_field, run_single_field_simulation,
%             get_default_config, resolve_noise_amp (in run_single_field_simulation)

    %% ======================== INPUT VALIDATION ========================

    if ~ischar(patient_id) && ~isstring(patient_id)
        error('calibrate_noise_amp:InvalidInput', ...
            'patient_id must be a string or character array.');
    end
    patient_id = char(patient_id);

    if ~ischar(session) && ~isstring(session)
        error('calibrate_noise_amp:InvalidInput', ...
            'session must be a string or character array.');
    end
    session = char(session);

    if ~isstruct(config) || ~isfield(config, 'working_dir')
        error('calibrate_noise_amp:MissingConfig', ...
            'config must be a struct containing working_dir.');
    end

    target_snr  = get_field(config, 'target_snr', 8);
    num_samples = get_field(config, 'num_calibration_fields', 5);
    if ~(isscalar(target_snr) && target_snr > 0)
        error('calibrate_noise_amp:BadTarget', 'target_snr must be a positive scalar.');
    end

    %% ======================== SAMPLE FIELDS ========================

    field_index = list_processed_field_doses(patient_id, session, config);
    n_fields    = numel(field_index);
    n_samp      = min(round(num_samples), n_fields);
    if n_samp < 1
        error('calibrate_noise_amp:NoFields', ...
            'No processed field doses found for %s / %s.', patient_id, session);
    end

    if isfield(config, 'calibration_seed') && ~isempty(config.calibration_seed)
        rng(config.calibration_seed);
    end
    sel = randperm(n_fields, n_samp);

    fprintf('[calibrate_noise_amp] Sampling %d of %d field(s) for a target mean SNR of %.2f.\n', ...
        n_samp, n_fields, target_snr);

    %% ======================== FORWARD-ONLY SIMS ========================
    %  Each sampled field runs its forward simulation only (forward_peak_only),
    %  returning the pre-noise sensor peak the noise amplitude is measured against.

    peaks = nan(n_samp, 1);
    for k = 1:n_samp
        fi = field_index(sel(k));
        fd = load_field_dose_file(fi.file);   % read ct_label to pick the CBCT geometry

        cfg = config;
        cfg.patient_id         = patient_id;
        cfg.session            = session;
        cfg.dose_file_override = fi.file;
        cfg.cbct_file_override = '';
        cfg.cbct_filename      = cbct_filename_for_label(fd);
        cfg.forward_peak_only  = true;    % engine returns after the peak; no recon
        cfg.blind_recon        = false;   % forward peak is on the field's own medium
        cfg.plot_results       = false;
        cfg.return_diagnostics = false;

        out = run_standalone_field(cfg);
        peaks(k) = read_peak(out);

        fprintf('  [%d/%d] %s : signal peak %.3e Pa\n', ...
            k, n_samp, fi.source_mat_filename, peaks(k));
    end

    %% ======================== CALIBRATE ========================

    valid = peaks(isfinite(peaks) & peaks > 0);
    if isempty(valid)
        error('calibrate_noise_amp:NoValidPeaks', ...
            ['None of the %d sampled fields produced a usable sensor peak ' ...
             '(check that the geometry/dose actually generate signal).'], n_samp);
    end
    mean_peak    = mean(valid);
    noise_amp_Pa = mean_peak / target_snr;
    snr_per_field = peaks / noise_amp_Pa;

    result = struct();
    result.noise_amp_Pa       = noise_amp_Pa;
    result.target_snr         = target_snr;
    result.mean_signal_peak   = mean_peak;
    result.num_fields_sampled = n_samp;
    result.sampled_files      = {field_index(sel).source_mat_filename}';
    result.signal_peaks       = peaks;
    result.snr_per_field      = snr_per_field;

    fprintf('[calibrate_noise_amp] Mean signal peak %.3e Pa over %d valid field(s).\n', ...
        mean_peak, numel(valid));
    fprintf('[calibrate_noise_amp] Calibrated noise amplitude: %.4e Pa (target mean SNR %.2f).\n', ...
        noise_amp_Pa, target_snr);
    fprintf('[calibrate_noise_amp] Per-field SNR at this amp: min %.2f, max %.2f.\n', ...
        min(snr_per_field), max(snr_per_field));
    fprintf('[calibrate_noise_amp] Set CONFIG.noise_amp_Pa = %.6e; then run the pipeline.\n', ...
        noise_amp_Pa);
end


%% ======================== LOCAL HELPERS ========================

function name = cbct_filename_for_label(fd)
%CBCT_FILENAME_FOR_LABEL  Resampled-CBCT filename for a field's ct_label.
%   CT_1 -> CBCT1; CT_2/CT_3 (synonyms) -> CBCT3. Defaults to CBCT1 when the
%   label is missing, matching the standalone default geometry.
    label = '';
    if isfield(fd, 'ct_label') && ~isempty(fd.ct_label)
        label = char(fd.ct_label);
    end
    if contains(label, '3') || contains(label, '2')
        name = 'CBCT3_resampled.mat';
    else
        name = 'CBCT1_resampled.mat';
    end
end


function peak = read_peak(out)
%READ_PEAK  Pull the pre-noise sensor peak out of a run_standalone_field result.
    peak = NaN;
    if isstruct(out) && isfield(out, 'sim_results') && isstruct(out.sim_results) ...
            && isfield(out.sim_results, 'sensor_signal_peak')
        peak = double(gather(out.sim_results.sensor_signal_peak));
    end
end


function v = get_field(s, f, d)
%GET_FIELD  Struct field with a default when absent or empty.
    if isfield(s, f) && ~isempty(s.(f))
        v = s.(f);
    else
        v = d;
    end
end
