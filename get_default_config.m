function CONFIG = get_default_config(varargin)
%GET_DEFAULT_CONFIG Default simulation CONFIG for the ETHOS k-Wave pipeline.
%
%   CONFIG = get_default_config()                    % all defaults
%   CONFIG = get_default_config(overrideStruct)      % struct of overrides
%   CONFIG = get_default_config('field', val, ...)   % name/value overrides
%
%   Single source of truth for the DEFAULT simulation parameters shared by the
%   standalone / sweep / study / test drivers. The defaults are the parameters
%   previously hard-coded at the top of run_standalone_simulation.m. Each driver
%   should call this and override only the fields it actually changes, e.g.
%
%       CONFIG = get_default_config();
%       CONFIG.reconstruction_method = 'das';
%       CONFIG.downscale_factor      = 2;
%
%   or equivalently  CONFIG = get_default_config('reconstruction_method','das', ...
%                                                 'downscale_factor', 2);
%
%   Overrides are applied over the defaults in the order given (name/value pairs)
%   or by field (struct). After overrides are applied, CONFIG.tissue_tables is
%   built from define_tissue_tables() and a '.uniform' row derived from the
%   (possibly overridden) CONFIG.uniform_* values -- unless the caller overrides
%   tissue_tables explicitly, in which case the supplied tables are kept.
%
%   This function covers SIMULATION config only. The pipeline orchestration
%   scripts (pipeline_setup / pipeline_simulate / step*) keep their own,
%   structurally different CONFIG blocks (patients, sessions, directories, MLC).
%
%   See also: define_tissue_tables, run_single_field_simulation,
%             run_standalone_field, create_acoustic_medium

    %% ======================== DEFAULTS ========================
    %  Ported from run_standalone_simulation.m (simulation parameters only).

    CONFIG = struct();

    % --- Machine / data selection (callers typically override these) ---
    CONFIG.working_dir    = '/mnt/weka/home/80030361/ETHOS_Simulations';
    CONFIG.patient_id     = '1194203';
    CONFIG.session        = 'Session_1';

    CONFIG.dose_filename  = 'dose_1194203_Session_1_reference_CT_1_B15_112.mat';
    % Standalone default geometry: CBCT1 (CT_1), regardless of which CBCT the
    % dose was computed on. Per-dose CBCT selection lives in the multi-field
    % driver, not here.
    CONFIG.cbct_filename  = 'CBCT1_resampled.mat';

    CONFIG.dose_file_override = '';
    CONFIG.cbct_file_override = '';

    % --- Sensor placement ---
    CONFIG.sensor_placement_method = 'determine_sensor_mask';
    CONFIG.sensor_x_index = 2;
    CONFIG.sensor_y_index = 4;

    % Physical 2D ultrasound array geometry (sparse element mask).
    % Kerf is derived inside determine_sensor_mask as (pitch - size).
    CONFIG.elements_per_side = 32;
    CONFIG.element_pitch_mm  = 4.35;   % Kerf ~= 0.7 mm
    CONFIG.element_size_mm   = 3.65;

    % --- Tissue / medium ---
    CONFIG.gruneisen_method = 'threshold_2';

    CONFIG.force_uniform_density     = false;
    CONFIG.force_uniform_sound_speed = false;
    CONFIG.force_uniform_attenuation = false;
    CONFIG.force_uniform_gruneisen   = false;

    CONFIG.uniform_density      = 1000;
    CONFIG.uniform_sound_speed  = 1540;
    CONFIG.uniform_alpha_coeff  = 0;
    CONFIG.uniform_alpha_power  = 1.1;
    CONFIG.uniform_gruneisen    = 1.0;

    % --- Dose / pulse / grid ---
    CONFIG.dose_per_pulse_cGy = 0.16;
    CONFIG.meterset           = 140;
    CONFIG.pml_size           = 10;
    CONFIG.cfl_number         = 0.3;
    CONFIG.Nt_scaling         = 6;      % >0: when air sets minC, divide Nt by this (0 = disabled)
    CONFIG.use_gpu            = true;

    % NOTE: run_standalone_simulation.m set correction_factor to 1.9 and then
    % immediately overrode it to 20 (the last assignment wins), so 20 is the
    % value it currently runs with; preserved here. The engine
    % (run_single_field_simulation) has its own internal default of 1.9.
    CONFIG.correction_factor  = 20;

    CONFIG.use_pressure_scale_correction = false;   % rescale by max(p0)/max(recon) before dose conversion
    CONFIG.mask_recon_to_dose_region     = true;    % zero recon dose outside the >1% dose mask

    % --- Reconstruction method: 'tr' (iterative time reversal) | 'das' | 'hybrid' ---
    CONFIG.reconstruction_method  = 'tr';
    CONFIG.num_time_reversal_iter = 1;
    CONFIG.convergence_tol        = 1e-3;

    % --- Pulse convolution / noise / deconvolution (set kernel 0 to disable) ---
    CONFIG.convolution_kernel  = 4e-6;   % Gaussian sigma (s)
    CONFIG.conv_noise_level    = 0.125;  % noise amplitude as fraction of peak sensor signal
    CONFIG.conv_deconv_lambda  = 1e-4;   % Wiener regularization for deconvolution

    CONFIG.downscale_factor = 1;
    CONFIG.use_grid_padding = true;

    % --- Output / display ---
    CONFIG.save_results = true;
    CONFIG.output_file  = 'standalone_recon_results.mat';
    CONFIG.plot_results = true;
    CONFIG.viz_smooth_sigma = 1.0;   % display-only dose smoothing (voxels); 0 disables

    % --- Normalization before comparison / gamma ---
    %   'max' : divide original and recon by their own max
    %   'lsq' : least-squares relative normalization (scale recon by
    %           g = sum(truth.*recon)/sum(recon.^2) over truth's >=10% region)
    CONFIG.normalize        = true;
    CONFIG.normalize_method = 'lsq';

    % --- Gamma logging ---
    CONFIG.log_gamma      = false;
    CONFIG.gamma_log_file = 'gamma_log.mat';

    % --- Diagnostics ---
    CONFIG.plot_exclusion_zone = false;

    %% ======================== APPLY OVERRIDES ========================

    provided = {};
    if nargin == 1 && isstruct(varargin{1})
        ov = varargin{1};
        provided = fieldnames(ov);
        for i = 1:numel(provided)
            CONFIG.(provided{i}) = ov.(provided{i});
        end
    elseif nargin >= 1
        if mod(nargin, 2) ~= 0
            error('get_default_config:BadArgs', ...
                'Name/value overrides must come in pairs (got %d arguments).', nargin);
        end
        for i = 1:2:nargin
            name = varargin{i};
            if ~(ischar(name) || (isstring(name) && isscalar(name)))
                error('get_default_config:BadArgs', ...
                    'Override names must be strings (argument %d).', i);
            end
            name = char(name);
            CONFIG.(name) = varargin{i+1};
            provided{end+1} = name; %#ok<AGROW>
        end
    end

    %% ======================== TISSUE TABLES ========================
    %  Built after overrides so the uniform row reflects any overridden
    %  uniform_* values. Skipped when the caller supplied tissue_tables.

    if ~ismember('tissue_tables', provided)
        CONFIG.tissue_tables = define_tissue_tables();
        CONFIG.tissue_tables.uniform = struct( ...
            'density',     CONFIG.uniform_density, ...
            'sound_speed', CONFIG.uniform_sound_speed, ...
            'alpha_coeff', CONFIG.uniform_alpha_coeff, ...
            'alpha_power', CONFIG.uniform_alpha_power, ...
            'gruneisen',   CONFIG.uniform_gruneisen);
    end
end
