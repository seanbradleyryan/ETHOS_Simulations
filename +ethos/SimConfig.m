classdef SimConfig
%SIMCONFIG  Value-object replacement for the copy-pasted CONFIG blocks.
%
%   PURPOSE:
%   One immutable, hashable configuration object that carries every parameter
%   the ETHOS photoacoustic pipeline reads. Replaces the ~100-line CONFIG
%   structs duplicated across run_*/study_*/test_* driver scripts. Static
%   factory presets (standalone/pipeline/convergenceSweep/gammaBatch) give the
%   common starting points; withOverrides() returns a modified copy (value
%   semantics keep it parfor-broadcast-safe); toStruct() emits a plain struct
%   so the existing physics functions (run_single_field_simulation,
%   create_acoustic_medium, determine_sensor_mask, load_recon_dose_data) keep
%   taking a bare `config` struct unchanged.
%
%   USAGE:
%       cfg = ethos.SimConfig.standalone();
%       cfg = cfg.withOverrides('num_time_reversal_iter', 20, 'use_gpu', false);
%       s   = cfg.toStruct();          % hand to any existing pipeline function
%       h   = cfg.hash();              % 8-char compute_sim_config_hash digest
%
%   NOTE ON THE CONFIG HASH:
%   hash() delegates to compute_sim_config_hash, which reads only an explicit
%   allow-list of sim-affecting fields. toStruct() emits the same field names a
%   pipeline CONFIG carries, so a SimConfig hashes identically to the CONFIG
%   that produced the on-disk recons. Analysis/plot/path knobs are ignored by
%   the hash by construction.
%
%   See also: compute_sim_config_hash, run_single_field_simulation,
%             ethos.Simulator, ethos.StudyRunner, ethos.DoseSource

    properties
        % ---- Identity / paths (not hashed) ----
        working_dir     char   = ''
        patient_id      char   = ''
        session         char   = ''
        treatment_site  char   = 'Pancreas'

        % ---- Grueneisen / tissue model ----
        gruneisen_method char  = 'threshold_2'
        tissue_tables          = []          % filled by toStruct() when empty

        % ---- LINAC / pulse ----
        dose_per_pulse_cGy     double = 0.16

        % ---- k-Wave numerics ----
        pml_size               double = 10
        cfl_number             double = 0.3
        Nt_scaling             double = 0
        use_gpu                logical = true

        % ---- Reconstruction ----
        reconstruction_method  char   = 'tr'    % 'tr' | 'das' | 'hybrid'
        num_time_reversal_iter double = 30
        convergence_tol        double = 1e-3
        correction_factor      double = 1.9
        use_pressure_scale_correction logical = false
        mask_recon_to_dose_region     logical = true
        use_attenuation        logical = true
        alpha_power            double = 1.1
        downscale_factor       double = 1
        use_grid_padding       logical = true
        normalize              logical = false

        % ---- Pulse convolution / noise / deconvolution ----
        convolution_kernel     double = 4e-6
        conv_noise_level       double = 0.01
        conv_deconv_lambda     double = 1e-4

        % ---- Force-uniform medium overrides ----
        force_uniform_density      logical = false
        force_uniform_sound_speed  logical = false
        force_uniform_attenuation  logical = false
        force_uniform_gruneisen    logical = false
        uniform_density        double = 1000
        uniform_sound_speed    double = 1540
        uniform_alpha_coeff    double = 0
        uniform_alpha_power    double = 1.1
        uniform_gruneisen      double = 1.0

        % ---- Sensor geometry ----
        sensor_placement_method char = 'full_plane_anterior'
        sensor_x_index         double = 1
        sensor_y_index         double = 1
        elements_per_side      double = 32
        element_pitch_mm       double = 4.35
        element_size_mm        double = 3.65

        % ---- Blind-recon (pipeline; hashed) ----
        blind_recon_ct3            logical = true
        blind_recon_reference_label char   = 'CT_1'

        % ---- Analysis (not hashed) ----
        dose_source            = ethos.DoseSource.Auto
        gamma_default          double = [3 3]   % [dose_pct, dist_mm]

        % ---- Visualization (not hashed) ----
        plot_results           logical = false
        viz_smooth_sigma       double = 1.0
        enable = struct('gammaMap', true, 'reconVsTruth', true, ...
            'sensorPlacement', true, 'doseDifference', true, ...
            'convergence', true, 'passRateCurve', true)
    end

    methods
        function obj = SimConfig(varargin)
        %SIMCONFIG  Construct with optional Name,Value overrides.
            if nargin > 0
                obj = obj.withOverrides(varargin{:});
            end
        end

        function obj = withOverrides(obj, varargin)
        %WITHOVERRIDES  Return a copy with the given properties changed.
        %   cfg2 = cfg.withOverrides('use_gpu', false, 'gamma_default', [2 2]);
            if mod(numel(varargin), 2) ~= 0
                error('ethos:SimConfig:BadOverrides', ...
                    'Overrides must be Name,Value pairs.');
            end
            for k = 1:2:numel(varargin)
                name = varargin{k};
                val  = varargin{k + 1};
                if ~ischar(name) && ~isstring(name)
                    error('ethos:SimConfig:BadOverrides', ...
                        'Override names must be strings.');
                end
                name = char(name);
                if ~isprop(obj, name)
                    error('ethos:SimConfig:UnknownField', ...
                        'SimConfig has no field ''%s''.', name);
                end
                if strcmp(name, 'dose_source')
                    val = ethos.DoseSource.fromChar(val);
                end
                obj.(name) = val;
            end
        end

        function obj = validate(obj)
        %VALIDATE  Check invariants; returns obj unchanged so it can chain.
            if ~ismember(lower(obj.reconstruction_method), {'tr', 'das', 'hybrid'})
                error('ethos:SimConfig:BadMethod', ...
                    'reconstruction_method must be tr | das | hybrid.');
            end
            if ~ismember(lower(obj.gruneisen_method), ...
                    {'uniform', 'threshold_1', 'threshold_2'})
                error('ethos:SimConfig:BadGruneisen', ...
                    'gruneisen_method must be uniform | threshold_1 | threshold_2.');
            end
            if numel(obj.gamma_default) ~= 2 || any(obj.gamma_default <= 0)
                error('ethos:SimConfig:BadGamma', ...
                    'gamma_default must be a 2-element positive [dose_pct dist_mm].');
            end
            if obj.downscale_factor <= 0
                error('ethos:SimConfig:BadDownscale', ...
                    'downscale_factor must be > 0.');
            end
        end

        function [h, canonical] = hash(obj)
        %HASH  8-char compute_sim_config_hash digest of the sim-affecting fields.
            [h, canonical] = compute_sim_config_hash(obj.toStruct());
        end

        function obj = withTissue(obj, method, tissue, field, value)
        %WITHTISSUE  Return a copy with one tissue-table entry overridden.
        %   Materializes the built-in tables when empty, then sets
        %   tables.(method).(field) for the named tissue. Supports the air
        %   sub-row (tissue='air' -> tables.(method).air.(field)) and any named
        %   tissue row (indexed by tables.(method).tissue_names). Used to sweep
        %   e.g. the in-body air sound speed without touching create_acoustic_medium.
            tables = obj.tissue_tables;
            if isempty(tables), tables = ethos.SimConfig.tissueTables(); end
            if ~isfield(tables, method)
                error('ethos:SimConfig:BadMethod', ...
                    'tissue_tables has no method ''%s''.', method);
            end
            if strcmpi(tissue, 'air')
                tables.(method).air.(field) = value;
            else
                idx = find(strcmpi(tables.(method).tissue_names, tissue), 1);
                if isempty(idx)
                    error('ethos:SimConfig:BadTissue', ...
                        'No tissue ''%s'' in method ''%s''.', tissue, method);
                end
                tables.(method).(field)(idx) = value;
            end
            obj.tissue_tables = tables;
        end

        function s = toStruct(obj)
        %TOSTRUCT  Plain struct for the existing pipeline functions.
        %   Every public property becomes a struct field. dose_source is
        %   flattened to char; tissue_tables is populated from the built-in
        %   tables when left empty so create_acoustic_medium has its lookup.
            s = struct();
            props = properties(obj);
            for i = 1:numel(props)
                p = props{i};
                v = obj.(p);
                if isa(v, 'ethos.DoseSource')
                    v = char(v);
                end
                s.(p) = v;
            end
            if isempty(s.tissue_tables)
                s.tissue_tables = ethos.SimConfig.tissueTables();
            end
        end
    end

    methods (Static)
        function cfg = standalone()
        %STANDALONE  Single-field standalone-simulation preset.
            cfg = ethos.SimConfig();
            cfg.sensor_placement_method = 'determine_sensor_mask';
            cfg.sensor_x_index          = 2;
            cfg.sensor_y_index          = 4;
            cfg.num_time_reversal_iter  = 10;
            cfg.conv_noise_level        = 0.125;
            cfg.normalize               = true;
            cfg.plot_results            = true;
            cfg.dose_source             = ethos.DoseSource.Auto;
        end

        function cfg = pipeline()
        %PIPELINE  Full-plan simulation preset (matches pipeline_simulate).
            cfg = ethos.SimConfig();
            cfg.sensor_placement_method = 'determine_sensor_mask';
            cfg.num_time_reversal_iter  = 30;
            cfg.blind_recon_ct3         = true;
            cfg.plot_results            = false;
            cfg.dose_source             = ethos.DoseSource.Auto;
        end

        function cfg = convergenceSweep()
        %CONVERGENCESWEEP  Preset for TR-iteration / numeric sweeps (no plots).
            cfg = ethos.SimConfig.standalone();
            cfg.plot_results = false;
            cfg.dose_source  = ethos.DoseSource.Simulate;
        end

        function cfg = gammaBatch()
        %GAMMABATCH  Per-dose gamma/SSIM batch preset (load-or-sim, 3%/3mm).
            cfg = ethos.SimConfig.standalone();
            cfg.plot_results  = false;
            cfg.gamma_default = [3 3];
            cfg.dose_source   = ethos.DoseSource.Auto;
        end

        function tables = tissueTables()
        %TISSUETABLES  Tissue-property lookup tables for HU thresholding.
        %   Single source of truth for the class layer, mirroring the canonical
        %   define_tissue_tables() inside pipeline_simulate (including the
        %   threshold_2 in-body air row: gruneisen 0 -> no PA signal).
            tables.threshold_1 = struct();
            tables.threshold_1.hu_boundaries = [-1000, -900, -500, -200, -50, 13, 50, 80, 300, 3000, Inf];
            tables.threshold_1.tissue_names  = {'Air','Lung','Fat','Water','Blood','Muscle','SoftTissue','Bone','Metal'};
            tables.threshold_1.density       = [1.2,  400,  920,  1000, 1060, 1050, 1040, 1900, 7800];
            tables.threshold_1.sound_speed   = [343,  600,  1450, 1480, 1575, 1580, 1540, 3200, 5900];
            tables.threshold_1.alpha_coeff   = [0,    0.5,  0.48, 0.0022, 0.2, 0.5, 0.5,  10,   0];
            tables.threshold_1.alpha_power   = [1.0,  1.5,  1.5,  2.0,  1.3,  1.0, 1.1,  1.0,  1.0];
            tables.threshold_1.gruneisen     = [0,    0.5,  0.7,  0.11, 0.15, 0.2, 1.0,  0,    0];

            tables.threshold_2 = struct();
            tables.threshold_2.hu_boundaries = [-1000, -200, -50, 100, Inf];
            tables.threshold_2.tissue_names  = {'Water','Fat','SoftTissue','Bone'};
            tables.threshold_2.density       = [1000,   920,  1040,        1900];
            tables.threshold_2.sound_speed   = [1480,   1450, 1540,        3200];
            tables.threshold_2.alpha_coeff   = [0.0022, 0.48, 0.5,         10];
            tables.threshold_2.alpha_power   = [2.0,    1.5,  1.1,         1.0];
            tables.threshold_2.gruneisen     = [0.11,   0.7,  1.0,         1.0];
            tables.threshold_2.air_hu_threshold = -300;
            tables.threshold_2.air = struct( ...
                'density',     1.2, ...
                'sound_speed', 343, ...
                'alpha_coeff', 0, ...
                'alpha_power', 1.0, ...
                'gruneisen',   0);

            tables.uniform = struct();
            tables.uniform.density      = 1000;
            tables.uniform.sound_speed  = 1540;
            tables.uniform.alpha_coeff  = 0;
            tables.uniform.alpha_power  = 1.1;
            tables.uniform.gruneisen    = 1.0;
        end
    end
end
