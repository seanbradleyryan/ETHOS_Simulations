classdef SimResult
%SIMRESULT  One reconstructed dose plus its truth doses, geometry and provenance.
%
%   PURPOSE:
%   The single currency the analysis/plot layer consumes. A SimResult is
%   byte-identical whether it was produced by k-Wave (ethos.Simulator.run) or
%   loaded from disk (load_recon_dose_data), so ethos.StudyRunner and
%   ethos.SimPlotter never branch on where the dose came from -- "some studies
%   simulate, some just analyze" is one flag, not two code paths.
%
%   All dose volumes (recon_dose, rs_truth, ethos_truth) share the dose grid in
%   .metadata.dimensions, so they drop straight into gamma/SSIM with no further
%   resampling.
%
%   FIELDS:
%       recon_dose      - 3D reconstructed dose (Gy).
%       rs_truth        - Paired RayStation "truth" field/total dose (Gy) or [].
%       ethos_truth     - ETHOS RTDOSE truth resampled to the grid (Gy) or [].
%       density         - CT density volume (kg/m^3) for plot backgrounds or [].
%       sensor_mask     - 3D logical ultrasound-array mask (may be [] when a
%                         loaded result has no cached sensor; StudyRunner can
%                         fill it on demand via ensureSensor()).
%       sensor_info     - Companion diagnostics struct from determine_sensor_mask.
%       metadata        - Struct: dimensions, spacing, origin, beam_metadata.
%       rtplan          - Struct: beam_num, seg_num, plan_type, ct_label,
%                         gantry_angle, meterset, isocenter, jaw_x, jaw_y.
%       p0              - Initial pressure volume (sim-only; [] when loaded).
%       reconPressure   - Reconstructed pressure volume (sim-only; [] when loaded).
%       provenance      - Struct: source ('Simulate'|'Load'), config_hash,
%                         recon_file, forward_time_s, tr_time_s, num_iters_done,
%                         conv_max_pressure, conv_rel_change.
%       patient_id, session, gruneisen_method, config_hash - identity.
%
%   See also: ethos.Simulator, ethos.StudyRunner, ethos.SimPlotter,
%             load_recon_dose_data

    properties
        patient_id       char = ''
        session          char = ''
        gruneisen_method char = ''
        config_hash      char = ''

        recon_dose    = []
        rs_truth      = []
        ethos_truth   = []
        density       = []

        sensor_mask   = []
        sensor_info   = struct()

        metadata      = struct()
        rtplan        = struct()

        p0            = []
        reconPressure = []

        provenance    = struct('source', '', 'config_hash', '', ...
            'recon_file', '', 'forward_time_s', NaN, 'tr_time_s', NaN, ...
            'num_iters_done', NaN)
    end

    methods
        function sp = spacing(obj)
        %SPACING  [dx dy dz] mm from the metadata (row vector).
            if isfield(obj.metadata, 'spacing') && ~isempty(obj.metadata.spacing)
                sp = obj.metadata.spacing(:)';
            else
                sp = [1 1 1];
            end
        end

        function d = dims(obj)
        %DIMS  Grid [nx ny nz]; falls back to size(recon_dose).
            if isfield(obj.metadata, 'dimensions') && ~isempty(obj.metadata.dimensions)
                d = obj.metadata.dimensions(:)';
            else
                d = size(obj.recon_dose);
            end
        end

        function tf = isSimulated(obj)
        %ISSIMULATED  True when this result came from a fresh k-Wave run.
            tf = strcmpi(obj.provenance.source, 'Simulate');
        end

        function obj = withSensor(obj, sensor_mask, sensor_info)
        %WITHSENSOR  Return a copy with the sensor mask/info attached.
            obj.sensor_mask = sensor_mask;
            if nargin > 2, obj.sensor_info = sensor_info; end
        end

        function plot(obj, varargin)
        %PLOT  Convenience: full standalone-style panel set for this result.
            ethos.SimPlotter.all(obj, varargin{:});
        end
    end

    methods (Static)
        function r = fromField(loadOut, idx)
        %FROMFIELD  Build from one field of a load_recon_dose_data set/single out.
            if nargin < 2, idx = 1; end
            F = loadOut.fields(idx);
            r = ethos.SimResult();
            r.patient_id       = loadOut.patient_id;
            r.session          = loadOut.session;
            r.gruneisen_method = loadOut.gruneisen_method;
            r.config_hash      = loadOut.config_hash;
            r.metadata         = loadOut.metadata;
            r.recon_dose       = F.recon_dose;
            r.rs_truth         = F.rs_dose;
            r.rtplan           = F.rtplan;
            if isfield(loadOut, 'ethos_truth')
                r.ethos_truth = loadOut.ethos_truth;
            end
            r.density = ethos.SimResult.densityFromCBCT(F);
            r.provenance = struct('source', 'Load', ...
                'config_hash', loadOut.config_hash, ...
                'recon_file', F.recon_file, 'forward_time_s', NaN, ...
                'tr_time_s', NaN, 'num_iters_done', NaN);
        end

        function r = fromTotal(loadOut)
        %FROMTOTAL  Build from a load_recon_dose_data total-mode output.
            r = ethos.SimResult();
            r.patient_id       = loadOut.patient_id;
            r.session          = loadOut.session;
            r.gruneisen_method = loadOut.gruneisen_method;
            r.config_hash      = loadOut.config_hash;
            r.metadata         = loadOut.metadata;
            r.recon_dose       = loadOut.recon_dose;
            if isfield(loadOut, 'rs_dose'),    r.rs_truth    = loadOut.rs_dose;    end
            if isfield(loadOut, 'ethos_truth'), r.ethos_truth = loadOut.ethos_truth; end
            if isfield(loadOut, 'cbct')
                if isfield(loadOut.cbct, 'CT_1') && isfield(loadOut.cbct.CT_1, 'cubeDensity')
                    r.density = loadOut.cbct.CT_1.cubeDensity;
                end
            end
            r.provenance = struct('source', 'Load', ...
                'config_hash', loadOut.config_hash, 'recon_file', '', ...
                'forward_time_s', NaN, 'tr_time_s', NaN, 'num_iters_done', NaN);
        end

        function r = fromSim(field_dose, cbct_resampled, medium, recon_dose, ...
                sim_results, cfg)
        %FROMSIM  Build from a fresh run_single_field_simulation output.
            r = ethos.SimResult();
            r.patient_id       = cfg.patient_id;
            r.session          = cfg.session;
            r.gruneisen_method = cfg.gruneisen_method;
            r.config_hash      = cfg.hash();
            r.recon_dose       = recon_dose;
            if isfield(field_dose, 'dose_Gy'), r.rs_truth = field_dose.dose_Gy; end

            meta = struct();
            meta.dimensions = size(field_dose.dose_Gy);
            if isfield(cbct_resampled, 'spacing'), meta.spacing = cbct_resampled.spacing;
            elseif isfield(field_dose, 'spacing'), meta.spacing = field_dose.spacing; end
            if isfield(cbct_resampled, 'origin'), meta.origin = cbct_resampled.origin; end
            r.metadata = meta;

            if ~isempty(medium) && isfield(medium, 'density')
                r.density = medium.density;
            elseif isfield(cbct_resampled, 'cubeDensity')
                r.density = cbct_resampled.cubeDensity;
            end

            r.rtplan = ethos.SimResult.rtplanFromField(field_dose);

            % Diagnostics / sensor from sim_results.
            if isfield(sim_results, 'sensor_info'), r.sensor_info = sim_results.sensor_info; end
            prov = struct('source', 'Simulate', 'config_hash', r.config_hash, ...
                'recon_file', '');
            prov = ethos.SimResult.copyIf(prov, sim_results, 'forward_time_s', 'forward_time_s');
            prov = ethos.SimResult.copyIf(prov, sim_results, 'tr_time_s', 'tr_time_s');
            prov = ethos.SimResult.copyIf(prov, sim_results, 'num_iters_done', 'num_iters_done');
            prov = ethos.SimResult.copyIf(prov, sim_results, 'conv_max_pressure', 'conv_max_pressure');
            prov = ethos.SimResult.copyIf(prov, sim_results, 'conv_rel_change', 'conv_rel_change');
            r.provenance = prov;
        end
    end

    methods (Static, Access = private)
        function d = densityFromCBCT(F)
            d = [];
            if isfield(F, 'cbct') && ~isempty(F.cbct)
                if isfield(F.cbct, 'cubeDensity') && ~isempty(F.cbct.cubeDensity)
                    d = F.cbct.cubeDensity;
                elseif isfield(F.cbct, 'cubeHU') && ~isempty(F.cbct.cubeHU)
                    d = F.cbct.cubeHU;   % HU fallback for a grayscale background
                end
            end
        end

        function rt = rtplanFromField(fd)
            rt = struct('beam_num', NaN, 'seg_num', NaN, 'plan_type', '', ...
                'ct_label', '', 'gantry_angle', NaN, 'meterset', NaN, ...
                'isocenter', [NaN NaN NaN], 'jaw_x', [NaN NaN], 'jaw_y', [NaN NaN]);
            map = {'beam_index','beam_num'; 'beam_num','beam_num'; ...
                   'segment','seg_num'; 'seg_num','seg_num'; ...
                   'plan_type','plan_type'; 'ct_label','ct_label'; ...
                   'gantry_angle','gantry_angle'; 'meterset','meterset'; ...
                   'isocenter','isocenter'; 'jaw_x','jaw_x'; 'jaw_y','jaw_y'};
            for k = 1:size(map, 1)
                rt = ethos.SimResult.copyIf(rt, fd, map{k, 1}, map{k, 2});
            end
        end

        function s = copyIf(s, src, srcField, dstField)
            if isstruct(src) && isfield(src, srcField) && ~isempty(src.(srcField))
                s.(dstField) = src.(srcField);
            end
        end
    end
end
