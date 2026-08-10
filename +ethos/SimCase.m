classdef SimCase
%SIMCASE  One immutable simulation-input bundle (dose + geometry + medium + sensor).
%
%   PURPOSE:
%   Packages everything run_single_field_simulation needs for one field into a
%   value object that is safe to broadcast into a parfor. Static factories build
%   cases from the processed pipeline outputs:
%       fromDoseFile - one field selected by Beam/Segment (materialized).
%       fromProcessed- every field of a session (wraps load_processed_data).
%       list         - a filtered, LIGHTWEIGHT array of cases (descriptors only;
%                      each worker materializes its own to keep host RAM down).
%
%   A lightweight case carries only identity + the dose file path; materialize()
%   loads the field dose, builds the acoustic medium (create_acoustic_medium),
%   and precomputes the sensor mask (determine_sensor_mask) on demand.
%
%   See also: load_processed_data, list_processed_field_doses,
%             load_field_dose_file, create_acoustic_medium, determine_sensor_mask,
%             ethos.Simulator, ethos.SimConfig

    properties
        patient_id  char = ''
        session     char = ''

        % Descriptor (always present)
        dose_file            char = ''
        source_mat_filename  char = ''
        beam    = []
        segment = []
        plan_type char = ''
        ct_label  char = ''

        % Bundle (populated by materialize)
        field_dose         = []
        cbct_resampled     = []
        medium             = []
        beam_metadata      = []
        precomputed_sensor = []
        materialized logical = false
    end

    methods
        function s = label(obj)
        %LABEL  Short human string for logs / plot titles.
            b = num2str_or(obj.beam);
            g = num2str_or(obj.segment);
            extra = strtrim(sprintf('%s %s', obj.plan_type, obj.ct_label));
            s = strtrim(sprintf('B%s_%s %s', b, g, extra));
        end

        function obj = materialize(obj, cfg)
        %MATERIALIZE  Load the dose, build the medium + sensor. No-op if already done.
            if obj.materialized
                return;
            end
            cfgS = cfg.toStruct();

            if isempty(obj.field_dose)
                if isempty(obj.dose_file)
                    error('ethos:SimCase:NoDose', ...
                        'SimCase has neither a materialized field_dose nor a dose_file.');
                end
                obj.field_dose = load_field_dose_file(obj.dose_file);
            end

            processed_dir = ethos.SimCase.processedDir(cfg, obj.patient_id, obj.session);

            % ---- Geometry: per-label CBCT (preferred) or session SCT ----
            if isempty(obj.cbct_resampled)
                obj.cbct_resampled = ethos.SimCase.resolveGeometry( ...
                    processed_dir, obj.ct_label, obj.field_dose);
            end

            % ---- Beam metadata (whole plan) for the exclusion zone ----
            if isempty(obj.beam_metadata)
                obj.beam_metadata = ethos.SimCase.loadBeamMetadata(processed_dir);
            end

            % ---- Deterministic placement dose (summed plan) if present ----
            if ~isfield(obj.field_dose, 'total_dose_Gy') || isempty(obj.field_dose.total_dose_Gy)
                total = ethos.SimCase.loadTotalRSDose(processed_dir, size(obj.field_dose.dose_Gy));
                if ~isempty(total)
                    obj.field_dose.total_dose_Gy = total;
                end
            end

            % ---- Acoustic medium ----
            obj.medium = create_acoustic_medium(obj.cbct_resampled, cfgS);

            % ---- Sensor (only the determine_sensor_mask method is precomputed) ----
            if strcmp(cfgS.sensor_placement_method, 'determine_sensor_mask')
                obj.precomputed_sensor = ethos.SimCase.buildSensor( ...
                    obj.cbct_resampled, obj.field_dose, obj.beam_metadata, cfgS);
            end

            obj.materialized = true;
        end
    end

    methods (Static)
        function obj = fromDoseFile(patient_id, session, cfg, varargin)
        %FROMDOSEFILE  One materialized case selected by Beam/Segment (+filters).
            cases = ethos.SimCase.list(patient_id, session, cfg, varargin{:});
            if isempty(cases)
                error('ethos:SimCase:NoMatch', ...
                    'No processed field matches the given filters.');
            end
            if numel(cases) > 1
                segs = zeros(1, numel(cases));
                for k = 1:numel(cases)
                    if ~isempty(cases(k).segment), segs(k) = cases(k).segment; end
                end
                error('ethos:SimCase:Ambiguous', ...
                    '%d fields match; narrow with Beam/Segment (segments: %s).', ...
                    numel(cases), mat2str(segs));
            end
            obj = cases.materialize(cfg);
        end

        function cases = list(patient_id, session, cfg, varargin)
        %LIST  Lightweight (unmaterialized) case array, filtered like load_recon_dose_data.
        %   Name/Value: 'Beam', 'Segment', 'PlanType' (any|adapted|reference),
        %   'CTLabel' (any|CT_1|CT_3).
            opts = ethos.SimCase.parseFilters(varargin);
            [field_index, metadata] = list_processed_field_doses(patient_id, session, cfg.toStruct());
            bm = [];
            if isfield(metadata, 'beam_metadata'), bm = metadata.beam_metadata; end

            keep = ethos.SimCase();  % template for preallocation type
            cases = repmat(keep, 1, 0);
            for i = 1:numel(field_index)
                [pt, ct] = ethos.SimCase.parseTokens(field_index(i).source_mat_filename);
                if ~isempty(opts.Beam)    && field_index(i).beam_index ~= opts.Beam,   continue; end
                if ~isempty(opts.Segment) && field_index(i).segment    ~= opts.Segment, continue; end
                if ~strcmp(opts.PlanType, 'any') && ~strcmpi(pt, opts.PlanType),        continue; end
                if ~strcmp(opts.CTLabel,  'any') && ~strcmpi(ct, opts.CTLabel),         continue; end

                c = ethos.SimCase();
                c.patient_id          = char(patient_id);
                c.session             = char(session);
                c.dose_file           = field_index(i).file;
                c.source_mat_filename = field_index(i).source_mat_filename;
                c.beam                = field_index(i).beam_index;
                c.segment             = field_index(i).segment;
                c.plan_type           = pt;
                c.ct_label            = ct;
                c.beam_metadata       = bm;   % small; shared descriptor
                cases(end + 1) = c; %#ok<AGROW>
            end
        end

        function cases = fromProcessed(patient_id, session, cfg)
        %FROMPROCESSED  Materialized case per field of a session (folds in load_processed_data).
            cfgS = cfg.toStruct();
            [field_doses, sct_resampled, total_rs_dose, metadata] = ...
                load_processed_data(patient_id, session, cfgS);

            processed_dir = ethos.SimCase.processedDir(cfg, patient_id, session);
            cbct_cache = ethos.SimCase.loadCBCTCache(processed_dir);
            bm = [];
            if isfield(metadata, 'beam_metadata'), bm = metadata.beam_metadata; end

            valid = ~cellfun(@isempty, field_doses);
            fds   = field_doses(valid);

            template = ethos.SimCase();
            cases = repmat(template, 1, numel(fds));
            for i = 1:numel(fds)
                fd = fds{i};
                if ~isempty(total_rs_dose) && ...
                        (~isfield(fd, 'total_dose_Gy') || isempty(fd.total_dose_Gy))
                    fd.total_dose_Gy = total_rs_dose;
                end
                [pt, ct] = ethos.SimCase.fieldTokens(fd);

                c = ethos.SimCase();
                c.patient_id          = char(patient_id);
                c.session             = char(session);
                c.source_mat_filename = ethos.SimCase.getName(fd);
                c.beam                = ethos.SimCase.getFieldNum(fd);
                c.segment             = ethos.SimCase.getSegNum(fd);
                c.plan_type           = pt;
                c.ct_label            = ct;
                c.field_dose          = fd;
                c.beam_metadata       = bm;

                % Prefer the per-label CBCT; fall back to session SCT.
                geom = ethos.SimCase.pickGeometry(cbct_cache, ct, sct_resampled);
                c.cbct_resampled = geom;
                c.medium = create_acoustic_medium(geom, cfgS);
                if strcmp(cfgS.sensor_placement_method, 'determine_sensor_mask')
                    c.precomputed_sensor = ethos.SimCase.buildSensor(geom, fd, bm, cfgS);
                end
                c.materialized = true;
                cases(i) = c;
            end
        end

        function precomp = buildSensor(geom, field_dose, beam_metadata, cfgS)
        %BUILDSENSOR  Wrap determine_sensor_mask; mirrors compute_plan_sensor_mask.
        %   Reused by the sim path (SimCase.materialize/fromProcessed) and by the
        %   load path (ethos.StudyRunner.ensureSensor) so both build the array
        %   identically. Returns [] when placement is not applicable.
            precomp = [];
            if isempty(geom) || ~isfield(geom, 'bodyMask') || isempty(geom.bodyMask)
                return;
            end

            sct = struct();
            if isfield(geom, 'cubeHU'),   sct.cubeHU   = geom.cubeHU;   end
            sct.bodyMask = geom.bodyMask;
            if isfield(geom, 'couchMask') && ~isempty(geom.couchMask)
                sct.couchMask = geom.couchMask;
            else
                sct.couchMask = false(size(geom.bodyMask));
            end
            if isfield(geom, 'origin') && ~isempty(geom.origin)
                sct.origin = geom.origin;
            else
                sct.origin = [0, 0, 0];
            end
            if isfield(geom, 'spacing'), sct.spacing = geom.spacing;
            elseif isfield(field_dose, 'spacing'), sct.spacing = field_dose.spacing; end

            fd = field_dose;
            if ~isfield(fd, 'origin')  || isempty(fd.origin),  fd.origin  = sct.origin;  end
            if ~isfield(fd, 'spacing') || isempty(fd.spacing), fd.spacing = sct.spacing; end
            if ~isfield(fd, 'dimensions') || isempty(fd.dimensions)
                fd.dimensions = size(fd.dose_Gy);
            end

            [sensor_mask, sensor_info] = determine_sensor_mask(sct, fd, beam_metadata, cfgS);
            precomp = struct('sensor_mask', sensor_mask, 'sensor_info', sensor_info, ...
                'base_dimensions', size(geom.bodyMask));
        end
    end

    methods (Static, Access = private)
        function opts = parseFilters(args)
            opts = struct('Beam', [], 'Segment', [], 'PlanType', 'any', 'CTLabel', 'any');
            if mod(numel(args), 2) ~= 0
                error('ethos:SimCase:BadOptions', 'Filters must be Name,Value pairs.');
            end
            for k = 1:2:numel(args)
                name = lower(char(args{k}));
                val  = args{k + 1};
                switch name
                    case 'beam',     opts.Beam    = val;
                    case 'segment',  opts.Segment = val;
                    case 'plantype', opts.PlanType = lower(char(val));
                    case 'ctlabel'
                        v = char(val);
                        if strcmpi(v, 'any')
                            opts.CTLabel = 'any';
                        else
                            tok = regexp(strrep(v, '-', '_'), '(\d+)', 'tokens', 'once');
                            if ~isempty(tok), opts.CTLabel = sprintf('CT_%s', tok{1}); else, opts.CTLabel = v; end
                        end
                    otherwise
                        error('ethos:SimCase:UnknownOption', 'Unknown filter ''%s''.', name);
                end
            end
        end

        function [pt, ct] = parseTokens(name)
            pt = ''; ct = '';
            tok = regexp(name, '_(adapted|reference)(?:_CT_(\d+))?_B\d+_\d+\.mat$', ...
                'tokens', 'once');
            if ~isempty(tok)
                pt = tok{1};
                if numel(tok) >= 2 && ~isempty(tok{2}), ct = sprintf('CT_%s', tok{2}); end
            end
        end

        function [pt, ct] = fieldTokens(fd)
            pt = ''; ct = '';
            if isfield(fd, 'plan_type') && ~isempty(fd.plan_type), pt = char(fd.plan_type); end
            if isfield(fd, 'ct_label') && ~isempty(fd.ct_label)
                ct = strrep(char(fd.ct_label), '-', '_');
            end
            if (isempty(pt) || isempty(ct)) && isfield(fd, 'source_mat_filename')
                [pt2, ct2] = ethos.SimCase.parseTokens(fd.source_mat_filename);
                if isempty(pt), pt = pt2; end
                if isempty(ct), ct = ct2; end
            end
        end

        function dir_ = processedDir(cfg, patient_id, session)
            dir_ = fullfile(cfg.working_dir, 'RayStationFiles', ...
                char(patient_id), char(session), 'processed');
        end

        function geom = resolveGeometry(processed_dir, ct_label, field_dose)
            cache = ethos.SimCase.loadCBCTCache(processed_dir);
            geom  = ethos.SimCase.pickGeometry(cache, ct_label, []);
            if isempty(geom)
                % Session SCT fallback.
                f = fullfile(processed_dir, 'sct_resampled.mat');
                if isfile(f)
                    L = load(f, 'sct_resampled');
                    geom = L.sct_resampled;
                end
            end
            if isempty(geom)
                error('ethos:SimCase:NoGeometry', ...
                    'No CBCT/SCT geometry found in %s for ct_label ''%s''.', ...
                    processed_dir, ct_label); %#ok<*CTPCT>
            end
            if ~isfield(geom, 'spacing') && isfield(field_dose, 'spacing')
                geom.spacing = field_dose.spacing;
            end
        end

        function geom = pickGeometry(cache, ct_label, sct_fallback)
            geom = [];
            lbl = ct_label;
            if isempty(lbl), lbl = 'CT_1'; end
            if isstruct(cache) && isfield(cache, lbl)
                geom = cache.(lbl);
            elseif isstruct(cache)
                fn = fieldnames(cache);
                if ~isempty(fn), geom = cache.(fn{1}); end
            end
            if isempty(geom) && ~isempty(sct_fallback)
                geom = sct_fallback;
            end
        end

        function cache = loadCBCTCache(processed_dir)
            specs = { 'CBCT1_resampled.mat', 'CBCT1_resampled'; ...
                      'CBCT3_resampled.mat', 'CBCT3_resampled' };
            cache = struct();
            for k = 1:size(specs, 1)
                fpath = fullfile(processed_dir, specs{k, 1});
                if ~isfile(fpath), continue; end
                L = load(fpath, specs{k, 2});
                c = L.(specs{k, 2});
                if isfield(c, 'ct_label') && ~isempty(c.ct_label)
                    lbl = strrep(char(c.ct_label), '-', '_');
                else
                    tok = regexp(specs{k, 1}, 'CBCT(\d+)', 'tokens', 'once');
                    lbl = sprintf('CT_%s', tok{1});
                end
                cache.(lbl) = c;
            end
        end

        function bm = loadBeamMetadata(processed_dir)
            bm = [];
            f = fullfile(processed_dir, 'metadata.mat');
            if isfile(f)
                L = load(f, 'metadata');
                if isfield(L, 'metadata') && isfield(L.metadata, 'beam_metadata')
                    bm = L.metadata.beam_metadata;
                end
            end
        end

        function d = loadTotalRSDose(processed_dir, target_size)
            d = [];
            f = fullfile(processed_dir, 'total_rs_dose.mat');
            if ~isfile(f), return; end
            L = load(f);
            if isfield(L, 'total_rs_dose_sparse')
                d = reshape(full(L.total_rs_dose_sparse), L.total_rs_dose_dims);
            elseif isfield(L, 'total_rs_dose')
                d = L.total_rs_dose;
            end
            if ~isempty(d) && nargin > 1 && ~isequal(size(d), target_size)
                d = [];   % size-mismatched; let placement fall back to field dose
            end
        end

        function n = getName(fd)
            n = '';
            if isfield(fd, 'source_mat_filename'), n = char(fd.source_mat_filename); end
        end

        function v = getFieldNum(fd)
            v = [];
            if isfield(fd, 'beam_index') && ~isempty(fd.beam_index), v = fd.beam_index;
            elseif isfield(fd, 'beam_num') && ~isempty(fd.beam_num), v = fd.beam_num; end
        end

        function v = getSegNum(fd)
            v = [];
            if isfield(fd, 'segment') && ~isempty(fd.segment), v = fd.segment;
            elseif isfield(fd, 'seg_num') && ~isempty(fd.seg_num), v = fd.seg_num; end
        end
    end
end


function s = num2str_or(v)
%NUM2STR_OR  num2str for a scalar, '?' for empty. (File-local helper.)
    if isempty(v)
        s = '?';
    else
        s = num2str(v);
    end
end
