function out = load_recon_dose_data(patient_id, session, config, varargin)
%LOAD_RECON_DOSE_DATA Load already-computed reconstructed doses for analysis.
%
%   out = load_recon_dose_data(patient_id, session, config, Name, Value, ...)
%
%   PURPOSE:
%   Load reconstructed (recon) doses produced by the main simulation pipeline
%   WITHOUT re-running k-Wave, so the standalone-simulation / comparison /
%   study analyses can be repeated cheaply. For one field, a filtered set, or
%   the summed total, it returns the reconstructed dose plus everything those
%   analyses need: the paired RayStation "truth" field dose, the ETHOS RTDOSE
%   truth (resampled onto the dose grid), the associated CBCT geometry, and the
%   RTPLAN statistics for each dose.
%
%   The recon files are written by pipeline_simulate and keyed by an 8-char
%   CONFIG hash (compute_sim_config_hash). The right hash is auto-discovered
%   from the simulation folder when possible, otherwise derived from a passed
%   CONFIG, otherwise overridable via the 'Hash' option.
%
%   INPUTS:
%       patient_id - Char/string, patient identifier (e.g., '1194203').
%       session    - Char/string, session name (e.g., 'Session_1').
%       config     - Struct with:
%           .working_dir      - REQUIRED. Base directory path.
%           .treatment_site   - Optional (default 'Pancreas'). Used to locate
%                               the ETHOS RTDOSE truth.
%           .gruneisen_method - Optional. Simulation-results subfolder. When
%                               omitted, the single existing method folder is
%                               used (errors if several exist).
%           May also carry the full simulation CONFIG so the hash can be
%           computed exactly via compute_sim_config_hash.
%
%   NAME/VALUE OPTIONS:
%       'Mode'         - 'total' (default) | 'set' | 'single'
%       'Beam'         - Scalar beam number. Required for 'single'; filters
%                        'set'.
%       'Segment'      - Scalar segment number. Optional filter; for 'single'
%                        it narrows to one field when a beam has several
%                        segments.
%       'PlanType'     - 'any' (default) | 'adapted' | 'reference'
%       'CTLabel'      - 'any' (default) | 'CT_1' | 'CT_3'
%       'Hash'         - Explicit 8-char config hash override (e.g. 'a1b2c3d4').
%       'IncludeRS'    - Logical, default true. Load the paired RayStation dose.
%       'IncludeEthos' - Logical, default true. Load the ETHOS truth dose.
%       'IncludeCBCT'  - Logical, default true. Load the CBCT geometry.
%
%   OUTPUTS:
%       out - Struct with common fields:
%           .patient_id, .session, .gruneisen_method, .config_hash
%           .metadata     - Grid metadata (dimensions, spacing, origin,
%                           beam_metadata).
%           .ethos_truth  - 3D ETHOS truth dose (Gy) resampled to the dose grid
%                           (only when IncludeEthos).
%         Mode 'total':
%           .recon_dose   - Total reconstructed dose (3D Gy).
%           .rs_dose      - Total RayStation dose (3D Gy) (when IncludeRS).
%           .cbct         - Struct with .CT_1 / .CT_3 geometry (when IncludeCBCT).
%         Mode 'single' / 'set':
%           .fields       - Struct array (1 element per matched field; 'single'
%                           yields a 1-element array), each with:
%               .recon_dose          - 3D reconstructed dose (Gy)
%               .rs_dose             - 3D RayStation dose (Gy) (when IncludeRS)
%               .cbct                - CBCT struct for this field (when IncludeCBCT)
%               .rtplan              - Struct: beam_num, seg_num, plan_type,
%                                      ct_label, beam_name, gantry_angle,
%                                      meterset, isocenter, jaw_x, jaw_y
%               .source_mat_filename - Input dose filename
%               .recon_file          - Full path to the loaded recon .mat
%
%   ALGORITHM:
%       1. Validate inputs; parse options.
%       2. Resolve gruneisen_method and the config hash (auto-discover -> config
%          -> 'Hash' override).
%       3. total: load total_recon_dose_<hash>.mat + total_rs_dose + both CBCTs.
%          set/single: index dose files, filter, load <base>_recon_<hash>.mat
%          per field plus its RS dose, RTPLAN stats and CBCT.
%       4. When IncludeEthos, load RTDOSE truth and resample to the dose grid.
%
%   EXAMPLE:
%       config.working_dir      = '/mnt/weka/home/80030361/ETHOS_Simulations';
%       config.treatment_site   = 'Pancreas';
%       config.gruneisen_method = 'threshold_2';
%       T = load_recon_dose_data('1194203','Session_1',config);
%       S = load_recon_dose_data('1194203','Session_1',config, ...
%               'Mode','single','Beam',15,'Segment',112);
%       A = load_recon_dose_data('1194203','Session_1',config, ...
%               'Mode','set','PlanType','adapted','CTLabel','CT_3');
%
%   DEPENDENCIES:
%       - compute_sim_config_hash, list_processed_field_doses,
%         load_field_dose_file
%       - Pipeline outputs from pipeline_simulate and step15_process_doses.
%
%   See also: compute_sim_config_hash, list_processed_field_doses,
%             load_field_dose_file, load_processed_data, step3_analysis

    %% ======================== INPUT VALIDATION ========================

    if ~ischar(patient_id) && ~isstring(patient_id)
        error('load_recon_dose_data:InvalidInput', ...
            'patient_id must be a string or character array. Received: %s', ...
            class(patient_id));
    end
    patient_id = char(patient_id);

    if ~ischar(session) && ~isstring(session)
        error('load_recon_dose_data:InvalidInput', ...
            'session must be a string or character array. Received: %s', ...
            class(session));
    end
    session = char(session);

    if ~isstruct(config) || ~isfield(config, 'working_dir')
        error('load_recon_dose_data:MissingConfig', ...
            'config must be a struct containing working_dir.');
    end

    opts = parse_options(varargin);

    %% ======================== RESOLVE PATHS + HASH ========================

    processed_dir = fullfile(config.working_dir, 'RayStationFiles', ...
        patient_id, session, 'processed');

    gruneisen_method = resolve_gruneisen_method(patient_id, session, config);
    sim_dir = fullfile(config.working_dir, 'SimulationResults', ...
        patient_id, session, gruneisen_method);
    if ~isfolder(sim_dir)
        error('load_recon_dose_data:NoSimDir', ...
            'Simulation results directory not found: %s', sim_dir);
    end

    hash8 = resolve_config_hash(sim_dir, opts.Mode, config, opts.Hash);

    %% ======================== ASSEMBLE OUTPUT ========================

    out                  = struct();
    out.patient_id       = patient_id;
    out.session          = session;
    out.gruneisen_method = gruneisen_method;
    out.config_hash      = hash8;
    out.metadata         = struct();

    fprintf(['[load_recon_dose_data] %s / %s | method=%s | hash=%s | ' ...
             'mode=%s\n'], patient_id, session, gruneisen_method, hash8, ...
             opts.Mode);

    switch opts.Mode
        case 'total'
            out = load_total_branch(out, sim_dir, processed_dir, hash8, opts);
        case {'set', 'single'}
            out = load_set_branch(out, patient_id, session, config, ...
                sim_dir, processed_dir, hash8, opts);
    end

    if opts.IncludeEthos
        out.ethos_truth = load_ethos_truth_resampled( ...
            patient_id, session, config, out.metadata);
    end

    fprintf('[load_recon_dose_data] Done.\n');
end


%% =========================================================================
%  OPTION PARSING
%% =========================================================================

function opts = parse_options(args)
    opts = struct('Mode', 'total', 'Beam', [], 'Segment', [], ...
        'PlanType', 'any', 'CTLabel', 'any', 'Hash', '', ...
        'IncludeRS', true, 'IncludeEthos', true, 'IncludeCBCT', true);

    if mod(numel(args), 2) ~= 0
        error('load_recon_dose_data:BadOptions', ...
            'Name/value options must be provided in pairs.');
    end

    for k = 1:2:numel(args)
        name = args{k};
        val  = args{k + 1};
        if ~ischar(name) && ~isstring(name)
            error('load_recon_dose_data:BadOptions', ...
                'Option names must be strings.');
        end
        switch lower(char(name))
            case 'mode',         opts.Mode    = validate_mode(val);
            case 'beam',         opts.Beam    = val;
            case 'segment',      opts.Segment = val;
            case 'plantype',     opts.PlanType = validate_plan_type(val);
            case 'ctlabel',      opts.CTLabel = normalize_ct_label(val);
            case 'hash',         opts.Hash    = char(val);
            case 'includers',    opts.IncludeRS    = logical(val);
            case 'includeethos', opts.IncludeEthos = logical(val);
            case 'includecbct',  opts.IncludeCBCT  = logical(val);
            otherwise
                error('load_recon_dose_data:UnknownOption', ...
                    'Unknown option: %s', char(name));
        end
    end
end

function m = validate_mode(v)
    m = lower(char(v));
    if ~ismember(m, {'total', 'set', 'single'})
        error('load_recon_dose_data:BadOption', ...
            'Mode must be one of: total | set | single.');
    end
end

function pt = validate_plan_type(v)
    pt = lower(char(v));
    if ~ismember(pt, {'any', 'adapted', 'reference'})
        error('load_recon_dose_data:BadOption', ...
            'PlanType must be one of: any | adapted | reference.');
    end
end

function lbl = normalize_ct_label(v)
    lbl = char(v);
    if strcmpi(lbl, 'any')
        lbl = 'any';
        return;
    end
    lbl = strrep(lbl, '-', '_');
    tok = regexp(lbl, '(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        lbl = sprintf('CT_%s', tok{1});
    end
end


%% =========================================================================
%  METHOD + HASH RESOLUTION
%% =========================================================================

function gm = resolve_gruneisen_method(patient_id, session, config)
    if isfield(config, 'gruneisen_method') && ~isempty(config.gruneisen_method)
        gm = char(config.gruneisen_method);
        return;
    end

    base = fullfile(config.working_dir, 'SimulationResults', patient_id, session);
    if ~isfolder(base)
        error('load_recon_dose_data:NoSimResults', ...
            'No SimulationResults directory: %s', base);
    end

    d     = dir(base);
    names = {d([d.isdir]).name};
    names = names(~ismember(names, {'.', '..'}));

    if isempty(names)
        error('load_recon_dose_data:NoMethodFolder', ...
            'No gruneisen-method subfolders under: %s', base);
    elseif numel(names) > 1
        error('load_recon_dose_data:AmbiguousMethod', ...
            ['Multiple gruneisen-method folders under %s: %s\n', ...
             'Set config.gruneisen_method to disambiguate.'], ...
            base, strjoin(names, ', '));
    end
    gm = names{1};
end

function hash8 = resolve_config_hash(sim_dir, mode, config, hash_override)
    if ~isempty(hash_override)
        hash8 = lower(hash_override);
        return;
    end

    % Try the CONFIG-derived hash first; use it only if its outputs exist.
    candidate = '';
    try
        candidate = compute_sim_config_hash(config);
    catch
        candidate = '';
    end
    if ~isempty(candidate) && hash_has_outputs(sim_dir, mode, candidate)
        hash8 = candidate;
        return;
    end

    % Fall back to scanning what is actually on disk.
    hashes = scan_hashes(sim_dir, mode);
    if isempty(hashes)
        error('load_recon_dose_data:NoReconFiles', ...
            'No reconstruction outputs found in %s for mode ''%s''.', ...
            sim_dir, mode);
    elseif numel(hashes) > 1
        error('load_recon_dose_data:AmbiguousHash', ...
            ['Multiple config hashes present in %s:\n  %s\n', ...
             'Pass ''Hash'' or a CONFIG that hashes to one of them ', ...
             '(see config_registry.json).'], ...
            sim_dir, strjoin(hashes, ', '));
    end
    hash8 = hashes{1};
end

function tf = hash_has_outputs(sim_dir, mode, hash8)
    if strcmp(mode, 'total')
        tf = isfile(fullfile(sim_dir, ...
            sprintf('total_recon_dose_%s.mat', hash8)));
    else
        tf = ~isempty(dir(fullfile(sim_dir, ...
            sprintf('*_recon_%s.mat', hash8))));
    end
end

function hashes = scan_hashes(sim_dir, mode)
    if strcmp(mode, 'total')
        files = dir(fullfile(sim_dir, 'total_recon_dose_*.mat'));
        pat   = '^total_recon_dose_([0-9a-fA-F]{8})\.mat$';
    else
        files = dir(fullfile(sim_dir, '*_recon_*.mat'));
        pat   = '_recon_([0-9a-fA-F]{8})\.mat$';
    end

    found = {};
    for i = 1:numel(files)
        tok = regexp(files(i).name, pat, 'tokens', 'once');
        if ~isempty(tok)
            found{end + 1} = lower(tok{1}); %#ok<AGROW>
        end
    end
    hashes = unique(found);
end


%% =========================================================================
%  TOTAL-MODE LOADING
%% =========================================================================

function out = load_total_branch(out, sim_dir, processed_dir, hash8, opts)
    recon_file = fullfile(sim_dir, sprintf('total_recon_dose_%s.mat', hash8));
    if ~isfile(recon_file)
        error('load_recon_dose_data:FileNotFound', ...
            'Total reconstruction not found: %s', recon_file);
    end

    data = load(recon_file);
    warn_hash_mismatch(data, hash8, recon_file);

    out.recon_dose = data.total_recon;
    if isfield(data, 'metadata')
        out.metadata = data.metadata;
    end

    if opts.IncludeRS
        out.rs_dose = load_total_rs_dose(processed_dir);
    end

    if opts.IncludeCBCT
        out.cbct = load_cbct_cache(processed_dir, {});
    end
end

function d = load_total_rs_dose(processed_dir)
    f = fullfile(processed_dir, 'total_rs_dose.mat');
    if ~isfile(f)
        error('load_recon_dose_data:FileNotFound', ...
            'total_rs_dose.mat not found: %s', f);
    end
    L = load(f);
    if isfield(L, 'total_rs_dose_sparse')
        d = reshape(full(L.total_rs_dose_sparse), L.total_rs_dose_dims);
    else
        d = L.total_rs_dose;
    end
end


%% =========================================================================
%  SET / SINGLE-MODE LOADING
%% =========================================================================

function out = load_set_branch(out, patient_id, session, config, ...
        sim_dir, processed_dir, hash8, opts)

    [field_index, metadata] = list_processed_field_doses(patient_id, session, config);
    out.metadata = metadata;

    % ---- Filter the index (cheap; parses filename tokens, no array loads) ----
    n_all   = numel(field_index);
    keep    = false(1, n_all);
    ptypes  = cell(1, n_all);
    ctlbls  = cell(1, n_all);
    for i = 1:n_all
        [pt, ct]  = parse_dose_filename_tokens(field_index(i).source_mat_filename);
        ptypes{i} = pt;
        ctlbls{i} = ct;

        ok = true;
        if ~isempty(opts.Beam)    && field_index(i).beam_index ~= opts.Beam, ok = false; end
        if ~isempty(opts.Segment) && field_index(i).segment    ~= opts.Segment, ok = false; end
        if ~strcmp(opts.PlanType, 'any') && ~strcmpi(pt, opts.PlanType), ok = false; end
        if ~strcmp(opts.CTLabel,  'any') && ~strcmpi(ct, opts.CTLabel),  ok = false; end
        keep(i) = ok;
    end
    sel = find(keep);

    % ---- Enforce mode-specific selection semantics ----
    if strcmp(opts.Mode, 'single')
        if isempty(opts.Beam)
            error('load_recon_dose_data:BeamRequired', ...
                'Mode ''single'' requires the ''Beam'' option.');
        end
        if isempty(sel)
            error('load_recon_dose_data:NoMatch', ...
                'No field matches Beam=%d with the given filters.', opts.Beam);
        end
        if numel(sel) > 1
            segs = [field_index(sel).segment];
            error('load_recon_dose_data:Ambiguous', ...
                ['Beam %d matches %d fields (segments: %s). ', ...
                 'Specify ''Segment'' to select one.'], ...
                opts.Beam, numel(sel), mat2str(segs));
        end
    elseif isempty(sel)
        error('load_recon_dose_data:NoMatch', ...
            'No fields match the given filters.');
    end

    % ---- Shared resources ----
    have_bm    = isfield(metadata, 'beam_metadata') && ~isempty(metadata.beam_metadata);
    cbct_cache = struct();
    if opts.IncludeCBCT
        cbct_cache = load_cbct_cache(processed_dir, {});
    end

    % ---- Build the per-field struct array ----
    fields_out = repmat(empty_field_struct(), 1, numel(sel));
    for n = 1:numel(sel)
        i    = sel(n);
        base = stem_of(field_index(i).source_mat_filename);

        recon_file = fullfile(sim_dir, sprintf('%s_recon_%s.mat', base, hash8));
        if ~isfile(recon_file)
            error('load_recon_dose_data:ReconMissing', ...
                ['Recon file missing for field %s at hash %s:\n  %s\n', ...
                 'Run the simulation for this CONFIG, or pass the matching ', ...
                 '''Hash''.'], ...
                field_index(i).source_mat_filename, hash8, recon_file);
        end

        rdata = load(recon_file);
        warn_hash_mismatch(rdata, hash8, recon_file);

        F                     = empty_field_struct();
        F.recon_dose          = rdata.recon_dose;
        F.source_mat_filename = field_index(i).source_mat_filename;
        F.recon_file          = recon_file;
        F.rtplan              = build_rtplan(metadata, have_bm, ...
            field_index(i).beam_index, field_index(i).segment, ...
            ptypes{i}, ctlbls{i});

        fd = [];
        if opts.IncludeRS
            fd         = load_field_dose_file(field_index(i).file);
            F.rs_dose  = fd.dose_Gy;
            F.rtplan   = merge_rtplan_from_fd(F.rtplan, fd);  % authoritative
        end

        if opts.IncludeCBCT
            lbl = ctlbls{i};
            if isempty(lbl) && ~isempty(fd) && isfield(fd, 'ct_label') ...
                    && ~isempty(fd.ct_label)
                lbl = strrep(char(fd.ct_label), '-', '_');
            end
            if isempty(lbl) || ~isfield(cbct_cache, lbl)
                error('load_recon_dose_data:CBCTNotFound', ...
                    'Could not resolve CBCT for field %s (ct_label=''%s'').', ...
                    field_index(i).source_mat_filename, lbl);
            end
            F.cbct = cbct_cache.(lbl);
        end

        fields_out(n) = F;
    end

    out.fields = fields_out;
end

function F = empty_field_struct()
    F = struct('recon_dose', [], 'rs_dose', [], 'cbct', [], ...
        'rtplan', struct(), 'source_mat_filename', '', 'recon_file', '');
end

function [plan_type, ct_label] = parse_dose_filename_tokens(name)
    plan_type = '';
    ct_label  = '';
    tok = regexp(name, ...
        '_(adapted|reference)(?:_CT_(\d+))?_B\d+_\d+\.mat$', 'tokens', 'once');
    if ~isempty(tok)
        plan_type = tok{1};
        if numel(tok) >= 2 && ~isempty(tok{2})
            ct_label = sprintf('CT_%s', tok{2});
        end
    end
end

function s = stem_of(name)
    [~, s, ~] = fileparts(name);
end


%% =========================================================================
%  RTPLAN STATISTICS
%% =========================================================================

function rt = build_rtplan(metadata, have_bm, beam_num, seg_num, plan_type, ct_label)
    rt = struct('beam_num', beam_num, 'seg_num', seg_num, ...
        'plan_type', plan_type, 'ct_label', ct_label, 'beam_name', '', ...
        'gantry_angle', NaN, 'meterset', NaN, ...
        'isocenter', [NaN NaN NaN], 'jaw_x', [NaN NaN], 'jaw_y', [NaN NaN]);

    if ~have_bm || isnan(beam_num)
        return;
    end

    bm = metadata.beam_metadata;
    for k = 1:numel(bm)
        if isfield(bm(k), 'beam_number') && ~isempty(bm(k).beam_number) ...
                && bm(k).beam_number == beam_num
            rt = copy_if(rt, bm(k), 'beam_name',    'beam_name');
            rt = copy_if(rt, bm(k), 'gantry_angle', 'gantry_angle');
            rt = copy_if(rt, bm(k), 'meterset',     'meterset');
            rt = copy_if(rt, bm(k), 'isocenter',    'isocenter');
            rt = copy_if(rt, bm(k), 'jaw_x',        'jaw_x');
            rt = copy_if(rt, bm(k), 'jaw_y',        'jaw_y');
            return;
        end
    end
end

function rt = merge_rtplan_from_fd(rt, fd)
    % The loaded field_dose carries values already matched to the RTPLAN at
    % processing time; prefer them over the beam_metadata lookup.
    rt = copy_if(rt, fd, 'gantry_angle', 'gantry_angle');
    rt = copy_if(rt, fd, 'meterset',     'meterset');
    rt = copy_if(rt, fd, 'isocenter',    'isocenter');
    rt = copy_if(rt, fd, 'jaw_x',        'jaw_x');
    rt = copy_if(rt, fd, 'jaw_y',        'jaw_y');
    rt = copy_if(rt, fd, 'plan_type',    'plan_type');
    rt = copy_if(rt, fd, 'ct_label',     'ct_label');
    rt = copy_if(rt, fd, 'beam_num',     'beam_num');
    rt = copy_if(rt, fd, 'seg_num',      'seg_num');
end

function s = copy_if(s, src, src_field, dst_field)
    if isfield(src, src_field) && ~isempty(src.(src_field))
        s.(dst_field) = src.(src_field);
    end
end


%% =========================================================================
%  CBCT GEOMETRY
%% =========================================================================

function cache = load_cbct_cache(processed_dir, needed_labels)
    specs = { 'CBCT1_resampled.mat', 'CBCT1_resampled'; ...
              'CBCT3_resampled.mat', 'CBCT3_resampled' };

    cache = struct();
    for k = 1:size(specs, 1)
        fpath = fullfile(processed_dir, specs{k, 1});
        if ~isfile(fpath)
            continue;
        end
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

    for i = 1:numel(needed_labels)
        if ~isfield(cache, needed_labels{i})
            error('load_recon_dose_data:CBCTNotFound', ...
                'Required CBCT geometry %s not found in %s.', ...
                needed_labels{i}, processed_dir);
        end
    end
end


%% =========================================================================
%  ETHOS TRUTH DOSE
%% =========================================================================

function ethos = load_ethos_truth_resampled(patient_id, session, config, metadata)
    treatment_site = 'Pancreas';
    if isfield(config, 'treatment_site') && ~isempty(config.treatment_site)
        treatment_site = char(config.treatment_site);
    end

    sct_dir  = fullfile(config.working_dir, 'EthosExports', patient_id, ...
        treatment_site, session, 'sct');
    rd_files = dir(fullfile(sct_dir, 'RTDOSE*.dcm'));
    if isempty(rd_files)
        error('load_recon_dose_data:EthosNotFound', ...
            'No RTDOSE*.dcm found in: %s', sct_dir);
    end

    rd_path = fullfile(sct_dir, rd_files(1).name);
    info    = dicominfo(rd_path);
    raw     = double(squeeze(dicomread(rd_path)));
    if isfield(info, 'DoseGridScaling')
        raw = raw * info.DoseGridScaling;
    end
    fprintf('  [ethos] Loaded truth dose: %s\n', rd_files(1).name);

    if isfield(metadata, 'dimensions') && ~isempty(metadata.dimensions)
        ethos = resample_to_grid(raw, metadata.dimensions(:)');
    else
        ethos = raw;
    end
end

function out_dose = resample_to_grid(dose, target_size)
    % Size-based trilinear resample, mirroring step3_analysis/resample_dose_to_grid
    % so recon, RS and ETHOS doses share one grid for gamma/SSIM.
    src = size(dose);
    if numel(src) < 3
        src(end + 1:3) = 1;
    end
    if isequal(src, target_size)
        out_dose = dose;
        return;
    end

    if exist('imresize3', 'file')
        out_dose = imresize3(double(dose), target_size, 'linear');
    else
        [Xi, Yi, Zi] = ndgrid( ...
            linspace(1, src(1), target_size(1)), ...
            linspace(1, src(2), target_size(2)), ...
            linspace(1, src(3), target_size(3)));
        out_dose = interp3(double(dose), Yi, Xi, Zi, 'linear', 0);
    end
    out_dose(out_dose < 0) = 0;
end


%% =========================================================================
%  SHARED HELPERS
%% =========================================================================

function warn_hash_mismatch(data, hash8, file)
    if isfield(data, 'config_hash') && ~strcmp(char(data.config_hash), hash8)
        warning('load_recon_dose_data:HashMismatch', ...
            ['File embeds config_hash %s but the resolved hash is %s:\n  %s\n', ...
             'This usually means the file was produced by a different CONFIG ', ...
             'that shares a filename.'], ...
            char(data.config_hash), hash8, file);
    end
end
