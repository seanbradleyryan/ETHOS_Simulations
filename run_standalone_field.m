function out = run_standalone_field(config)
%RUN_STANDALONE_FIELD Load one dose + CBCT, build the medium, run the engine.
%
%   out = run_standalone_field(config)
%
%   Single reusable setup path for every standalone / sweep / study / test
%   driver. It performs the load + medium-build + call sequence that each of
%   those scripts used to duplicate (~300 lines apiece), then delegates the
%   actual physics to run_single_field_simulation (the one simulation engine).
%   The driver keeps only its own analysis (gamma / plots / saving), built from
%   the arrays returned here.
%
%   WHAT IT DOES
%     1. Resolve + load the dose file (reuses load_field_dose_file for the
%        step15 field_dose format; falls back to total-dose formats).
%     2. Resolve + load the CBCT geometry (CBCT1_resampled / CBCT3_resampled).
%     3. Load plan beam_metadata (for determine_sensor_mask exclusion zones).
%     4. Build the acoustic medium via create_acoustic_medium and force the
%        coupling bath (outside body / couch) to the uniform medium -- the one
%        standalone-only medium step the engine does not do. Whole-grid
%        force-uniform toggles are left to the engine (config.force_uniform_*),
%        so nothing is double-applied.
%     5. Optionally build a reference-geometry reconstruction medium for a blind
%        reconstruction (config.blind_recon + config.recon_cbct_filename).
%     6. Call run_single_field_simulation with return_diagnostics = true and
%        normalize = false (drivers normalize themselves at the gamma stage).
%
%   INPUTS:
%       config - Simulation CONFIG (see get_default_config). Fields used here:
%           .working_dir, .patient_id, .session         - path resolution
%           .dose_filename / .dose_file_override         - which dose to load
%           .cbct_filename / .cbct_file_override         - forward geometry
%           .meterset                                    - fallback MU
%           .gruneisen_method, .tissue_tables            - medium build
%           .uniform_*                                   - coupling bath
%           .beam_metadata                               - optional; else loaded
%           (blind) .blind_recon, .recon_cbct_filename / .recon_cbct_file_override
%           plus all engine knobs (sensor, recon, pulse, downscale, ...).
%
%   OUTPUT (struct):
%       .recon_dose     - engine reconstruction (Gy), input-dose grid
%       .doseGrid       - original dose (Gy), body-masked, same grid
%       .reconPressure  - reconstructed pressure (Pa)
%       .p0             - initial pressure (Pa)
%       .sensor_mask    - logical sensor mask on the recon grid
%       .density        - reconstruction-geometry density background
%       .spacing        - [dx dy dz] mm actually used (post-downscale)
%       .gridSize       - size(recon_dose)
%       .sim_results    - full diagnostics struct from the engine
%       .field_dose     - the field_dose struct passed to the engine
%       .sct            - the CBCT geometry struct used (forward)
%       .beam_metadata  - beam metadata used for sensor placement
%
%   See also: run_single_field_simulation, get_default_config,
%             create_acoustic_medium, load_field_dose_file, load_processed_data

    if ~isstruct(config)
        error('run_standalone_field:InvalidInput', 'config must be a struct.');
    end

    processed_dir = fullfile(config.working_dir, 'RayStationFiles', ...
        config.patient_id, config.session, 'processed');

    %% ---- Load dose ----
    if isfield(config, 'dose_file_override') && ~isempty(config.dose_file_override)
        dose_filepath = config.dose_file_override;
    else
        dose_filepath = fullfile(processed_dir, config.dose_filename);
    end
    [doseGrid, fd] = load_dose_any(dose_filepath);

    %% ---- Load forward CBCT geometry ----
    if isfield(config, 'cbct_file_override') && ~isempty(config.cbct_file_override)
        cbct_filepath = config.cbct_file_override;
    else
        cbct_filepath = fullfile(processed_dir, config.cbct_filename);
    end
    sct = load_cbct_struct(cbct_filepath);

    % Spacing: prefer the CBCT field, fall back to the dose file's spacing.
    if isfield(sct, 'spacing') && ~isempty(sct.spacing)
        spacing_mm = sct.spacing(:)';
    elseif isfield(fd, 'spacing') && ~isempty(fd.spacing)
        spacing_mm = fd.spacing(:)';
        sct.spacing = spacing_mm;
    else
        error('run_standalone_field:NoSpacing', ...
            'Neither the CBCT struct nor the dose file carry a spacing.');
    end

    if ~isfield(sct, 'couchMask') || isempty(sct.couchMask)
        if isfield(sct, 'bodyMask')
            sct.couchMask = false(size(sct.bodyMask));
        end
    end
    if ~isfield(sct, 'origin') || isempty(sct.origin)
        sct.origin = [0, 0, 0];
    end

    gridSize = size(doseGrid);
    if ~isequal(gridSize, size(sct.cubeHU))
        error('run_standalone_field:GridMismatch', ...
            'Dose grid [%d %d %d] does not match CBCT grid [%d %d %d].', ...
            gridSize, size(sct.cubeHU));
    end

    % Body-mask the dose (standalone convention).
    if isfield(sct, 'bodyMask')
        doseGrid = doseGrid .* double(sct.bodyMask);
    end

    %% ---- Beam metadata (for determine_sensor_mask) ----
    beam_metadata = load_beam_metadata(config, processed_dir);

    %% ---- field_dose struct for the engine ----
    field_dose = struct();
    field_dose.dose_Gy    = doseGrid;
    field_dose.spacing    = spacing_mm;
    field_dose.dimensions = gridSize;
    field_dose.origin     = pick_field(fd, 'origin', sct.origin);
    field_dose.gantry_angle = pick_field(fd, 'gantry_angle', 0);
    % Meterset: file value if present and positive, else CONFIG fallback.
    ms = safe_get(config, 'meterset', 100);
    if isfield(fd, 'meterset') && ~isempty(fd.meterset) && fd.meterset > 0
        ms = fd.meterset;
    end
    field_dose.meterset = ms;
    field_dose = copy_if_present(fd, field_dose, {'isocenter', 'jaw_x', 'jaw_y'});

    %% ---- Build forward medium (+ coupling-bath override) ----
    medium = build_medium_with_bath(sct, config);

    %% ---- Optional blind reconstruction reference medium ----
    medium_recon = [];
    if safe_get(config, 'blind_recon', false)
        if isfield(config, 'recon_cbct_file_override') && ~isempty(config.recon_cbct_file_override)
            recon_cbct_filepath = config.recon_cbct_file_override;
        elseif isfield(config, 'recon_cbct_filename') && ~isempty(config.recon_cbct_filename)
            recon_cbct_filepath = fullfile(processed_dir, config.recon_cbct_filename);
        else
            error('run_standalone_field:NoReconCBCT', ...
                'config.blind_recon is true but no recon CBCT filename/override was given.');
        end
        sct_recon = load_cbct_struct(recon_cbct_filepath);
        if ~isfield(sct_recon, 'spacing') || isempty(sct_recon.spacing)
            sct_recon.spacing = spacing_mm;
        end
        medium_recon = build_medium_with_bath(sct_recon, config);
    end

    %% ---- Call the engine ----
    % return_diagnostics -> hand back the pressure/geometry arrays. normalize
    % is forced off here: the engine's 'normalize' means max-normalize, but the
    % standalone drivers normalize (lsq/max) themselves at the gamma stage.
    engine_config = config;
    engine_config.return_diagnostics = true;
    engine_config.normalize          = false;

    [recon_dose, sim_results] = run_single_field_simulation( ...
        field_dose, sct, medium, beam_metadata, engine_config, [], medium_recon);

    %% ---- Assemble output ----
    out = struct();
    out.recon_dose    = recon_dose;
    out.sim_results   = sim_results;
    out.field_dose    = field_dose;
    out.sct           = sct;
    out.beam_metadata = beam_metadata;
    out.gridSize      = size(recon_dose);
    out.spacing       = getfield_default(sim_results, 'spacing', spacing_mm);
    out.doseGrid      = getfield_default(sim_results, 'doseGrid', doseGrid);
    out.reconPressure = getfield_default(sim_results, 'reconPressure', []);
    out.p0            = getfield_default(sim_results, 'p0', []);
    out.sensor_mask   = getfield_default(sim_results, 'sensor_mask', []);
    out.density       = getfield_default(sim_results, 'density', []);
end


%% ======================== LOCAL HELPERS ========================

function [doseGrid, fd] = load_dose_any(dose_filepath)
%LOAD_DOSE_ANY Load a dose file in any of the step15 formats. Returns the dense
%  3D dose and the field_dose struct (empty struct for non-field_dose formats).
    if ~isfile(dose_filepath)
        error('run_standalone_field:DoseNotFound', 'Dose file not found: %s', dose_filepath);
    end
    dose_data = load(dose_filepath);
    fd = struct();
    if isfield(dose_data, 'field_dose')
        fd = load_field_dose_file(dose_filepath);   % reuse shared sparse handling
        doseGrid = double(fd.dose_Gy);
    elseif isfield(dose_data, 'total_rs_dose_sparse')
        if ~isfield(dose_data, 'total_rs_dose_dims')
            error('run_standalone_field:NoDims', 'total_rs_dose_dims missing for sparse total dose.');
        end
        doseGrid = double(reshape(full(dose_data.total_rs_dose_sparse), dose_data.total_rs_dose_dims));
    elseif isfield(dose_data, 'total_rs_dose')
        doseGrid = double(dose_data.total_rs_dose);
    elseif isfield(dose_data, 'dose_Gy')
        doseGrid = double(dose_data.dose_Gy);
    else
        dfn = fieldnames(dose_data);
        if numel(dfn) == 1
            doseGrid = double(dose_data.(dfn{1}));
        else
            error('run_standalone_field:UnknownDose', ...
                'Cannot auto-detect dose variable. Fields: %s', strjoin(dfn, ', '));
        end
    end
    if ndims(doseGrid) ~= 3
        error('run_standalone_field:BadDose', 'Dose data must be a 3D numeric array.');
    end
end


function sct = load_cbct_struct(cbct_filepath)
%LOAD_CBCT_STRUCT Load a resampled CBCT struct (CBCT1 or CBCT3) from a .mat.
    if ~isfile(cbct_filepath)
        error('run_standalone_field:CBCTNotFound', 'CBCT file not found: %s', cbct_filepath);
    end
    cbct_data = load(cbct_filepath);
    if isfield(cbct_data, 'CBCT1_resampled')
        sct = cbct_data.CBCT1_resampled;
    elseif isfield(cbct_data, 'CBCT3_resampled')
        sct = cbct_data.CBCT3_resampled;
    else
        error('run_standalone_field:NoCBCTStruct', ...
            'CBCT1_resampled / CBCT3_resampled not found in %s', cbct_filepath);
    end
    if ~isfield(sct, 'cubeHU')
        error('run_standalone_field:NoCubeHU', 'CBCT struct missing cubeHU: %s', cbct_filepath);
    end
end


function beam_metadata = load_beam_metadata(config, processed_dir)
%LOAD_BEAM_METADATA Prefer config.beam_metadata; else load metadata.mat.
    beam_metadata = [];
    if isfield(config, 'beam_metadata') && ~isempty(config.beam_metadata)
        beam_metadata = config.beam_metadata;
        return;
    end
    metadata_filepath = fullfile(processed_dir, 'metadata.mat');
    if isfile(metadata_filepath)
        try
            md = load(metadata_filepath, 'metadata');
            if isfield(md, 'metadata') && isfield(md.metadata, 'beam_metadata') ...
                    && ~isempty(md.metadata.beam_metadata)
                beam_metadata = md.metadata.beam_metadata;
            end
        catch ME
            warning('run_standalone_field:MetadataLoadFail', ...
                'Failed to load %s: %s', metadata_filepath, ME.message);
        end
    end
end


function medium = build_medium_with_bath(sct, config)
%BUILD_MEDIUM_WITH_BATH create_acoustic_medium + force the coupling bath
%  (outside body / couch) to the uniform medium. The whole-grid force-uniform
%  toggles are intentionally NOT applied here -- the engine owns those.
    medium = create_acoustic_medium(sct, config);
    ud = safe_get(config, 'uniform_density',     1000);
    uc = safe_get(config, 'uniform_sound_speed', 1540);
    ua = safe_get(config, 'uniform_alpha_coeff', 0);
    ug = safe_get(config, 'uniform_gruneisen',   1.0);
    if isfield(sct, 'bodyMask')
        outside = ~logical(sct.bodyMask);
        if isfield(sct, 'couchMask') && ~isempty(sct.couchMask)
            outside = outside | logical(sct.couchMask);
        end
        medium.density(outside)     = ud;
        medium.sound_speed(outside) = uc;
        if numel(medium.alpha_coeff) > 1
            medium.alpha_coeff(outside) = ua;
        end
        medium.gruneisen(outside)   = ug;
    end
end


function v = safe_get(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end

function v = getfield_default(s, f, d)
    if isstruct(s) && isfield(s, f); v = s.(f); else; v = d; end
end

function v = pick_field(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end

function dst = copy_if_present(src, dst, names)
    for i = 1:numel(names)
        if isfield(src, names{i}) && ~isempty(src.(names{i}))
            dst.(names{i}) = src.(names{i});
        end
    end
end
