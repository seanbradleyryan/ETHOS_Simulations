function results = verify_pipeline_simulate(patient_id, session, config, config_hash)
%VERIFY_PIPELINE_SIMULATE Check that pipeline_simulate ran correctly; save debug plots.
%
%   results = verify_pipeline_simulate()
%       Standalone: run the file directly (F5). Edit STANDALONE DEFAULTS below.
%   results = verify_pipeline_simulate(patient_id, session, CONFIG, CONFIG_HASH)
%       Called by pipeline_simulate at the end of a run (CONFIG.run_verify).
%
%   PURPOSE:
%   Confirms that the Step 2 k-Wave simulations finished and that their outputs
%   make physical sense, and saves plots for a human to review.
%
%   TESTS (each prints PASS, FAIL, INFO (figure / numbers to review), SKIP or
%   ERROR (the test itself crashed); FAIL/ERROR lines are repeated at the end):
%     1  Completion: beam/segment pairs in the RTPLAN x {CT_1, CT_3} = processed
%        doses = recon files for this config hash, and the saved total recon
%        equals the sum of the per-field recons.
%     2  Sensor placement on CBCT1 and CBCT3: axial/coronal/sagittal slices
%        through the sensor center.                                    [figure]
%     3  (added) The sensor is clear of the body and couch on BOTH CBCTs. The
%        plan sensor is placed once on CBCT1 and reused for every CT_3 field, so
%        an anatomy change can put it inside the CT_3 patient.
%     4  Random beam/segment pairs: RS truth, recon and gamma map with pass rate,
%        on CT_1 and CT_3.                                            [figures]
%     5  Totals per CT: ETHOS truth, RS total and recon total, with gamma maps,
%        pass rates and gamma histograms.                             [figures]
%     6  Sensor signal after each processing stage (raw, pulse convolved,
%        band-limited, noise added, deconvolved), averaged over all sensor
%        points, on a log axis. The pipeline does not save these signals, so
%        this re-runs ONE forward k-Wave simulation (no reconstruction). [figure]
%     7  CBCT1 vs CBCT3: both images and their HU difference; FAIL when they are
%        the same image.                                               [figure]
%     8  (added) Silent failures: run_single_field_simulation returns an
%        all-zero recon when k-Wave fails, and pipeline_simulate hides its
%        warning (evalc), so every recon is checked for zero / NaN / wrong size.
%        Also plots each field's recon/truth dose ratio and SNR so outliers
%        stand out.                                                    [figure]
%
%   INPUTS:
%       patient_id  - char, e.g. '1194203'
%       session     - char, e.g. 'Session_1'
%       config      - pipeline_simulate CONFIG. Uses .working_dir,
%                     .gruneisen_method, .treatment_site, the simulation fields
%                     (test 6), and the Step 2.5 gamma settings gamma_dose_pct,
%                     gamma_dist_mm, gamma_dose_cutoff_pct and metrics_normalize
%                     (so pass rates match Step 2.5). Optional, default in []:
%           .verify_num_pairs              [5]          random pairs (test 4)
%           .verify_random_seed            ['shuffle']  integer = same pairs each run
%           .verify_ct_window_hu           [-500 500]   CT display window (HU)
%           .verify_hu_diff_window         [300]        HU difference color range
%           .verify_min_cbct_hu_diff       [1]          test 7: mean |dHU| in body below this = same image
%           .verify_sensor_in_body_tol_pct [1]          test 3: % of sensor voxels allowed in body/couch
%           .verify_ratio_outlier_factor   [2]          test 8: flag ratios beyond x or 1/x of the median
%       config_hash - 8-char hash of the run. '' = the only hash found on disk;
%                     the simulation fields are then read from config_registry.json.
%
%   OUTPUTS:
%       results - struct: .patient_id, .session, .config_hash, .out_dir and
%                 .tests (struct array with .name, .status, .detail)
%       PNG figures in SimulationResults/<id>/<session>/<method>/simulation_debug/
%       named <plot>_<hash>.png. They are also shown when MATLAB has a desktop.
%
%   ALGORITHM:
%       1. Resolve the config hash; index the processed doses; load both CBCTs.
%       2. One parfor pass loads every RS field dose and its recon once, checks
%          the recon and adds both to per-CT sums (used by tests 1, 2, 5, 8).
%       3. Tests 1-8, each in its own try/catch, then a summary.
%
%   EXAMPLE:
%       verify_pipeline_simulate('1194203', 'Session_1', CONFIG, CONFIG_HASH);
%
%   DEPENDENCIES:
%       utils/ (get_default_config, hashtoconfig, list_processed_field_doses,
%       load_field_dose_file, least_squares_gain, compute_gamma, CalcGamma,
%       create_acoustic_medium), run_single_field_simulation + k-Wave (test 6),
%       Image Processing Toolbox (dicominfo, dicomread), Parallel Computing
%       Toolbox (optional; parfor runs serially without it).
%
%   See also: pipeline_simulate, verify_pipeline_compress_output,
%             run_single_field_simulation, step25_segment_metrics

    %% ========================= STANDALONE DEFAULTS ===========================
    % Used only when this file is run directly with no inputs.
    if nargin == 0
        addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'utils')));
        patient_id  = '1194203';
        session     = 'Session_1';
        config      = get_default_config();
        config_hash = '';   % '' = the only config hash on disk (errors if several)
    elseif nargin ~= 4
        error('verify_pipeline_simulate:InvalidInput', ...
            'Call with no inputs, or with (patient_id, session, config, config_hash).');
    end

    %% ========================= INPUT VALIDATION ==============================
    if ~ischar(patient_id) && ~isstring(patient_id)
        error('verify_pipeline_simulate:InvalidInput', ...
            'patient_id must be a string or character array. Received: %s', class(patient_id));
    end
    patient_id = char(patient_id);
    if ~ischar(session) && ~isstring(session)
        error('verify_pipeline_simulate:InvalidInput', ...
            'session must be a string or character array. Received: %s', class(session));
    end
    session     = char(session);
    config_hash = char(config_hash);

    %% ========================= SETTINGS ======================================
    numPairs         = config_value(config, 'verify_num_pairs', 5);
    randomSeed       = config_value(config, 'verify_random_seed', 'shuffle');
    ctWindowHu       = config_value(config, 'verify_ct_window_hu', [-500 500]);
    huDiffWindow     = config_value(config, 'verify_hu_diff_window', 300);
    minCbctHuDiff    = config_value(config, 'verify_min_cbct_hu_diff', 1);
    sensorBodyTolPct = config_value(config, 'verify_sensor_in_body_tol_pct', 1);
    ratioOutlier     = config_value(config, 'verify_ratio_outlier_factor', 2);
    totalRelTol      = 1e-6;   % test 1: saved total vs summed recons may differ by rounding only
    treatmentSite    = config_value(config, 'treatment_site', 'Pancreas');

    % Gamma criteria, cutoff and recon gain as in Step 2.5, so pass rates are comparable
    gammaPct       = config_value(config, 'gamma_dose_pct', 3);
    gammaDta       = config_value(config, 'gamma_dist_mm', 3);
    gammaCriteria  = {gammaPct, gammaDta, sprintf('%g%%/%gmm', gammaPct, gammaDta)};
    cutoffFraction = config_value(config, 'gamma_dose_cutoff_pct', 10) / 100;
    normalizeRecon = config_value(config, 'metrics_normalize', true);   % least-squares recon gain

    ctLabels  = {'CT_1', 'CT_3'};
    cbctNames = {'CBCT1', 'CBCT3'};
    viewNames = {'axial', 'coronal', 'sagittal'};
    % Gamma colors: blue -> green for 0..1 (pass), yellow -> red for 1..2 (fail)
    gammaCmap = [winter(128); flipud(autumn(128))];
    % HU difference colors: blue (negative) - white (0) - red (positive)
    divergingCmap = [linspace(0, 1, 128)', linspace(0, 1, 128)', ones(128, 1); ...
                     ones(128, 1), linspace(1, 0, 128)', linspace(1, 0, 128)'];
    % Show figures only when there is a screen; always save them
    showFigures = usejava('desktop');   % false in -batch / -nodisplay runs

    %% ========================= PATHS + CONFIG HASH ===========================
    rsDir        = fullfile(config.working_dir, 'RayStationFiles', patient_id, session);
    processedDir = fullfile(rsDir, 'processed');
    simDir       = fullfile(config.working_dir, 'SimulationResults', patient_id, session, ...
        config.gruneisen_method);
    % Folders that may hold the RTPLAN / ETHOS RTDOSE DICOMs, in search order
    dicomDirs = {rsDir, ...
        fullfile(config.working_dir, 'Raystation_Input', patient_id, session), ...
        fullfile(config.working_dir, 'EthosExports', patient_id, treatmentSite, session, 'sct')};

    if isempty(config_hash)
        % Standalone: use the one hash present in the recon file names
        reconFiles = dir(fullfile(simDir, '*_recon_*.mat'));
        hashes = {};
        for i = 1:numel(reconFiles)
            tok = regexp(reconFiles(i).name, '_recon_([0-9a-f]{8})\.mat$', 'tokens', 'once');
            if ~isempty(tok)
                hashes{end + 1} = tok{1}; %#ok<AGROW>
            end
        end
        hashes = unique(hashes);
        if numel(hashes) ~= 1
            error('verify_pipeline_simulate:AmbiguousHash', ...
                'Expected one config hash in %s, found %d (%s). Set config_hash.', ...
                simDir, numel(hashes), strjoin(hashes(:)', ', '));
        end
        config_hash = hashes{1};
        % Take the simulation fields that produced this hash from config_registry.json
        % ('__unset__' = the field was absent in that run, so keep the default)
        registryConfig = hashtoconfig(config_hash, config);
        names = fieldnames(registryConfig);
        for i = 1:numel(names)
            if ~isequal(registryConfig.(names{i}), '__unset__')
                config.(names{i}) = registryConfig.(names{i});
            end
        end
    end

    outDir = fullfile(simDir, 'simulation_debug');
    if ~isfolder(outDir)
        mkdir(outDir);
    end

    fprintf('\n=========================================================\n');
    fprintf('  Verify pipeline_simulate: %s / %s (config %s)\n', patient_id, session, config_hash);
    fprintf('  Figures: %s\n', outDir);
    fprintf('=========================================================\n');

    tests = struct('name', {}, 'status', {}, 'detail', {});

    %% ========================= SHARED INPUTS =================================
    [fieldIndex, metadata] = list_processed_field_doses(patient_id, session, config);
    gridDims   = metadata.dimensions(:)';   % [ny nx nz] = (rows, cols, slices)
    spacing    = metadata.spacing(:)';      % [dx dy dz] mm
    origin     = metadata.origin(:)';       % [x y z] mm of voxel (1,1,1)
    gammaWidth = [spacing(2) spacing(1) spacing(3)];   % CalcGamma widths follow the array: rows = y

    cbcts  = cell(1, 2);
    loaded = load(fullfile(processedDir, 'CBCT1_resampled.mat'), 'CBCT1_resampled');
    cbcts{1} = loaded.CBCT1_resampled;
    loaded = load(fullfile(processedDir, 'CBCT3_resampled.mat'), 'CBCT3_resampled');
    cbcts{2} = loaded.CBCT3_resampled;
    for c = 1:2
        if ~isfield(cbcts{c}, 'couchMask') || isempty(cbcts{c}.couchMask)
            cbcts{c}.couchMask = false(size(cbcts{c}.bodyMask));
        end
    end

    % Plan type, CT label and expected recon file of every processed dose.
    % Names: dose_<id>_<session>_<plan>_<CT_n>_B<beam>_<seg>.mat
    nFields   = numel(fieldIndex);
    planType  = cell(nFields, 1);   % 'reference' | 'adapted' | '' (name not parsable)
    ctLabel   = cell(nFields, 1);   % 'CT_1' | 'CT_3' | ''
    reconName = cell(nFields, 1);
    reconPath = cell(nFields, 1);
    for i = 1:nFields
        tok = regexp(fieldIndex(i).source_mat_filename, ...
            '_(adapted|reference)_(CT_\d+)_B\d+_\d+\.mat$', 'tokens', 'once');
        if isempty(tok)
            tok = {'', ''};
        end
        planType{i}  = tok{1};
        ctLabel{i}   = tok{2};
        [~, stem]    = fileparts(fieldIndex(i).source_mat_filename);
        reconName{i} = sprintf('%s_recon_%s.mat', stem, config_hash);
        reconPath{i} = fullfile(simDir, reconName{i});
    end

    % Pair the CT_1 and CT_3 dose of every beam/segment.
    % pairCt1(k) / pairCt3(k) = field index of that CT's dose, 0 when missing.
    pairKeys = {};
    pairCt1  = [];
    pairCt3  = [];
    for i = 1:nFields
        if isempty(planType{i})
            continue;
        end
        key = sprintf('%s_B%d_S%d', planType{i}, fieldIndex(i).beam_index, fieldIndex(i).segment);
        k = find(strcmp(pairKeys, key), 1);
        if isempty(k)
            pairKeys{end + 1} = key; %#ok<AGROW>
            pairCt1(end + 1)  = 0;   %#ok<AGROW>
            pairCt3(end + 1)  = 0;   %#ok<AGROW>
            k = numel(pairKeys);
        end
        if strcmp(ctLabel{i}, 'CT_1')
            pairCt1(k) = i;
        elseif strcmp(ctLabel{i}, 'CT_3')
            pairCt3(k) = i;
        end
    end

    %% ========================= ONE PASS OVER EVERY FIELD =====================
    % Loads each RS field dose and its recon exactly once. The four sums are
    % parfor "reduction" variables: each worker keeps its own sum and MATLAB
    % adds them together at the end.
    fprintf('\nLoading %d field doses and their recons (one pass)...\n', nFields);
    hasRecon    = false(nFields, 1);   % recon file on disk
    reconOk     = false(nFields, 1);   % ...and right size, finite, not all zero
    loadFailed  = false(nFields, 1);   % a dose or recon file could not be read
    doseRatio   = nan(nFields, 1);     % sum(recon) / sum(truth) inside the truth's cutoff region
    fieldSnr    = nan(nFields, 1);     % SNR the field saw (noise_stats saved with the recon)
    reconSumCt1 = zeros(gridDims);
    reconSumCt3 = zeros(gridDims);
    rsSumCt1    = zeros(gridDims);
    rsSumCt3    = zeros(gridDims);
    parfor i = 1:nFields
        warning('off', 'MATLAB:load:variableNotFound');   % older recons have no noise_stats
        onDisk = isfile(reconPath{i});
        hasRecon(i) = onDisk;
        try
            fd = load_field_dose_file(fieldIndex(i).file);
            rs = double(fd.dose_Gy);
            if strcmp(ctLabel{i}, 'CT_1')
                rsSumCt1 = rsSumCt1 + rs;
            elseif strcmp(ctLabel{i}, 'CT_3')
                rsSumCt3 = rsSumCt3 + rs;
            end
            if onDisk
                loadedRecon = load(reconPath{i}, 'recon_dose', 'noise_stats');
                recon = double(loadedRecon.recon_dose);
                if isfield(loadedRecon, 'noise_stats')
                    fieldSnr(i) = loadedRecon.noise_stats.snr;
                end
                ok = isequal(size(recon), size(rs)) && all(isfinite(recon(:))) && max(recon(:)) > 0;
                reconOk(i) = ok;
                if ok
                    inRegion = rs >= cutoffFraction * max(rs(:));
                    doseRatio(i) = sum(recon(inRegion)) / sum(rs(inRegion));
                    if strcmp(ctLabel{i}, 'CT_1')
                        reconSumCt1 = reconSumCt1 + recon;
                    elseif strcmp(ctLabel{i}, 'CT_3')
                        reconSumCt3 = reconSumCt3 + recon;
                    end
                end
            end
        catch
            loadFailed(i) = true;
        end
    end
    reconSums = {reconSumCt1, reconSumCt3};
    rsSums    = {rsSumCt1, rsSumCt3};

    %% ---------------- TEST 1: completion ----------------
    testName = '1. Completion: plan pairs x 2 CTs = doses = recons';
    fprintf('\n[TEST 1] %s\n', testName);
    try
        % Expected doses: every segment of every RTPLAN beam, on both CTs. Step 0.6
        % makes one segment per control-point interval, numbered from 0.
        planKeys = {};
        plans = unique(planType(~cellfun(@isempty, planType)));
        for p = 1:numel(plans)
            rtplan = dicominfo(find_dicom(sprintf('RTPLAN_%s_adjusted_mlc.dcm', plans{p}), dicomDirs));
            beamItems = fieldnames(rtplan.BeamSequence);
            for b = 1:numel(beamItems)
                beam = rtplan.BeamSequence.(beamItems{b});
                numSegments = numel(fieldnames(beam.ControlPointSequence)) - 1;
                for seg = 0:numSegments - 1
                    for c = 1:2
                        planKeys{end + 1} = sprintf('%s_%s_B%d_S%d', plans{p}, ctLabels{c}, ...
                            beam.BeamNumber, seg); %#ok<AGROW>
                    end
                end
            end
        end
        doseKeys = cell(1, nFields);
        for i = 1:nFields
            doseKeys{i} = sprintf('%s_%s_B%d_S%d', planType{i}, ctLabel{i}, ...
                fieldIndex(i).beam_index, fieldIndex(i).segment);
        end
        missingDoses  = setdiff(planKeys, doseKeys);
        extraDoses    = setdiff(doseKeys, planKeys);
        missingRecons = reconName(~hasRecon);
        % Recon files for this hash that no current processed dose produced
        hashRecons  = dir(fullfile(simDir, sprintf('*_recon_%s.mat', config_hash)));
        strayRecons = setdiff({hashRecons.name}, reconName);

        totalFile = fullfile(simDir, sprintf('total_recon_dose_%s.mat', config_hash));
        totalErr  = Inf;
        totalMsg  = 'total_recon_dose file MISSING';
        if isfile(totalFile)
            loaded = load(totalFile, 'total_recon');
            summed = reconSums{1} + reconSums{2};
            totalMsg = 'saved total recon has a different size than the field recons';
            if isequal(size(loaded.total_recon), size(summed))
                totalErr = max(abs(loaded.total_recon(:) - summed(:))) / max(abs(loaded.total_recon(:)));
                totalMsg = sprintf('max|saved total - summed recons| / max = %.1e', totalErr);
            end
        end

        passed = isempty(missingDoses) && isempty(extraDoses) && isempty(missingRecons) ...
            && isempty(strayRecons) && totalErr < totalRelTol;
        tests = add_result(tests, testName, pass_fail(passed), sprintf( ...
            ['RTPLAN: %d beam/segment pairs -> %d doses expected; found %d doses ' ...
             '(%d missing, %d not in plan); %d recons (%d missing, %d stray); %s'], ...
            numel(planKeys) / 2, numel(planKeys), nFields, numel(missingDoses), ...
            numel(extraDoses), nnz(hasRecon), numel(missingRecons), numel(strayRecons), totalMsg));
        print_examples('missing dose', missingDoses);
        print_examples('dose not in plan', extraDoses);
        print_examples('missing recon', missingRecons);
        print_examples('stray recon', strayRecons);
    catch ME
        tests = add_result(tests, testName, 'ERROR', ME.message);
    end

    %% ---------------- TESTS 2 + 3: sensor placement / sensor clear of the body ----------------
    testName2 = '2. Sensor placement on CBCT1 and CBCT3';
    testName3 = '3. (added) Sensor clear of body and couch on both CBCTs';
    fprintf('\n[TESTS 2-3] %s\n', testName2);
    sensorFile = fullfile(simDir, sprintf('sensor_mask_%s.mat', config_hash));
    try
        if ~isfile(sensorFile)
            msg = sprintf(['no sensor_mask_%s.mat: only the determine_sensor_mask* methods save ' ...
                'a plan sensor (this run: %s)'], config_hash, config.sensor_placement_method);
            tests = add_result(tests, testName2, 'SKIP', msg);
            tests = add_result(tests, testName3, 'SKIP', msg);
        else
            loaded     = load(sensorFile, 'precomputed_sensor');
            sensorInfo = loaded.precomputed_sensor.sensor_info;
            sensor     = loaded.precomputed_sensor.sensor_mask > 0;
            % determine_sensor_mask may pad the grid with water so the array fits.
            % grid_pad y/x/z = array dims 1/2/3. Put each CBCT on that bigger grid.
            gp     = sensorInfo.grid_pad;
            offset = [gp.y_pre, gp.x_pre, gp.z_pre];
            % Summed plan dose: what the sensor placement steered away from
            planDose   = rsSums{1} + rsSums{2};
            planDose10 = embed_in_grid(planDose >= cutoffFraction * max(planDose(:)), ...
                size(sensor), offset, false);
            [rows, cols, slices] = ind2sub(size(sensor), find(sensor));
            center = round([mean(rows), mean(cols), mean(slices)]);

            fig    = new_figure('Sensor placement', showFigures);
            inBody = zeros(1, 2);
            for c = 1:2
                hu    = embed_in_grid(cbcts{c}.cubeHU, size(sensor), offset, -1000);
                body  = embed_in_grid(cbcts{c}.bodyMask, size(sensor), offset, false);
                couch = embed_in_grid(cbcts{c}.couchMask, size(sensor), offset, false);
                inBody(c) = nnz(sensor & (body | couch));
                for v = 1:3
                    ax = subplot(2, 3, (c - 1) * 3 + v, 'Parent', fig);
                    show_ct(ax, hu, v, center, spacing, ctWindowHu);
                    add_overlay(ax, double(sensor), sensor, v, center, [0 1], [1 0 0]);
                    add_contour(ax, body, v, center, 'g');
                    add_contour(ax, planDose10, v, center, 'c');
                    title(ax, sprintf('%s %s', cbctNames{c}, viewNames{v}));
                end
            end
            sgtitle(fig, sprintf(['%s %s | sensor (red), body (green), plan dose >= %g%% (cyan) | ' ...
                'slices through sensor center [row %d, col %d, slice %d] | sensor in body: ' ...
                'CBCT1 %d, CBCT3 %d voxels'], patient_id, session, 100 * cutoffFraction, ...
                center, inBody), 'Interpreter', 'none');
            save_figure(fig, outDir, ['sensor_placement_' config_hash], showFigures);

            tests = add_result(tests, testName2, 'INFO', sprintf( ...
                ['%d sensor voxels, center %s mm, tilt %.1f deg, placement_valid = %d, ' ...
                 'grid expanded = %d; figure sensor_placement_%s.png'], nnz(sensor), ...
                mat2str(sensorInfo.sensor_center_mm, 4), sensorInfo.tilt_angle_deg, ...
                sensorInfo.placement_valid, gp.expanded, config_hash));
            pctInBody = 100 * inBody / nnz(sensor);
            tests = add_result(tests, testName3, pass_fail(all(pctInBody <= sensorBodyTolPct)), sprintf( ...
                'sensor voxels inside body/couch: CBCT1 %d (%.2f%%), CBCT3 %d (%.2f%%); tolerance %g%%', ...
                inBody(1), pctInBody(1), inBody(2), pctInBody(2), sensorBodyTolPct));
        end
    catch ME
        tests = add_result(tests, testName2, 'ERROR', ME.message);
        tests = add_result(tests, testName3, 'ERROR', ME.message);
    end

    %% ---------------- TEST 4: random beam/segment pairs ----------------
    testName = sprintf('4. Random pairs: truth, recon, %s gamma on CT_1 and CT_3', gammaCriteria{3});
    fprintf('\n[TEST 4] %s\n', testName);
    picks = [];
    try
        rng(randomSeed);
        complete = find(pairCt1 > 0 & pairCt3 > 0);
        complete = complete(reconOk(pairCt1(complete)) & reconOk(pairCt3(complete)));
        picks    = complete(randperm(numel(complete), min(numPairs, numel(complete))));
        passRates = nan(numel(picks), 2);   % columns: CT_1, CT_3
        for n = 1:numel(picks)
            k = picks(n);
            pairFields = [pairCt1(k), pairCt3(k)];
            fig = new_figure(['Pair ' pairKeys{k}], showFigures);
            for c = 1:2
                fd     = load_field_dose_file(fieldIndex(pairFields(c)).file);
                truth  = double(fd.dose_Gy);
                loaded = load(reconPath{pairFields(c)}, 'recon_dose');
                gain   = 1;
                if normalizeRecon
                    gain = least_squares_gain(truth, loaded.recon_dose);   % same scaling as Step 2.5
                end
                recon = gain * double(loaded.recon_dose);
                g = compute_gamma(truth, recon, gammaWidth, 'Criteria', gammaCriteria, ...
                    'Cutoff', cutoffFraction);
                passRates(n, c) = g.pass_rates(1);

                % Both rows: axial slice through the CT_1 truth max, CT_1 color scale
                if c == 1
                    [~, iMax] = max(truth(:));
                    [row, col, slc] = ind2sub(size(truth), iMax);
                    center  = [row, col, slc];
                    doseMax = max(truth(:));
                end
                hu = cbcts{c}.cubeHU;
                ax = subplot(2, 3, (c - 1) * 3 + 1, 'Parent', fig);
                show_ct(ax, hu, 1, center, spacing, ctWindowHu);
                add_overlay(ax, truth, truth >= cutoffFraction * doseMax, 1, center, [0 doseMax], jet(256));
                colorbar(ax);
                title(ax, sprintf('%s RS truth (Gy)', ctLabels{c}), 'Interpreter', 'none');

                ax = subplot(2, 3, (c - 1) * 3 + 2, 'Parent', fig);
                show_ct(ax, hu, 1, center, spacing, ctWindowHu);
                add_overlay(ax, recon, recon >= cutoffFraction * doseMax, 1, center, [0 doseMax], jet(256));
                colorbar(ax);
                title(ax, sprintf('%s recon x %.3g (Gy)', ctLabels{c}, gain), 'Interpreter', 'none');

                ax = subplot(2, 3, (c - 1) * 3 + 3, 'Parent', fig);
                show_ct(ax, hu, 1, center, spacing, ctWindowHu);
                add_overlay(ax, g.maps{1}, g.eval_mask, 1, center, [0 2], gammaCmap);
                colorbar(ax);
                title(ax, sprintf('%s gamma %s: %.1f%% pass', ctLabels{c}, gammaCriteria{3}, ...
                    passRates(n, c)), 'Interpreter', 'none');
            end
            sgtitle(fig, sprintf('%s %s | %s | axial slice %d (CT_1 truth max) | recon x least-squares gain: %d', ...
                patient_id, session, pairKeys{k}, center(3), normalizeRecon), 'Interpreter', 'none');
            save_figure(fig, outDir, sprintf('pair_%s_%s', pairKeys{k}, config_hash), showFigures);
            fprintf('    %s: gamma pass CT_1 %.1f%%, CT_3 %.1f%%\n', pairKeys{k}, passRates(n, 1), passRates(n, 2));
        end

        if isempty(picks)
            tests = add_result(tests, testName, 'FAIL', ...
                'no beam/segment has a valid CT_1 AND CT_3 recon');
        else
            tests = add_result(tests, testName, 'INFO', sprintf( ...
                '%d pairs, gamma pass: CT_1 mean %.1f%% (min %.1f%%), CT_3 mean %.1f%% (min %.1f%%); figures pair_*.png', ...
                numel(picks), mean(passRates(:, 1)), min(passRates(:, 1)), ...
                mean(passRates(:, 2)), min(passRates(:, 2))));
        end
    catch ME
        tests = add_result(tests, testName, 'ERROR', ME.message);
    end

    %% ---------------- TEST 5: totals (ETHOS, RS, recon) ----------------
    testName = sprintf('5. Totals per CT: ETHOS vs RS vs recon, %s gamma', gammaCriteria{3});
    fprintf('\n[TEST 5] %s\n', testName);
    try
        plans = unique(planType(~cellfun(@isempty, planType)));
        if numel(plans) ~= 1
            tests = add_result(tests, testName, 'SKIP', sprintf( ...
                'expected one plan type, found %d (%s); the totals would mix plans', ...
                numel(plans), strjoin(plans(:)', ', ')));
        else
            % ETHOS truth in Gy per fraction on the dose grid (same method as
            % verify_pipeline_compress_output test 4). RS doses are ONE fraction;
            % a PLAN-summed RTDOSE covers every fraction, so divide it.
            rtdoseFile = find_dicom(sprintf('RTDOSE_%s.dcm', plans{1}), dicomDirs);
            info       = dicominfo(rtdoseFile);
            ethosRaw   = double(squeeze(dicomread(rtdoseFile))) * info.DoseGridScaling;
            nFractions = 1;
            if strcmpi(info.DoseSummationType, 'PLAN')
                rtplan = dicominfo(find_dicom(sprintf('RTPLAN_%s_adjusted_mlc.dcm', plans{1}), dicomDirs));
                nFractions = double(rtplan.FractionGroupSequence.Item_1.NumberOfFractionsPlanned);
            end
            ethosRaw = ethosRaw / nFractions;
            % Resample by patient position (mm). DICOM PixelSpacing = [row (y), column (x)].
            ethosX = info.ImagePositionPatient(1) + (0:size(ethosRaw, 2) - 1) * info.PixelSpacing(2);
            ethosY = info.ImagePositionPatient(2) + (0:size(ethosRaw, 1) - 1) * info.PixelSpacing(1);
            ethosZ = info.GridFrameOffsetVector(:)';
            if ethosZ(1) == 0   % offsets relative to the first frame (the usual case)
                ethosZ = ethosZ + info.ImagePositionPatient(3);
            end
            [qx, qy, qz] = meshgrid(origin(1) + (0:gridDims(2) - 1) * spacing(1), ...
                                    origin(2) + (0:gridDims(1) - 1) * spacing(2), ...
                                    origin(3) + (0:gridDims(3) - 1) * spacing(3));
            ethosOnGrid = interp3(ethosX, ethosY, ethosZ, ethosRaw, qx, qy, qz, 'linear', 0);
            clear qx qy qz ethosRaw;

            shortNames = {'ETHOS', 'RS', 'recon'};
            compRef    = [1 2 1];   % gamma comparisons: reference volume ...
            compTgt    = [2 3 3];   % ... vs target volume (indices into shortNames)
            passRates  = nan(2, 3);
            for c = 1:2
                % Zero ETHOS outside this CBCT's body / inside the couch, as Step 1.5 did to the RS doses
                ethos = ethosOnGrid;
                ethos(~(cbcts{c}.bodyMask & ~cbcts{c}.couchMask)) = 0;
                gain = 1;
                if normalizeRecon
                    gain = least_squares_gain(rsSums{c}, reconSums{c});
                end
                volumes  = {ethos, rsSums{c}, gain * reconSums{c}};
                volNames = {sprintf('ETHOS / %d fx', nFractions), ...
                            sprintf('RS total %s', ctLabels{c}), ...
                            sprintf('recon total %s x %.3g', ctLabels{c}, gain)};
                gammas = cell(1, 3);
                for m = 1:3
                    gammas{m} = compute_gamma(volumes{compRef(m)}, volumes{compTgt(m)}, gammaWidth, ...
                        'Criteria', gammaCriteria, 'Cutoff', cutoffFraction);
                    passRates(c, m) = gammas{m}.pass_rates(1);
                end

                % Axial slice through the ETHOS max, one color scale for all three doses
                [~, iMax] = max(ethos(:));
                [row, col, slc] = ind2sub(size(ethos), iMax);
                center  = [row, col, slc];
                doseMax = max([max(volumes{1}(:)), max(volumes{2}(:)), max(volumes{3}(:))]);
                hu      = cbcts{c}.cubeHU;
                fig     = new_figure(['Totals ' ctLabels{c}], showFigures);
                for m = 1:3   % row 1: the three doses
                    ax = subplot(3, 3, m, 'Parent', fig);
                    show_ct(ax, hu, 1, center, spacing, ctWindowHu);
                    add_overlay(ax, volumes{m}, volumes{m} >= cutoffFraction * doseMax, 1, center, ...
                        [0 doseMax], jet(256));
                    colorbar(ax);
                    title(ax, [volNames{m} ' (Gy)'], 'Interpreter', 'none');
                end
                for m = 1:3   % row 2: gamma maps over each reference's cutoff region
                    ax = subplot(3, 3, 3 + m, 'Parent', fig);
                    show_ct(ax, hu, 1, center, spacing, ctWindowHu);
                    add_overlay(ax, gammas{m}.maps{1}, gammas{m}.eval_mask, 1, center, [0 2], gammaCmap);
                    colorbar(ax);
                    title(ax, sprintf('gamma %s vs %s: %.1f%% pass', shortNames{compRef(m)}, ...
                        shortNames{compTgt(m)}, passRates(c, m)));
                end
                ax = subplot(3, 3, 7:9, 'Parent', fig);   % row 3: gamma histograms
                hold(ax, 'on');
                for m = 1:3
                    gammaVals = min(gammas{m}.maps{1}(gammas{m}.eval_mask > 0), 3);
                    histogram(ax, gammaVals, 0:0.05:3, 'Normalization', 'probability', ...
                        'DisplayStyle', 'stairs', 'LineWidth', 1.5, 'DisplayName', sprintf( ...
                        '%s vs %s (%.1f%% pass)', shortNames{compRef(m)}, shortNames{compTgt(m)}, passRates(c, m)));
                end
                xline(ax, 1, 'k--', 'HandleVisibility', 'off');
                hold(ax, 'off');
                legend(ax);
                xlabel(ax, 'gamma index (clipped at 3); <= 1 passes');
                ylabel(ax, 'fraction of evaluated voxels');
                nCt = nnz(strcmp(ctLabel, ctLabels{c}));
                sgtitle(fig, sprintf('%s %s | %s totals (%d of %d recons valid) | axial slice %d (ETHOS max)', ...
                    patient_id, session, ctLabels{c}, nnz(reconOk & strcmp(ctLabel, ctLabels{c})), ...
                    nCt, center(3)), 'Interpreter', 'none');
                save_figure(fig, outDir, sprintf('totals_%s_%s', ctLabels{c}, config_hash), showFigures);
            end
            tests = add_result(tests, testName, 'INFO', sprintf( ...
                ['gamma pass, ETHOS-RS / RS-recon / ETHOS-recon: CT_1 %.1f / %.1f / %.1f%%, ' ...
                 'CT_3 %.1f / %.1f / %.1f%%; figures totals_CT_*.png'], passRates(1, :), passRates(2, :)));
        end
    catch ME
        tests = add_result(tests, testName, 'ERROR', ME.message);
    end

    %% ---------------- TEST 6: sensor signal after each processing stage ----------------
    testName = '6. Sensor signal after each processing stage';
    fprintf('\n[TEST 6] %s\n', testName);
    try
        if ~exist('kWaveGrid', 'file')
            tests = add_result(tests, testName, 'SKIP', 'k-Wave is not on the MATLAB path');
        else
            % Re-run the forward simulation of one CT_1 field (a test-4 pick when there is one)
            stageField = find(strcmp(ctLabel, 'CT_1') & reconOk, 1);
            if ~isempty(picks)
                stageField = pairCt1(picks(1));
            end
            fd     = load_field_dose_file(fieldIndex(stageField).file);
            medium = create_acoustic_medium(cbcts{1}, config);
            precomputedSensor = [];
            if isfile(sensorFile)
                loaded = load(sensorFile, 'precomputed_sensor');
                precomputedSensor = loaded.precomputed_sensor;
            end
            beamMetadata = [];
            if isfield(metadata, 'beam_metadata')
                beamMetadata = metadata.beam_metadata;
            end
            simConfig = config;
            simConfig.sensor_stages_only = true;    % stop before the reconstruction
            simConfig.plot_results       = false;
            [~, simResults] = run_single_field_simulation(fd, cbcts{1}, medium, ...
                beamMetadata, simConfig, precomputedSensor);

            if ~isfield(simResults, 'sensor_stages')
                tests = add_result(tests, testName, 'FAIL', sprintf( ...
                    'no sensor stages returned for %s (forward sim failed or no dose; see warnings above)', ...
                    fieldIndex(stageField).source_mat_filename));
            else
                stages = simResults.sensor_stages;
                timeUs = (0:size(stages.mean_abs, 2) - 1) * stages.dt * 1e6;
                fig = new_figure('Sensor signal per stage', showFigures);
                ax  = axes('Parent', fig);
                semilogy(ax, timeUs, stages.mean_abs', 'LineWidth', 1.2);
                grid(ax, 'on');
                legend(ax, stages.names, 'Location', 'best');
                xlabel(ax, 'time (\mus)');
                ylabel(ax, 'mean |p| over all sensor points (Pa)');
                snrText = '';
                if isfield(simResults, 'snr')
                    snrText = sprintf(' | SNR %.2f', simResults.snr);
                end
                title(ax, sprintf('%s %s | %s%s', patient_id, session, ...
                    fieldIndex(stageField).source_mat_filename, snrText), 'Interpreter', 'none');
                save_figure(fig, outDir, ['sensor_stages_' config_hash], showFigures);

                peakText = '';
                for s = 1:numel(stages.names)
                    peakText = sprintf('%s%s %.2e, ', peakText, stages.names{s}, max(stages.mean_abs(s, :)));
                end
                tests = add_result(tests, testName, 'INFO', sprintf( ...
                    'field %s; peak mean|p| (Pa): %sfigure sensor_stages_%s.png', ...
                    fieldIndex(stageField).source_mat_filename, peakText, config_hash));
            end
        end
    catch ME
        tests = add_result(tests, testName, 'ERROR', ME.message);
    end

    %% ---------------- TEST 7: CBCT1 and CBCT3 are different images ----------------
    testName = '7. CBCT1 and CBCT3 are different images';
    fprintf('\n[TEST 7] %s\n', testName);
    try
        huDiff      = double(cbcts{2}.cubeHU) - double(cbcts{1}.cubeHU);
        bodyUnion   = cbcts{1}.bodyMask | cbcts{2}.bodyMask;
        meanAbsDiff = mean(abs(huDiff(bodyUnion)));
        sameUid     = isfield(cbcts{1}, 'series_uid') && isfield(cbcts{2}, 'series_uid') ...
            && strcmp(char(cbcts{1}.series_uid), char(cbcts{2}.series_uid));

        [rows, cols, slices] = ind2sub(size(bodyUnion), find(cbcts{1}.bodyMask));
        center = round([mean(rows), mean(cols), mean(slices)]);
        fig = new_figure('CBCT1 vs CBCT3', showFigures);
        for v = 1:3
            for c = 1:2   % rows 1-2: each CBCT with its body contour
                ax = subplot(3, 3, (c - 1) * 3 + v, 'Parent', fig);
                show_ct(ax, cbcts{c}.cubeHU, v, center, spacing, ctWindowHu);
                add_contour(ax, cbcts{c}.bodyMask, v, center, 'g');
                title(ax, sprintf('%s %s', cbctNames{c}, viewNames{v}));
            end
            ax = subplot(3, 3, 6 + v, 'Parent', fig);   % row 3: difference inside the body
            show_ct(ax, cbcts{1}.cubeHU, v, center, spacing, ctWindowHu);
            add_overlay(ax, huDiff, bodyUnion, v, center, [-huDiffWindow huDiffWindow], divergingCmap);
            colorbar(ax);
            title(ax, sprintf('CBCT3 - CBCT1 (HU) %s', viewNames{v}));
        end
        sgtitle(fig, sprintf('%s %s | slices through the CBCT1 body centroid | mean |HU3 - HU1| in body = %.1f HU', ...
            patient_id, session, meanAbsDiff), 'Interpreter', 'none');
        save_figure(fig, outDir, ['cbct_difference_' config_hash], showFigures);

        tests = add_result(tests, testName, pass_fail(~sameUid && meanAbsDiff > minCbctHuDiff), sprintf( ...
            'same SeriesInstanceUID = %d; mean |HU_CBCT3 - HU_CBCT1| inside body = %.1f HU; figure cbct_difference_%s.png', ...
            sameUid, meanAbsDiff, config_hash));
    catch ME
        tests = add_result(tests, testName, 'ERROR', ME.message);
    end

    %% ---------------- TEST 8: silent failures, per-field dose scale and SNR ----------------
    testName = '8. (added) No silent k-Wave failures; per-field dose scale and SNR';
    fprintf('\n[TEST 8] %s\n', testName);
    try
        badRecon    = (hasRecon & ~reconOk) | loadFailed;   % zero, NaN/Inf, wrong size or unreadable
        okIdx       = find(reconOk);
        medianRatio = median(doseRatio(okIdx));
        isOutlier   = doseRatio(okIdx) > ratioOutlier * medianRatio | doseRatio(okIdx) < medianRatio / ratioOutlier;
        outliers    = okIdx(isOutlier);

        fig = new_figure('Per-field dose scale and SNR', showFigures);
        ax  = subplot(2, 1, 1, 'Parent', fig);
        hold(ax, 'on');
        for c = 1:2
            idx = find(strcmp(ctLabel, ctLabels{c}) & reconOk);
            plot(ax, idx, doseRatio(idx), '.', 'MarkerSize', 10, 'DisplayName', ctLabels{c});
        end
        plot(ax, find(badRecon), zeros(nnz(badRecon), 1), 'kx', 'MarkerSize', 8, ...
            'DisplayName', 'zero / invalid recon');
        if ~isnan(medianRatio)
            yline(ax, medianRatio, 'k--', 'DisplayName', 'median');
        end
        hold(ax, 'off');
        lgd = legend(ax);
        lgd.Interpreter = 'none';
        ylabel(ax, 'sum(recon) / sum(truth)');
        title(ax, sprintf('Recon / RS truth dose inside the truth''s %g%% region (raw recon, no gain)', ...
            100 * cutoffFraction));

        ax = subplot(2, 1, 2, 'Parent', fig);
        hold(ax, 'on');
        for c = 1:2
            idx = find(strcmp(ctLabel, ctLabels{c}) & reconOk);
            plot(ax, idx, fieldSnr(idx), '.', 'MarkerSize', 10, 'DisplayName', ctLabels{c});
        end
        hold(ax, 'off');
        set(ax, 'YScale', 'log');
        lgd = legend(ax);
        lgd.Interpreter = 'none';
        xlabel(ax, 'field index (sorted by beam, then segment)');
        ylabel(ax, 'SNR (signal peak / noise amp)');
        title(ax, 'SNR each field saw (noise\_stats saved with the recon)');
        sgtitle(fig, sprintf('%s %s | %d recons, %d zero/invalid', patient_id, session, ...
            nnz(hasRecon), nnz(badRecon)), 'Interpreter', 'none');
        save_figure(fig, outDir, ['field_scale_snr_' config_hash], showFigures);

        tests = add_result(tests, testName, pass_fail(~any(badRecon)), sprintf( ...
            ['%d of %d recons zero / NaN / wrong size / unreadable; recon/truth ratio median %.3g, ' ...
             '%d field(s) beyond x%g of it; median SNR %.2f (%d fields have SNR); figure field_scale_snr_%s.png'], ...
            nnz(badRecon), nnz(hasRecon), medianRatio, numel(outliers), ratioOutlier, ...
            median(fieldSnr, 'omitnan'), nnz(~isnan(fieldSnr)), config_hash));
        print_examples('zero/invalid recon', {fieldIndex(badRecon).source_mat_filename});
        print_examples('ratio outlier', {fieldIndex(outliers).source_mat_filename});
    catch ME
        tests = add_result(tests, testName, 'ERROR', ME.message);
    end

    %% ========================= SUMMARY =======================================
    statuses = {tests.status};
    fprintf('\n=========================================================\n');
    fprintf('  Verification %s / %s: %d PASS, %d FAIL, %d ERROR, %d INFO, %d SKIP\n', ...
        patient_id, session, sum(strcmp(statuses, 'PASS')), sum(strcmp(statuses, 'FAIL')), ...
        sum(strcmp(statuses, 'ERROR')), sum(strcmp(statuses, 'INFO')), sum(strcmp(statuses, 'SKIP')));
    for j = 1:numel(tests)
        if any(strcmp(tests(j).status, {'FAIL', 'ERROR'}))
            fprintf('  [%s] %s - %s\n', tests(j).status, tests(j).name, tests(j).detail);
        end
    end
    fprintf('  Figures: %s\n', outDir);
    fprintf('=========================================================\n');

    results = struct('patient_id', patient_id, 'session', session, ...
        'config_hash', config_hash, 'out_dir', outDir, 'tests', tests);
end


%% =========================================================================
%  HELPER FUNCTIONS
%% =========================================================================

function tests = add_result(tests, name, status, detail)
%ADD_RESULT Append one test result to the list and print it.
    tests(end + 1) = struct('name', name, 'status', status, 'detail', detail);
    fprintf('  [%s] %s - %s\n', status, name, detail);
end

function status = pass_fail(passed)
%PASS_FAIL 'PASS' when passed is true, otherwise 'FAIL'.
    if passed
        status = 'PASS';
    else
        status = 'FAIL';
    end
end

function print_examples(label, names)
%PRINT_EXAMPLES Print up to 5 offending file names / keys under a result line.
    for i = 1:min(5, numel(names))
        fprintf('      %s: %s\n', label, names{i});
    end
    if numel(names) > 5
        fprintf('      ... and %d more\n', numel(names) - 5);
    end
end

function value = config_value(config, name, default)
%CONFIG_VALUE config.(name) when it is set, otherwise default.
    if isfield(config, name) && ~isempty(config.(name))
        value = config.(name);
    else
        value = default;
    end
end

function filePath = find_dicom(fileName, folders)
%FIND_DICOM Full path of fileName in the first folder that has it.
    for k = 1:numel(folders)
        filePath = fullfile(folders{k}, fileName);
        if isfile(filePath)
            return;
        end
    end
    error('verify_pipeline_simulate:FileNotFound', '%s not found in: %s', ...
        fileName, strjoin(folders, ', '));
end

function out = embed_in_grid(vol, gridSize, offset, fillValue)
%EMBED_IN_GRID Put vol inside a larger gridSize array after offset = [rows cols
%   slices] of padding; every other voxel is fillValue.
    out = repmat(cast(fillValue, 'like', vol), gridSize);
    out(offset(1) + (1:size(vol, 1)), offset(2) + (1:size(vol, 2)), ...
        offset(3) + (1:size(vol, 3))) = vol;
end

function fig = new_figure(name, showFigures)
%NEW_FIGURE Blank figure; hidden when MATLAB has no desktop to show it on.
    visible = 'off';
    if showFigures
        visible = 'on';
    end
    fig = figure('Name', name, 'Color', 'w', 'Visible', visible, 'Position', [50 50 1500 900]);
end

function save_figure(fig, outDir, fileStem, showFigures)
%SAVE_FIGURE Save fig as <outDir>/<fileStem>.png; close it when it is not shown.
    drawnow;
    exportgraphics(fig, fullfile(outDir, [fileStem '.png']), 'Resolution', 150);
    if ~showFigures
        close(fig);
    end
end

function show_ct(ax, hu, viewNum, center, spacing, ctWindow)
%SHOW_CT Grayscale CT slice through center = [row col slice] in true proportions.
%   viewNum 1 = axial (anterior at top, patient right on the left), 2 = coronal,
%   3 = sagittal (superior at top; sagittal has anterior on the left).
%   Leaves hold on so overlays and contours can be drawn on top.
    ctSlice = get_slice(hu, viewNum, center);
    ctGray  = min(max((ctSlice - ctWindow(1)) / diff(ctWindow), 0), 1);
    image(ax, repmat(ctGray, 1, 1, 3));   % truecolor, so the colormap only colors overlays
    hold(ax, 'on');
    % Voxel size (mm) along each view's [horizontal, vertical] screen axis
    viewSpacing = [spacing(1) spacing(2); spacing(1) spacing(3); spacing(2) spacing(3)];
    daspect(ax, [viewSpacing(viewNum, 2) viewSpacing(viewNum, 1) 1]);
    if viewNum > 1
        set(ax, 'YDir', 'normal');   % higher slice index (superior) at the top
    end
    axis(ax, 'off');
end

function add_overlay(ax, vol, showMask, viewNum, center, range, cmap)
%ADD_OVERLAY Color-mapped slice of vol over the CT: opaque where showMask is
%   true, transparent elsewhere, colors fixed to range = [low high].
    h = imagesc(ax, get_slice(vol, viewNum, center), range);
    set(h, 'AlphaData', 0.7 * get_slice(showMask, viewNum, center));
    colormap(ax, cmap);
    if viewNum > 1
        set(ax, 'YDir', 'normal');   % imagesc resets YDir; keep superior at the top
    end
end

function add_contour(ax, mask, viewNum, center, color)
%ADD_CONTOUR Outline of a 3D mask on the slice (skipped when empty there).
    maskSlice = get_slice(mask, viewNum, center);
    if any(maskSlice(:))
        contour(ax, maskSlice, [0.5 0.5], 'LineColor', color, 'LineWidth', 1);
    end
end

function s = get_slice(vol, viewNum, center)
%GET_SLICE 2D slice of a (rows = y, cols = x, slices = z) volume through center.
%   1 = axial (y down, x across), 2 = coronal (z up, x across),
%   3 = sagittal (z up, y across).
    switch viewNum
        case 1
            s = vol(:, :, center(3));
        case 2
            s = squeeze(vol(center(1), :, :))';
        case 3
            s = squeeze(vol(:, center(2), :))';
    end
    s = double(s);
end
