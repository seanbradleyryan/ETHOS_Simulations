%% =========================================================================
%  VERIFY_PIPELINE_COMPRESS_OUTPUT.m
%  ETHOS Photoacoustic Pipeline - Checks on the Step 1.4 / 1.5 outputs
%  =========================================================================
%
%  PURPOSE:
%  Run after pipeline_compress.m and before pipeline_simulate.m. Confirms that
%  RayStationFiles/<id>/<session>/processed/ is complete and physically sane,
%  and draws figures for a visual check.
%
%  TESTS (per patient/session):
%    1  10 random beam/segment pairs: the CT_1 dose on CBCT1 and the CT_3 dose
%       on CBCT3 (orthogonal views through the CT_1 max dose), and a 2%/2mm
%       gamma of CT_1 vs CT_3.                                 [figures, INFO]
%    2  Body contour drawn on each CBCT.                        [figure, INFO]
%    3  Every dose_*.npz has its Step 1.4 .mat and its Step 1.5 processed .mat
%       (and every processed dose traces back to an NPZ).
%    4  2%/2mm gamma of the summed RS dose on CT_1 vs the ETHOS RTDOSE truth.
%    5  CBCT1 and CBCT3 are different images.
%    6  The CT_1 and CT_3 doses of every beam/segment are different.
%    7  No empty doses or CBCTs.
%    8a Every field dose carries its RTPLAN beam data (gantry, MU, isocenter,
%       jaws). pipeline_simulate uses these for the gantry angle and sensor
%       placement; Step 1.5 silently writes 0 / [] when the beam is not found.
%    8b The saved totals equal the sum of the per-field doses, and every dose
%       and CBCT is on the metadata.mat grid (catches stale files left behind
%       by a resumed run with skip_completed = true).
%
%  Each test prints PASS, FAIL, INFO (figure / numbers to review), SKIP, or
%  ERROR (the test itself crashed). FAIL/ERROR lines are repeated at the end.
%
%  WHY TEST 4 USES total_dose_CT_1.mat AND NOT total_rs_dose.mat:
%  total_rs_dose is the sum of the CT_1 AND CT_3 field doses, i.e. the same
%  plan delivered twice, so it is ~2x the ETHOS dose. The CT_1 total is one
%  delivery of the plan, which is what the ETHOS RTDOSE describes. RS doses
%  are for ONE fraction (Step 0.6 divides the MU by the fraction count), so a
%  DoseSummationType = PLAN RTDOSE is divided by NumberOfFractionsPlanned.
%
%  INPUTS (set in the CONFIGURATION section):
%    CONFIG.patients / .sessions / .treatment_site - what to verify
%    CONFIG.num_random_pairs, .random_seed          - test 1 sampling
%    CONFIG.gamma_*, .ethos_min_pass_pct            - gamma (tests 1 and 4)
%    CONFIG.*_tol, .min_*                           - pass/fail thresholds
%    CONFIG.ct_window_hu, .close_figures            - plotting
%
%  OUTPUTS:
%    Console report, and PNG figures in
%    AnalysisResults/<id>/<session>/compress_verification/
%      pair_<plan>_B<beam>_S<seg>.png   (test 1)
%      body_contours.png                (test 2)
%
%  ALGORITHM:
%    1. Load metadata.mat and both CBCT*_resampled.mat; index the processed
%       dose files and pair each CT_1 file with its CT_3 file.
%    2. Tests 1-5, each in its own try/catch.
%    3. One parfor pass loads every field dose once for tests 6, 7, 8a, 8b and
%       sums the doses for comparison with the saved totals.
%
%  EXAMPLE:
%    Set CONFIG.patients / CONFIG.sessions below, then run the script.
%
%  DEPENDENCIES:
%    utils/ (get_repo_root, load_field_dose_file, compute_gamma, CalcGamma),
%    Image Processing Toolbox (dicominfo, dicomread),
%    Parallel Computing Toolbox (optional; parfor runs serially without it).
%
%  See also: pipeline_compress, step14_npz_to_mat, step15_process_doses,
%            pipeline_simulate, compute_gamma
%  =========================================================================

clear; clc; close all;

%% ========================= CONFIGURATION =================================

% --- Patient and Session Selection (copy the lists from pipeline_compress.m) ---
CONFIG.patients       = {'1194203'};
CONFIG.sessions       = {'Session_1'};
CONFIG.treatment_site = 'Pancreas';

% --- Directory Paths ---
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'utils')));
CONFIG.working_dir = get_repo_root();

% --- Test 1: random beam/segment pairs ---
CONFIG.num_random_pairs = 10;
CONFIG.random_seed      = 'shuffle';   % set to an integer to redraw the same pairs

% --- Gamma (tests 1 and 4): global gamma ---
CONFIG.gamma_dose_pct        = 2;      % dose difference, % of the reference max
CONFIG.gamma_dist_mm         = 2;      % distance to agreement, mm
CONFIG.gamma_cutoff_fraction = 0.10;   % ignore voxels below 10% of the reference max
CONFIG.ethos_min_pass_pct    = 90;     % test 4 passes at or above this pass rate

% --- Pass/fail thresholds ---
CONFIG.identical_rel_tol  = 1e-6;   % test 6: CT_1 vs CT_3 relative difference below this = identical
CONFIG.min_cbct_hu_diff   = 1;      % test 5: mean |HU1 - HU3| in body below this = same image
CONFIG.min_body_median_hu = -500;   % test 7: body median HU below this = empty or misaligned CBCT
CONFIG.total_rel_tol      = 1e-6;   % test 8b: saved total vs summed field doses
CONFIG.geometry_tol_mm    = 0.01;   % test 8b: CBCT origin/spacing vs metadata.mat

% --- Plotting ---
CONFIG.ct_window_hu  = [-500 500];
% Keep figures open for a single session; close them (PNG still saved) for many
CONFIG.close_figures = numel(CONFIG.patients) * numel(CONFIG.sessions) > 1;

%% ========================= INITIALIZATION ================================

fprintf('=========================================================\n');
fprintf('  ETHOS Pipeline - Verify pipeline_compress output\n');
fprintf('=========================================================\n');
fprintf('  Started: %s\n', datetime('now'));
fprintf('  Working directory: %s\n', CONFIG.working_dir);
fprintf('=========================================================\n');

rng(CONFIG.random_seed);
gammaCriteria = {CONFIG.gamma_dose_pct, CONFIG.gamma_dist_mm, ...
    sprintf('%g%%/%gmm', CONFIG.gamma_dose_pct, CONFIG.gamma_dist_mm)};
allResults = struct('patient_id', {}, 'session', {}, 'tests', {});

%% ========================= MAIN LOOP =====================================

for pIdx = 1:numel(CONFIG.patients)
    patientID = CONFIG.patients{pIdx};

    for sIdx = 1:numel(CONFIG.sessions)
        session = CONFIG.sessions{sIdx};
        fprintf('\n=== Verifying: Patient %s, %s ===\n', patientID, session);
        tests = struct('name', {}, 'status', {}, 'detail', {});

        try
            %% ---------------- Paths and shared inputs ----------------
            rsDir        = fullfile(CONFIG.working_dir, 'RayStationFiles', patientID, session);
            processedDir = fullfile(rsDir, 'processed');
            sctDir       = fullfile(CONFIG.working_dir, 'EthosExports', patientID, ...
                CONFIG.treatment_site, session, 'sct');
            outDir       = fullfile(CONFIG.working_dir, 'AnalysisResults', patientID, session, ...
                'compress_verification');
            if ~isfolder(outDir)
                mkdir(outDir);
            end

            loaded   = load(fullfile(processedDir, 'metadata.mat'), 'metadata');
            metadata = loaded.metadata;
            loaded   = load(fullfile(processedDir, 'CBCT1_resampled.mat'), 'CBCT1_resampled');
            cbct1    = loaded.CBCT1_resampled;
            loaded   = load(fullfile(processedDir, 'CBCT3_resampled.mat'), 'CBCT3_resampled');
            cbct3    = loaded.CBCT3_resampled;
            clear loaded;

            gridDims = metadata.dimensions(:)';   % [ny nx nz] = (rows, cols, slices)
            origin   = metadata.origin(:)';       % [x y z] mm of voxel (1,1,1)
            spacing  = metadata.spacing(:)';      % [dx dy dz] mm
            % CalcGamma voxel widths follow the array dimensions: rows = y, cols = x
            gammaWidth = [spacing(2) spacing(1) spacing(3)];

            %% ---------------- Pair every CT_1 dose with its CT_3 dose ----------------
            % Processed names: dose_<id>_<session>_<plan>_<CT_n>_B<beam>_<seg>.mat
            % pairKeys{k} = '<plan>_B<beam>_S<seg>'; pairCt1{k} / pairCt3{k} = file
            % path, or '' when that CT's file is missing. Files without a plan type
            % and CT_1/CT_3 label are skipped here and reported by test 3.
            doseFiles = dir(fullfile(processedDir, 'dose_*.mat'));
            pairKeys  = {};
            pairPlan  = {};
            pairCt1   = {};
            pairCt3   = {};
            for i = 1:numel(doseFiles)
                name = doseFiles(i).name;
                tok  = regexp(name, '_(adapted|reference)_(CT_\d+)_B(\d+)_(\d+)\.mat$', ...
                    'tokens', 'once', 'ignorecase');
                if isempty(tok) || ~any(strcmpi(tok{2}, {'CT_1', 'CT_3'}))
                    continue;
                end
                key = sprintf('%s_B%s_S%s', lower(tok{1}), tok{3}, tok{4});
                k = find(strcmp(pairKeys, key), 1);
                if isempty(k)
                    pairKeys{end + 1} = key;
                    pairPlan{end + 1} = lower(tok{1});
                    pairCt1{end + 1}  = '';
                    pairCt3{end + 1}  = '';
                    k = numel(pairKeys);
                end
                if strcmpi(tok{2}, 'CT_1')
                    pairCt1{k} = fullfile(processedDir, name);
                else
                    pairCt3{k} = fullfile(processedDir, name);
                end
            end
            fprintf('  %d processed dose files -> %d beam/segment keys\n', ...
                numel(doseFiles), numel(pairKeys));

            %% ---------------- TEST 1: random pairs, dose on each CBCT + gamma ----------------
            testName = '1. Random pairs: dose on each CBCT + CT_1 vs CT_3 gamma';
            fprintf('\n[TEST 1] %s\n', testName);
            try
                complete  = find(~cellfun(@isempty, pairCt1) & ~cellfun(@isempty, pairCt3));
                nPick     = min(CONFIG.num_random_pairs, numel(complete));
                picks     = complete(randperm(numel(complete), nPick));
                passRates = nan(1, nPick);
                for n = 1:nPick
                    k   = picks(n);
                    fd1 = load_field_dose_file(pairCt1{k});
                    fd3 = load_field_dose_file(pairCt3{k});
                    doseMax = max(max(fd1.dose_Gy(:)), max(fd3.dose_Gy(:)));
                    [~, iMax] = max(fd1.dose_Gy(:));
                    [row, col, slc] = ind2sub(size(fd1.dose_Gy), iMax);

                    g = compute_gamma(fd1.dose_Gy, fd3.dose_Gy, gammaWidth, ...
                        'Criteria', gammaCriteria, 'Cutoff', CONFIG.gamma_cutoff_fraction);
                    passRates(n) = g.pass_rates(1);

                    fig = figure('Name', pairKeys{k}, 'Position', [50 50 1500 850]);
                    plot_three_views(1, cbct1.cubeHU, cbct1.bodyMask, fd1.dose_Gy, ...
                        [row col slc], spacing, doseMax, CONFIG, 'CT_1 dose on CBCT1');
                    plot_three_views(2, cbct3.cubeHU, cbct3.bodyMask, fd3.dose_Gy, ...
                        [row col slc], spacing, doseMax, CONFIG, 'CT_3 dose on CBCT3');
                    sgtitle(fig, sprintf('%s %s %s | slices through CT_1 max | CT_1 vs CT_3 gamma %s: %.1f%% pass', ...
                        patientID, session, pairKeys{k}, gammaCriteria{3}, passRates(n)), ...
                        'Interpreter', 'none');
                    exportgraphics(fig, fullfile(outDir, sprintf('pair_%s.png', pairKeys{k})), ...
                        'Resolution', 150);
                    if CONFIG.close_figures
                        close(fig);
                    end
                    fprintf('    %s: gamma %s pass %.1f%%\n', pairKeys{k}, gammaCriteria{3}, passRates(n));
                end

                if nPick == 0
                    tests = add_result(tests, testName, 'FAIL', ...
                        'no beam/segment has both a CT_1 and a CT_3 dose');
                else
                    tests = add_result(tests, testName, 'INFO', sprintf( ...
                        '%d pairs, gamma %s pass: mean %.1f%% (min %.1f%%, max %.1f%%); figures in %s', ...
                        nPick, gammaCriteria{3}, mean(passRates), min(passRates), max(passRates), outDir));
                end
            catch ME
                tests = add_result(tests, testName, 'ERROR', ME.message);
            end

            %% ---------------- TEST 2: body contour on each CBCT ----------------
            testName = '2. Body contour on each CBCT';
            fprintf('\n[TEST 2] %s\n', testName);
            try
                [row, col, slc] = ind2sub(size(cbct1.bodyMask), find(cbct1.bodyMask));
                bodyCenter = round([mean(row) mean(col) mean(slc)]);

                fig = figure('Name', 'Body contours', 'Position', [50 50 1500 850]);
                plot_three_views(1, cbct1.cubeHU, cbct1.bodyMask, [], bodyCenter, spacing, 0, ...
                    CONFIG, 'CBCT1 (CT_1)');
                plot_three_views(2, cbct3.cubeHU, cbct3.bodyMask, [], bodyCenter, spacing, 0, ...
                    CONFIG, 'CBCT3 (CT_3)');
                sgtitle(fig, sprintf('%s %s | body contour (green), slices through the CBCT1 body centroid', ...
                    patientID, session), 'Interpreter', 'none');
                exportgraphics(fig, fullfile(outDir, 'body_contours.png'), 'Resolution', 150);
                if CONFIG.close_figures
                    close(fig);
                end

                tests = add_result(tests, testName, 'INFO', sprintf( ...
                    'body voxels: CBCT1 %d, CBCT3 %d; figure body_contours.png', ...
                    nnz(cbct1.bodyMask), nnz(cbct3.bodyMask)));
            catch ME
                tests = add_result(tests, testName, 'ERROR', ME.message);
            end

            %% ---------------- TEST 3: every NPZ became a .mat ----------------
            testName = '3. Every NPZ converted to .mat (Step 1.4) and processed (Step 1.5)';
            fprintf('\n[TEST 3] %s\n', testName);
            try
                npzFiles = dir(fullfile(rsDir, 'dose_*.npz'));
                if isempty(npzFiles)
                    tests = add_result(tests, testName, 'SKIP', ...
                        'no dose_*.npz in the RayStation folder (legacy DICOM input)');
                else
                    missingRaw       = {};
                    missingProcessed = {};
                    expectedNames    = repmat({''}, 1, numel(npzFiles));
                    for i = 1:numel(npzFiles)
                        name = npzFiles(i).name;
                        [~, stem] = fileparts(name);

                        % Step 1.4 keeps the stem: dose_...npz -> dose_...mat in the same folder
                        if ~isfile(fullfile(rsDir, [stem '.mat']))
                            missingRaw{end + 1} = name;
                        end

                        % Step 1.5 renames: dose_<id>_<session>_<plan>_<CT_n>_B<beam>_<seg>.mat,
                        % with the segment written without zero padding
                        tok = regexp(name, '_(adapted|reference)_(CT_\d+)_B(\d+)_(\d+)\.npz$', ...
                            'tokens', 'once', 'ignorecase');
                        if isempty(tok)
                            missingProcessed{end + 1} = name;   % Step 1.5 cannot parse this name
                            continue;
                        end
                        expectedNames{i} = sprintf('dose_%s_%s_%s_%s_B%d_%d.mat', patientID, session, ...
                            lower(tok{1}), tok{2}, str2double(tok{3}), str2double(tok{4}));
                        if ~isfile(fullfile(processedDir, expectedNames{i}))
                            missingProcessed{end + 1} = name;
                        end
                    end
                    % Processed doses that no current NPZ produced (stale files pipeline_simulate would still run)
                    orphans = setdiff({doseFiles.name}, expectedNames);

                    tests = add_result(tests, testName, ...
                        pass_fail(isempty(missingRaw) && isempty(missingProcessed) && isempty(orphans)), ...
                        sprintf('%d NPZ: %d without a Step 1.4 .mat, %d without a processed .mat; %d processed file(s) with no NPZ', ...
                        numel(npzFiles), numel(missingRaw), numel(missingProcessed), numel(orphans)));
                    print_examples('no Step 1.4 .mat', missingRaw);
                    print_examples('no processed .mat', missingProcessed);
                    print_examples('no NPZ', orphans);
                end
            catch ME
                tests = add_result(tests, testName, 'ERROR', ME.message);
            end

            %% ---------------- TEST 4: RS dose on CT_1 vs ETHOS truth ----------------
            testName = '4. RS dose on CT_1 vs ETHOS truth (gamma)';
            fprintf('\n[TEST 4] %s\n', testName);
            try
                planTypes = unique(pairPlan);
                if numel(planTypes) ~= 1
                    tests = add_result(tests, testName, 'SKIP', sprintf( ...
                        'expected one plan type in processed/, found %d (%s); total_dose_CT_1 would mix plans', ...
                        numel(planTypes), strjoin(planTypes, ', ')));
                else
                    rsCt1Total = load_total_dose(fullfile(processedDir, 'total_dose_CT_1.mat'));

                    % ETHOS truth for the same plan type, in Gy
                    rtdoseFile = fullfile(sctDir, sprintf('RTDOSE_%s.dcm', planTypes{1}));
                    info       = dicominfo(rtdoseFile);
                    ethosRaw   = double(squeeze(dicomread(rtdoseFile))) * info.DoseGridScaling;

                    % RS doses are one fraction; a PLAN-summed RTDOSE covers every fraction
                    nFractions = 1;
                    if strcmpi(info.DoseSummationType, 'PLAN')
                        rtplan     = dicominfo(fullfile(sctDir, sprintf('RTPLAN_%s.dcm', planTypes{1})));
                        nFractions = double(rtplan.FractionGroupSequence.Item_1.NumberOfFractionsPlanned);
                    end
                    ethosRaw = ethosRaw / nFractions;

                    % Resample ETHOS onto the RS dose grid by patient position (mm), not
                    % by array size. DICOM PixelSpacing = [row (y), column (x)] spacing.
                    ethosX = info.ImagePositionPatient(1) + (0:size(ethosRaw, 2) - 1) * info.PixelSpacing(2);
                    ethosY = info.ImagePositionPatient(2) + (0:size(ethosRaw, 1) - 1) * info.PixelSpacing(1);
                    ethosZ = info.GridFrameOffsetVector(:)';
                    if ethosZ(1) == 0   % offsets relative to the first frame (the usual case)
                        ethosZ = ethosZ + info.ImagePositionPatient(3);
                    end
                    [qx, qy, qz] = meshgrid(origin(1) + (0:gridDims(2) - 1) * spacing(1), ...
                                            origin(2) + (0:gridDims(1) - 1) * spacing(2), ...
                                            origin(3) + (0:gridDims(3) - 1) * spacing(3));
                    ethos = interp3(ethosX, ethosY, ethosZ, ethosRaw, qx, qy, qz, 'linear', 0);
                    clear qx qy qz ethosRaw;

                    % Zero ETHOS outside the CT_1 body / inside the couch, as Step 1.5 did to the RS doses
                    ethos(~(cbct1.bodyMask & ~cbct1.couchMask)) = 0;

                    g = compute_gamma(ethos, rsCt1Total, gammaWidth, ...
                        'Criteria', gammaCriteria, 'Cutoff', CONFIG.gamma_cutoff_fraction);
                    % A centroid shift points at a grid-origin problem; a max-dose ratio at a scaling problem
                    shiftMm = dose_centroid_mm(rsCt1Total, origin, spacing) - dose_centroid_mm(ethos, origin, spacing);

                    tests = add_result(tests, testName, ...
                        pass_fail(g.pass_rates(1) >= CONFIG.ethos_min_pass_pct), sprintf( ...
                        ['%s pass %.1f%% (threshold %.0f%%); max dose ETHOS %.3f Gy/fx (%s, %d fx), ' ...
                         'RS %.3f Gy; RS - ETHOS dose centroid [%.1f %.1f %.1f] mm'], ...
                        gammaCriteria{3}, g.pass_rates(1), CONFIG.ethos_min_pass_pct, max(ethos(:)), ...
                        info.DoseSummationType, nFractions, max(rsCt1Total(:)), shiftMm));
                    clear ethos rsCt1Total;
                end
            catch ME
                tests = add_result(tests, testName, 'ERROR', ME.message);
            end

            %% ---------------- TEST 5: CBCT1 and CBCT3 are different images ----------------
            testName = '5. CBCT1 and CBCT3 are different images';
            fprintf('\n[TEST 5] %s\n', testName);
            try
                sameUid    = strcmp(char(cbct1.series_uid), char(cbct3.series_uid));
                bodyUnion  = cbct1.bodyMask | cbct3.bodyMask;
                meanHuDiff = mean(abs(cbct1.cubeHU(bodyUnion) - cbct3.cubeHU(bodyUnion)));
                tests = add_result(tests, testName, ...
                    pass_fail(~sameUid && meanHuDiff > CONFIG.min_cbct_hu_diff), sprintf( ...
                    'same SeriesInstanceUID = %d; mean |HU_CBCT1 - HU_CBCT3| inside body = %.1f HU', ...
                    sameUid, meanHuDiff));
            catch ME
                tests = add_result(tests, testName, 'ERROR', ME.message);
            end

            %% ---------------- TESTS 6, 7, 8a, 8b: one pass over every field dose ----------------
            % Each iteration loads the CT_1 and CT_3 file of one beam/segment, checks
            % them, and adds them to running sums (sumCt1 / sumCt3 are parfor
            % "reduction" variables: each worker keeps its own sum and MATLAB adds
            % them together at the end). Every file is loaded exactly once.
            fprintf('\n[TESTS 6-8] Loading every processed dose file once...\n');
            try
                nPairs = numel(pairKeys);
                if nPairs == 0
                    error('verify_pipeline_compress_output:NoDoses', ...
                        'No labeled dose_*.mat files in %s', processedDir);
                end

                rtplanBeams = [];
                if isfield(metadata, 'beam_metadata') && ~isempty(metadata.beam_metadata)
                    rtplanBeams = [metadata.beam_metadata.beam_number];
                end
                % Voxels kept by the Step 1.5 masking on BOTH CBCTs. An identical RS dose
                % on both CTs would still differ at the body edge, so compare only here.
                validBoth = cbct1.bodyMask & ~cbct1.couchMask & cbct3.bodyMask & ~cbct3.couchMask;

                % Per pair, column 1 = CT_1 file, column 2 = CT_3 file
                pairTemplate = struct('has', [false false], 'dose_ok', [false false], ...
                    'grid_ok', [false false], 'meta_ok', [false false], 'rel_diff', NaN);
                pairResults = repmat(pairTemplate, nPairs, 1);
                sumCt1 = zeros(gridDims);
                sumCt3 = zeros(gridDims);

                parfor k = 1:nPairs
                    r     = pairTemplate;
                    files = {pairCt1{k}, pairCt3{k}};
                    doses = {[], []};
                    for ct = 1:2
                        if isempty(files{ct})
                            continue;
                        end
                        fd = load_field_dose_file(files{ct});
                        d  = fd.dose_Gy;
                        r.has(ct)     = true;
                        r.dose_ok(ct) = ~isempty(d) && all(isfinite(d(:))) && min(d(:)) >= 0 && max(d(:)) > 0;
                        r.grid_ok(ct) = isequal(size(d), gridDims);
                        r.meta_ok(ct) = ismember(fd.beam_num, rtplanBeams) ...
                            && isscalar(fd.gantry_angle) && isfinite(fd.gantry_angle) ...
                            && isscalar(fd.meterset) && fd.meterset > 0 ...
                            && numel(fd.isocenter) == 3 && numel(fd.jaw_x) == 2 && numel(fd.jaw_y) == 2;
                        if r.grid_ok(ct)
                            doses{ct} = d;
                        end
                    end
                    if ~isempty(doses{1})
                        sumCt1 = sumCt1 + doses{1};
                    end
                    if ~isempty(doses{2})
                        sumCt3 = sumCt3 + doses{2};
                    end
                    if ~isempty(doses{1}) && ~isempty(doses{2})
                        a = doses{1}(validBoth);
                        b = doses{2}(validBoth);
                        r.rel_diff = norm(a - b) / norm(a);
                    end
                    pairResults(k) = r;
                end

                has      = vertcat(pairResults.has);       % nPairs x 2
                doseOk   = vertcat(pairResults.dose_ok);
                gridOk   = vertcat(pairResults.grid_ok);
                metaOk   = vertcat(pairResults.meta_ok);
                relDiff  = vertcat(pairResults.rel_diff);
                allFiles = [pairCt1(:), pairCt3(:)];       % same layout as has

                % ---- Test 6: the CT_1 and CT_3 dose of each beam/segment differ ----
                both      = has(:, 1) & has(:, 2);
                unpaired  = pairKeys(xor(has(:, 1), has(:, 2)));
                identical = pairKeys(both & relDiff < CONFIG.identical_rel_tol);
                tests = add_result(tests, '6. CT_1 and CT_3 doses differ', ...
                    pass_fail(isempty(identical) && isempty(unpaired)), sprintf( ...
                    '%d pairs: %d identical, %d missing their CT_1/CT_3 partner; median relative difference %.1f%%', ...
                    nnz(both), numel(identical), numel(unpaired), 100 * median(relDiff(both), 'omitnan')));
                print_examples('identical', identical);
                print_examples('no partner', unpaired);

                % ---- CBCT checks used by tests 7 and 8b ----
                cbcts      = {cbct1, cbct3};
                cbctOk     = [false false];   % test 7: not empty, body contour sits on tissue
                geomOk     = [false false];   % test 8b: on the metadata.mat grid
                bodyMedian = [NaN NaN];
                for c = 1:2
                    hu     = cbcts{c}.cubeHU;
                    body   = cbcts{c}.bodyMask;
                    sizeOk = isequal(size(hu), gridDims) && isequal(size(body), gridDims);
                    geomOk(c) = sizeOk ...
                        && max(abs(cbcts{c}.origin(:)' - origin))   < CONFIG.geometry_tol_mm ...
                        && max(abs(cbcts{c}.spacing(:)' - spacing)) < CONFIG.geometry_tol_mm;
                    if sizeOk && any(body(:))
                        bodyMedian(c) = median(hu(body));
                        cbctOk(c) = all(isfinite(hu(:))) && std(hu(:)) > 0 ...
                            && bodyMedian(c) > CONFIG.min_body_median_hu;
                    end
                end

                % ---- Test 7: no empty doses or CBCTs ----
                badDoses = allFiles(has & ~doseOk);
                tests = add_result(tests, '7. No empty doses or CBCTs', ...
                    pass_fail(isempty(badDoses) && all(cbctOk)), sprintf( ...
                    '%d of %d field doses empty or invalid (NaN/Inf/negative); CBCT1 ok=%d (body median %.0f HU), CBCT3 ok=%d (body median %.0f HU)', ...
                    numel(badDoses), nnz(has), cbctOk(1), bodyMedian(1), cbctOk(2), bodyMedian(2)));
                print_examples('empty/invalid', badDoses);

                % ---- Test 8a: RTPLAN beam data on every field (needed by pipeline_simulate) ----
                badMeta = allFiles(has & ~metaOk);
                tests = add_result(tests, '8a. RTPLAN beam data on every field dose', ...
                    pass_fail(isempty(badMeta)), sprintf( ...
                    '%d of %d field doses lack a matching RTPLAN beam (gantry, MU > 0, isocenter, jaws)', ...
                    numel(badMeta), nnz(has)));
                print_examples('no RTPLAN data', badMeta);

                % ---- Test 8b: saved totals = sum of fields; everything on one grid ----
                badGrid = allFiles(has & ~gridOk);
                errCt1  = rel_max_diff(sumCt1, load_total_dose(fullfile(processedDir, 'total_dose_CT_1.mat')));
                errCt3  = rel_max_diff(sumCt3, load_total_dose(fullfile(processedDir, 'total_dose_CT_3.mat')));
                errAll  = rel_max_diff(sumCt1 + sumCt3, load_total_dose(fullfile(processedDir, 'total_rs_dose.mat')));
                tests = add_result(tests, '8b. Totals match field doses; one shared grid', ...
                    pass_fail(isempty(badGrid) && all(geomOk) ...
                        && all([errCt1 errCt3 errAll] < CONFIG.total_rel_tol)), sprintf( ...
                    ['%d field dose(s) off the metadata grid; CBCT grid ok: CBCT1=%d CBCT3=%d; ' ...
                     'max|summed - saved| / max(saved): CT_1 %.1e, CT_3 %.1e, total_rs %.1e'], ...
                    numel(badGrid), geomOk(1), geomOk(2), errCt1, errCt3, errAll));
                print_examples('off grid', badGrid);
            catch ME
                tests = add_result(tests, '6-8. Full pass over field doses', 'ERROR', ME.message);
            end

        catch ME
            tests = add_result(tests, 'Session setup (metadata / CBCT / file listing)', 'ERROR', ME.message);
        end

        allResults(end + 1) = struct('patient_id', patientID, 'session', session, 'tests', tests);
    end  % session loop
end  % patient loop

%% ========================= SUMMARY =======================================

fprintf('\n=========================================================\n');
fprintf('  Verification Summary (FAIL / ERROR lines repeated)\n');
fprintf('=========================================================\n');
for i = 1:numel(allResults)
    statuses = {allResults(i).tests.status};
    fprintf('\n%s / %s: %d PASS, %d FAIL, %d ERROR, %d INFO, %d SKIP\n', ...
        allResults(i).patient_id, allResults(i).session, ...
        sum(strcmp(statuses, 'PASS')), sum(strcmp(statuses, 'FAIL')), ...
        sum(strcmp(statuses, 'ERROR')), sum(strcmp(statuses, 'INFO')), ...
        sum(strcmp(statuses, 'SKIP')));
    for j = 1:numel(allResults(i).tests)
        t = allResults(i).tests(j);
        if any(strcmp(t.status, {'FAIL', 'ERROR'}))
            fprintf('  [%s] %s - %s\n', t.status, t.name, t.detail);
        end
    end
end


%% =========================================================================
%  HELPER FUNCTIONS
%% =========================================================================

function tests = add_result(tests, name, status, detail)
%ADD_RESULT Append one test result to the session's list and print it.
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
%PRINT_EXAMPLES Print up to 5 offending file names / pair keys under a result line.
    for i = 1:min(5, numel(names))
        [~, nm, ext] = fileparts(names{i});
        fprintf('      %s: %s%s\n', label, nm, ext);
    end
    if numel(names) > 5
        fprintf('      ... and %d more\n', numel(names) - 5);
    end
end

function dose = load_total_dose(filePath)
%LOAD_TOTAL_DOSE Load a Step 1.5 total dose (sparse or dense storage) as a dense 3D array.
    L = load(filePath);
    if isfield(L, 'total_rs_dose_sparse')
        dose = reshape(full(L.total_rs_dose_sparse), L.total_rs_dose_dims);
    elseif isfield(L, 'ct_total_sparse')
        dose = reshape(full(L.ct_total_sparse), L.ct_total_dims);
    elseif isfield(L, 'total_rs_dose')
        dose = L.total_rs_dose;
    else
        dose = L.ct_total;
    end
end

function err = rel_max_diff(summed, saved)
%REL_MAX_DIFF max|summed - saved| / max|saved|, or Inf when the arrays differ in size.
    if ~isequal(size(summed), size(saved))
        err = Inf;
        return;
    end
    err = max(abs(summed(:) - saved(:))) / max(abs(saved(:)));
end

function c = dose_centroid_mm(dose, origin, spacing)
%DOSE_CENTROID_MM Dose-weighted center of mass in patient coordinates [x y z] (mm).
%   Array layout is (rows = y, cols = x, slices = z).
    total = sum(dose(:));
    xIdx = sum(squeeze(sum(sum(dose, 1), 3)) .* (1:size(dose, 2))) / total;
    yIdx = sum(squeeze(sum(sum(dose, 2), 3)) .* (1:size(dose, 1))') / total;
    zIdx = sum(squeeze(sum(sum(dose, 1), 2)) .* (1:size(dose, 3))') / total;
    c = origin + ([xIdx yIdx zIdx] - 1) .* spacing;
end

function plot_three_views(row, hu, bodyMask, dose, center, spacing, doseMax, CONFIG, rowLabel)
%PLOT_THREE_VIEWS Fill one row of a 2x3 figure: transverse, coronal, sagittal.
%
%   Slices pass through center = [row col slice]. CBCT in grayscale, dose (skipped
%   when empty) as a translucent color wash above the gamma cutoff on a fixed
%   [0 doseMax] scale, body contour in green.
%   Transverse: anterior at top, patient right on the left.
%   Coronal / sagittal: superior at top; sagittal has anterior on the left.
    viewNames = {'Transverse', 'Coronal', 'Sagittal'};
    % Voxel size (mm) along each view's [horizontal, vertical] screen axis
    viewSpacing = [spacing(1) spacing(2); spacing(1) spacing(3); spacing(2) spacing(3)];

    for v = 1:3
        ax = subplot(2, 3, (row - 1) * 3 + v);

        ctSlice = get_slice(hu, v, center);
        ctGray  = min(max((ctSlice - CONFIG.ct_window_hu(1)) / diff(CONFIG.ct_window_hu), 0), 1);
        image(ax, repmat(ctGray, 1, 1, 3));   % truecolor, so the colormap only colors the dose
        hold(ax, 'on');

        if ~isempty(dose)
            doseSlice = get_slice(dose, v, center);
            h = imagesc(ax, doseSlice, [0 doseMax]);
            set(h, 'AlphaData', 0.5 * (doseSlice >= CONFIG.gamma_cutoff_fraction * doseMax));
            colormap(ax, jet);
            if v == 3
                cb = colorbar(ax);
                cb.Label.String = 'Dose (Gy)';
            end
        end

        bodySlice = get_slice(bodyMask, v, center);
        if any(bodySlice(:))
            contour(ax, bodySlice, [0.5 0.5], 'g', 'LineWidth', 1);
        end
        hold(ax, 'off');

        daspect(ax, [viewSpacing(v, 2) viewSpacing(v, 1) 1]);   % true physical proportions
        if v > 1
            set(ax, 'YDir', 'normal');   % higher slice index (superior) at the top
        end
        axis(ax, 'off');
        title(ax, sprintf('%s: %s', rowLabel, viewNames{v}), 'Interpreter', 'none');
    end
end

function s = get_slice(vol, viewNum, center)
%GET_SLICE 2D slice of a (rows = y, cols = x, slices = z) volume through center.
%   1 = transverse (y down, x across), 2 = coronal (z up, x across),
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
