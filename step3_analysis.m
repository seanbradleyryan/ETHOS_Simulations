function results = step3_analysis(patient_id, session, config)
%% STEP3_ANALYSIS - Gamma analysis, SSIM, and visualization
%
%   results = step3_analysis(patient_id, session, config)
%
%   PURPOSE:
%   Quantify agreement between dose distributions using gamma analysis
%   (3%/3mm default) and Structural Similarity Index (SSIM). Performs two
%   comparisons:
%     1. ETHOS truth dose vs Raystation recalculation
%     2. Raystation dose vs Photoacoustic reconstruction
%   Optionally generates visualization figures for quality review.
%
%   INPUTS:
%       patient_id  - String, patient identifier (e.g., '1194203')
%       session     - String, session name (e.g., 'Session_1')
%       config      - Struct with configuration parameters:
%           .working_dir           - Base directory path (REQUIRED)
%           .treatment_site        - Subfolder name (default: 'Pancreas')
%           .gruneisen_method      - Simulation method subfolder (default: 'threshold_2')
%           .gamma_dose_pct        - Dose difference % (default: 3.0)
%           .gamma_dist_mm         - Distance-to-agreement mm (default: 3.0)
%           .gamma_dose_cutoff_pct - Low-dose cutoff, % of max (default: 10.0)
%           .analysis_compute_ssim - Enable SSIM computation (default: true)
%           .analysis_plot_results - Enable figure generation (default: false)
%           .analysis_plot_slices  - Slice indices or 'auto' (default: 'auto')
%
%   OUTPUTS:
%       results - Struct containing:
%           .ethos_vs_rs.gamma     - Gamma analysis struct
%           .ethos_vs_rs.ssim      - SSIM analysis struct (if enabled)
%           .rs_vs_recon.gamma     - Gamma analysis struct
%           .rs_vs_recon.ssim      - SSIM analysis struct (if enabled)
%           .metadata              - Analysis metadata
%
%   FILES CREATED (in AnalysisResults/[PatientID]/[Session]/):
%       - gamma_ethos_vs_rs.mat
%       - gamma_rs_vs_recon.mat
%       - analysis_results.mat (combined results)
%       - figures/ (if analysis_plot_results = true)
%           - gamma_ethos_vs_rs.png
%           - gamma_rs_vs_recon.png
%           - dose_comparison_ethos_rs.png
%           - dose_comparison_rs_recon.png
%
%   PIPELINE INTEGRATION:
%       This function is called from ethos_master_pipeline.m at Step 3.
%       It uses the same CONFIG struct and directory conventions as all
%       other pipeline steps. Helper functions defined here
%       (load_ethos_truth_dose, resample_dose_to_grid,
%       compute_gamma_analysis, save_gamma_results) match the stub
%       signatures in the master pipeline pseudocode.
%
%   DEPENDENCIES:
%       - CalcGamma.m (Mark Geurts gamma calculation)
%       - load_processed_data.m (pipeline Step 1.5 loader)
%       - Image Processing Toolbox (optional, for built-in ssim)
%
%   EXAMPLE:
%       config.working_dir = '/mnt/weka/home/80030361/ETHOS_Simulations';
%       config.treatment_site = 'Pancreas';
%       config.gamma_dose_pct = 3.0;
%       config.gamma_dist_mm = 3.0;
%       config.analysis_plot_results = true;
%       results = step3_analysis('1194203', 'Session_1', config);
%
%   AUTHOR: ETHOS Pipeline Team
%   DATE: February 2026
%   VERSION: 1.0
%
%   See also: CalcGamma, load_processed_data, ethos_master_pipeline

%% ======================== INPUT VALIDATION ========================

if ~ischar(patient_id) && ~isstring(patient_id)
    error('step3_analysis:InvalidInput', ...
        'patient_id must be a string or character array. Received: %s', class(patient_id));
end
patient_id = char(patient_id);

if ~ischar(session) && ~isstring(session)
    error('step3_analysis:InvalidInput', ...
        'session must be a string or character array. Received: %s', class(session));
end
session = char(session);

if ~isstruct(config)
    error('step3_analysis:InvalidInput', ...
        'config must be a struct. Received: %s', class(config));
end

if ~isfield(config, 'working_dir')
    error('step3_analysis:MissingConfig', ...
        'config.working_dir is required but not provided.');
end

%% ======================== SET DEFAULTS ========================

if ~isfield(config, 'treatment_site'),         config.treatment_site = 'Pancreas'; end
if ~isfield(config, 'gruneisen_method'),        config.gruneisen_method = 'threshold_2'; end
if ~isfield(config, 'gamma_dose_pct'),          config.gamma_dose_pct = 3.0; end
if ~isfield(config, 'gamma_dist_mm'),           config.gamma_dist_mm = 3.0; end
if ~isfield(config, 'gamma_dose_cutoff_pct'),   config.gamma_dose_cutoff_pct = 10.0; end
if ~isfield(config, 'analysis_compute_ssim'),    config.analysis_compute_ssim = true; end
if ~isfield(config, 'analysis_plot_results'),    config.analysis_plot_results = false; end
if ~isfield(config, 'analysis_plot_slices'),     config.analysis_plot_slices = 'auto'; end

%% ======================== PRINT HEADER ========================

fprintf('\n=========================================================\n');
fprintf('  [STEP 3] Gamma Analysis, SSIM & Visualization\n');
fprintf('  Patient: %s | Session: %s\n', patient_id, session);
fprintf('  Gamma criteria: %.1f%% / %.1f mm (cutoff: %.1f%% of max)\n', ...
    config.gamma_dose_pct, config.gamma_dist_mm, config.gamma_dose_cutoff_pct);
fprintf('  SSIM: %s | Plotting: %s\n', ...
    bool2str(config.analysis_compute_ssim), bool2str(config.analysis_plot_results));
fprintf('=========================================================\n\n');

analysis_timer = tic;

%% ======================== CREATE OUTPUT DIRECTORY ========================

analysis_dir = get_analysis_directory(patient_id, session, config);
if ~exist(analysis_dir, 'dir')
    mkdir(analysis_dir);
    fprintf('  Created output directory: %s\n', analysis_dir);
end

%% ======================== LOAD DOSE DISTRIBUTIONS ========================

fprintf('  Loading dose distributions...\n');

% --- Load ETHOS truth dose (uses pipeline-matching signature) ---
fprintf('    Loading ETHOS truth dose...\n');
ethos_truth_dose = load_ethos_truth_dose(patient_id, session, config);
fprintf('      Size: [%s], Max: %.4f Gy\n', ...
    num2str(size(ethos_truth_dose), '%d x '), max(ethos_truth_dose(:)));

% --- Load Raystation processed data ---
fprintf('    Loading Raystation total dose and metadata...\n');
[~, ~, total_rs_dose, dose_metadata] = load_processed_data(patient_id, session, config);
fprintf('      Size: [%s], Max: %.4f Gy\n', ...
    num2str(size(total_rs_dose), '%d x '), max(total_rs_dose(:)));
fprintf('      Spacing: [%.2f, %.2f, %.2f] mm\n', ...
    dose_metadata.spacing(1), dose_metadata.spacing(2), dose_metadata.spacing(3));

% --- Load reconstructed dose ---
fprintf('    Loading reconstructed dose...\n');
[total_recon, ~] = load_total_reconstruction(patient_id, session, config);
fprintf('      Size: [%s], Max: %.4f Gy\n', ...
    num2str(size(total_recon), '%d x '), max(total_recon(:)));

%% ======================== DIMENSION MATCHING ========================
% Matches pipeline logic (lines 280-283): resample if sizes differ from RS grid

fprintf('\n  Checking dimension alignment...\n');

rs_size = size(total_rs_dose);

% --- Resample ETHOS if needed (matches pipeline line 280-283) ---
if ~isequal(size(ethos_truth_dose), rs_size)
    fprintf('    [WARNING] ETHOS dose size %s differs from RS dose size %s\n', ...
        mat2str(size(ethos_truth_dose)), mat2str(rs_size));
    fprintf('    Resampling ETHOS dose to match grid...\n');
    ethos_truth_dose = resample_dose_to_grid(ethos_truth_dose, rs_size, dose_metadata);
    fprintf('      Resampled size: %s, Max: %.4f Gy\n', ...
        mat2str(size(ethos_truth_dose)), max(ethos_truth_dose(:)));
end

% --- Resample reconstructed dose if needed ---
if ~isequal(size(total_recon), rs_size)
    fprintf('    [WARNING] Recon dose size %s differs from RS dose size %s\n', ...
        mat2str(size(total_recon)), mat2str(rs_size));
    fprintf('    Resampling recon dose to match grid...\n');
    total_recon = resample_dose_to_grid(total_recon, rs_size, dose_metadata);
    fprintf('      Resampled size: %s, Max: %.4f Gy\n', ...
        mat2str(size(total_recon)), max(total_recon(:)));
end

fprintf('    All distributions aligned: [%d x %d x %d]\n', rs_size(1), rs_size(2), rs_size(3));

%% ======================== GAMMA 1: ETHOS vs RAYSTATION ========================

fprintf('\n  Computing gamma: ETHOS vs Raystation...\n');
gamma_ethos_vs_rs = compute_gamma_analysis(...
    ethos_truth_dose, total_rs_dose, dose_metadata.spacing, config);
fprintf('    Pass rate: %.1f%% (%d/%d voxels)\n', ...
    gamma_ethos_vs_rs.pass_rate, gamma_ethos_vs_rs.num_passed, gamma_ethos_vs_rs.num_evaluated);
fprintf('    Mean gamma: %.3f, Max gamma: %.3f\n', ...
    gamma_ethos_vs_rs.mean_gamma, gamma_ethos_vs_rs.max_gamma);

results.ethos_vs_rs.gamma = gamma_ethos_vs_rs;

%% ======================== GAMMA 2: RS vs RECONSTRUCTED ========================

fprintf('\n  Computing gamma: Raystation vs Reconstructed...\n');
gamma_rs_vs_recon = compute_gamma_analysis(...
    total_rs_dose, total_recon, dose_metadata.spacing, config);
fprintf('    Pass rate: %.1f%% (%d/%d voxels)\n', ...
    gamma_rs_vs_recon.pass_rate, gamma_rs_vs_recon.num_passed, gamma_rs_vs_recon.num_evaluated);
fprintf('    Mean gamma: %.3f, Max gamma: %.3f\n', ...
    gamma_rs_vs_recon.mean_gamma, gamma_rs_vs_recon.max_gamma);

results.rs_vs_recon.gamma = gamma_rs_vs_recon;

%% ======================== SSIM (OPTIONAL) ========================

if config.analysis_compute_ssim
    fprintf('\n  Computing SSIM...\n');
    
    fprintf('    ETHOS vs Raystation...\n');
    results.ethos_vs_rs.ssim = compute_dose_ssim(ethos_truth_dose, total_rs_dose, config);
    fprintf('      3D SSIM: %.4f  |  Mean-slice SSIM: %.4f\n', ...
        results.ethos_vs_rs.ssim.ssim_3d, results.ethos_vs_rs.ssim.ssim_mean_slice);
    
    fprintf('    Raystation vs Reconstructed...\n');
    results.rs_vs_recon.ssim = compute_dose_ssim(total_rs_dose, total_recon, config);
    fprintf('      3D SSIM: %.4f  |  Mean-slice SSIM: %.4f\n', ...
        results.rs_vs_recon.ssim.ssim_3d, results.rs_vs_recon.ssim.ssim_mean_slice);
else
    results.ethos_vs_rs.ssim = [];
    results.rs_vs_recon.ssim = [];
end

%% ======================== PLOTTING (OPTIONAL) ========================

if config.analysis_plot_results
    fprintf('\n  Generating visualization figures...\n');
    
    fig_dir = fullfile(analysis_dir, 'figures');
    if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
    
    plot_analysis_results(ethos_truth_dose, total_rs_dose, ...
        results.ethos_vs_rs.gamma, results.ethos_vs_rs.ssim, ...
        dose_metadata.spacing, 'ETHOS_vs_RS', fig_dir, config);
    
    plot_analysis_results(total_rs_dose, total_recon, ...
        results.rs_vs_recon.gamma, results.rs_vs_recon.ssim, ...
        dose_metadata.spacing, 'RS_vs_Recon', fig_dir, config);
    
    fprintf('    Figures saved to: %s\n', fig_dir);
end

%% ======================== SAVE RESULTS ========================

fprintf('\n  Saving results...\n');

% Populate metadata
results.metadata.patient_id = patient_id;
results.metadata.session = session;
results.metadata.timestamp = datetime('now');
results.metadata.config = config;
results.metadata.spacing_mm = dose_metadata.spacing;
results.metadata.grid_size = rs_size;
results.metadata.ethos_dose_max_Gy = max(ethos_truth_dose(:));
results.metadata.rs_dose_max_Gy = max(total_rs_dose(:));
results.metadata.recon_dose_max_Gy = max(total_recon(:));

% Save individual gamma files (matches save_gamma_results stub signature)
save_gamma_results(gamma_ethos_vs_rs, gamma_rs_vs_recon, patient_id, session, config);

% Save combined results file
analysis_results = results; %#ok<NASGU>
save(fullfile(analysis_dir, 'analysis_results.mat'), 'analysis_results', '-v7.3');
fprintf('    Saved: analysis_results.mat\n');

%% ======================== SUMMARY ========================

elapsed = toc(analysis_timer);
fprintf('\n=========================================================\n');
fprintf('  [STEP 3] Analysis Complete (%.1f sec)\n', elapsed);
fprintf('  ETHOS vs RS:    Gamma pass = %.1f%%', gamma_ethos_vs_rs.pass_rate);
if config.analysis_compute_ssim
    fprintf('  |  SSIM = %.4f', results.ethos_vs_rs.ssim.ssim_3d);
end
fprintf('\n');
fprintf('  RS vs Recon:    Gamma pass = %.1f%%', gamma_rs_vs_recon.pass_rate);
if config.analysis_compute_ssim
    fprintf('  |  SSIM = %.4f', results.rs_vs_recon.ssim.ssim_3d);
end
fprintf('\n');
fprintf('  Output: %s\n', analysis_dir);
fprintf('=========================================================\n\n');

end


%% =========================================================================
%  HELPER FUNCTIONS - Directory / Path Utilities
%  (Match signatures from ethos_master_pipeline helper functions)
%  =========================================================================

function analysis_dir = get_analysis_directory(patient_id, session, config)
    % Return path to analysis output directory
    analysis_dir = fullfile(config.working_dir, 'AnalysisResults', patient_id, session);
    if ~exist(analysis_dir, 'dir')
        mkdir(analysis_dir);
    end
end

function sct_dir = get_sct_directory(patient_id, session, config)
    % Return path to SCT directory
    sct_dir = fullfile(config.working_dir, 'EthosExports', patient_id, ...
        config.treatment_site, session, 'sct');
end

function sim_dir = get_simulation_directory(patient_id, session, config)
    % Return path to simulation output directory
    sim_dir = fullfile(config.working_dir, 'SimulationResults', patient_id, ...
        session, config.gruneisen_method);
end

function s = bool2str(val)
    if val, s = 'enabled'; else, s = 'disabled'; end
end


%% =========================================================================
%  DATA LOADING FUNCTIONS
%  (Signatures match master pipeline stubs exactly)
%  =========================================================================

function ethos_dose = load_ethos_truth_dose(patient_id, session, config)
%% LOAD_ETHOS_TRUTH_DOSE - Load ETHOS RTDOSE (truth) from sct directory
%
%   ethos_dose = load_ethos_truth_dose(patient_id, session, config)
%
%   Matches master pipeline stub signature (line 498): single output.
%   Reads the first RTDOSE*.dcm in the sct directory and applies DoseGridScaling.

    sct_dir = get_sct_directory(patient_id, session, config);
    rd_files = dir(fullfile(sct_dir, 'RTDOSE*.dcm'));
    
    if isempty(rd_files)
        error('step3_analysis:FileNotFound', ...
            'No RTDOSE file (RTDOSE*.dcm) found in: %s', sct_dir);
    end
    
    rd_path = fullfile(sct_dir, rd_files(1).name);
    dose_info = dicominfo(rd_path);
    ethos_dose = double(squeeze(dicomread(rd_path)));
    
    % Apply DoseGridScaling to convert raw values to Gy
    if isfield(dose_info, 'DoseGridScaling')
        ethos_dose = ethos_dose * dose_info.DoseGridScaling;
    end
    
    fprintf('      Loaded: %s\n', rd_files(1).name);
end


function [total_recon, recon_metadata] = load_total_reconstruction(patient_id, session, config)
%% LOAD_TOTAL_RECONSTRUCTION - Load total reconstructed dose from Step 2
%
%   [total_recon, metadata] = load_total_reconstruction(patient_id, session, config)
%
%   Matches master pipeline stub signature (line 490).

    sim_dir = get_simulation_directory(patient_id, session, config);
    recon_file = fullfile(sim_dir, 'total_recon_dose.mat');
    
    if ~isfile(recon_file)
        error('step3_analysis:FileNotFound', ...
            'total_recon_dose.mat not found in: %s\nRun Step 2 (k-Wave simulation) first.', sim_dir);
    end
    
    data = load(recon_file);
    total_recon = data.total_recon;
    
    if isfield(data, 'metadata')
        recon_metadata = data.metadata;
    else
        recon_metadata = struct();
    end
end


%% =========================================================================
%  RESAMPLING
%  (Signature matches master pipeline stub at line 515:
%   resample_dose_to_grid(dose, target_size, metadata)  — 3 args)
%  =========================================================================

function resampled_dose = resample_dose_to_grid(dose, target_size, metadata) %#ok<INUSD>
%% RESAMPLE_DOSE_TO_GRID - Resample a 3D dose array to a target grid size
%
%   resampled_dose = resample_dose_to_grid(dose, target_size, metadata)
%
%   Uses trilinear interpolation. Matches the master pipeline stub signature.
%
%   INPUTS:
%       dose        - 3D source dose array
%       target_size - [nx, ny, nz] target dimensions (e.g. from size(total_rs_dose))
%       metadata    - Struct with .spacing and .dimensions of the target grid
%                     (available for future coordinate-aware resampling)
%
%   OUTPUT:
%       resampled_dose - 3D array of size target_size

    src_size = size(dose);
    
    % If already correct size, no-op
    if isequal(src_size, target_size)
        resampled_dose = dose;
        return;
    end
    
    % Use imresize3 if available (Image Processing Toolbox R2017a+)
    if exist('imresize3', 'file')
        resampled_dose = imresize3(double(dose), target_size, 'linear');
    else
        % Manual trilinear interpolation via interp3
        % interp3 expects (Y, X, Z) ordering for meshgrid convention
        [Xi, Yi, Zi] = ndgrid(...
            linspace(1, src_size(1), target_size(1)), ...
            linspace(1, src_size(2), target_size(2)), ...
            linspace(1, src_size(3), target_size(3)));
        resampled_dose = interp3(double(dose), Yi, Xi, Zi, 'linear', 0);
    end
    
    % Clamp negative values from interpolation
    resampled_dose(resampled_dose < 0) = 0;
end


%% =========================================================================
%  GAMMA ANALYSIS
%  (Signature matches master pipeline stub at line 520:
%   compute_gamma_analysis(reference, evaluated, spacing, config)  — 4 args)
%  =========================================================================

function gamma = compute_gamma_analysis(reference, evaluated, spacing, config)
%% COMPUTE_GAMMA_ANALYSIS - Gamma index between two 3D dose distributions
%
%   gamma = compute_gamma_analysis(reference, evaluated, spacing, config)
%
%   Wraps CalcGamma.m (Mark Geurts) using the required struct format.
%   Applies a low-dose cutoff mask to exclude clinically irrelevant voxels.
%
%   INPUTS:
%       reference - 3D reference dose array (Gy)
%       evaluated - 3D evaluated dose array (Gy), same size as reference
%       spacing   - [dx, dy, dz] voxel spacing in mm
%       config    - Config struct with:
%           .gamma_dose_pct        - Dose difference %
%           .gamma_dist_mm         - DTA in mm
%           .gamma_dose_cutoff_pct - Ignore voxels below this % of max
%
%   OUTPUT:
%       gamma - Struct with fields:
%           .pass_rate        - % of voxels with gamma <= 1
%           .mean_gamma       - Mean gamma in evaluated region
%           .max_gamma        - Maximum gamma value
%           .gamma_map        - 3D gamma map (NaN where masked)
%           .num_evaluated    - Number of evaluated voxels
%           .num_passed       - Number passing (gamma <= 1)
%           .computation_time_sec - CalcGamma wall time
%           .parameters       - Struct echoing input criteria

    % Validate sizes match
    if ~isequal(size(reference), size(evaluated))
        error('compute_gamma_analysis:SizeMismatch', ...
            'Reference size %s does not match evaluated size %s', ...
            mat2str(size(reference)), mat2str(size(evaluated)));
    end
    
    dose_pct   = config.gamma_dose_pct;
    dist_mm    = config.gamma_dist_mm;
    cutoff_pct = config.gamma_dose_cutoff_pct;
    
    % Low-dose cutoff mask: ignore voxels below cutoff_pct of reference max
    max_ref_dose    = max(reference(:));
    dose_threshold  = max_ref_dose * cutoff_pct / 100;
    mask = (reference >= dose_threshold) | (evaluated >= dose_threshold);
    
    fprintf('      Dose threshold: %.4f Gy (%.1f%% of %.4f Gy)\n', ...
        dose_threshold, cutoff_pct, max_ref_dose);
    fprintf('      Voxels above threshold: %d / %d (%.1f%%)\n', ...
        sum(mask(:)), numel(mask), 100 * sum(mask(:)) / numel(mask));
    
    % --- Build CalcGamma input structures ---
    % CalcGamma expects: struct with .start (1xN), .width (1xN), .data (ND)
    ref_struct.start = [0, 0, 0];   % Both on same grid, use relative coords
    ref_struct.width = spacing;      % [dx, dy, dz] in mm
    ref_struct.data  = reference;
    
    tar_struct.start = [0, 0, 0];
    tar_struct.width = spacing;
    tar_struct.data  = evaluated;
    
    % --- Call CalcGamma ---
    % Global gamma, restricted 3D search (axes only), res=20, limit=2
    fprintf('      Running CalcGamma (%.0f%%/%.0fmm, restricted 3D search)...\n', ...
        dose_pct, dist_mm);
    gamma_timer = tic;
    
    gamma_map_raw = CalcGamma(ref_struct, tar_struct, dose_pct, dist_mm, ...
        'local', 0, 'restrict', 1, 'res', 20, 'limit', 2);
    
    gamma_elapsed = toc(gamma_timer);
    fprintf('      CalcGamma completed in %.1f sec\n', gamma_elapsed);
    
    % --- Apply low-dose mask ---
    gamma_map_masked = gamma_map_raw;
    gamma_map_masked(~mask) = NaN;
    
    % --- Compute statistics on masked region ---
    valid_gammas = gamma_map_masked(mask);
    
    gamma.pass_rate     = 100 * sum(valid_gammas <= 1) / numel(valid_gammas);
    gamma.mean_gamma    = mean(valid_gammas);
    gamma.max_gamma     = max(valid_gammas);
    gamma.gamma_map     = gamma_map_masked;
    gamma.num_evaluated = numel(valid_gammas);
    gamma.num_passed    = sum(valid_gammas <= 1);
    gamma.computation_time_sec = gamma_elapsed;
    
    gamma.parameters.dose_pct              = dose_pct;
    gamma.parameters.dist_mm               = dist_mm;
    gamma.parameters.cutoff_pct            = cutoff_pct;
    gamma.parameters.dose_threshold_Gy     = dose_threshold;
    gamma.parameters.max_ref_dose_Gy       = max_ref_dose;
end


%% =========================================================================
%  SSIM COMPUTATION
%  =========================================================================

function ssim_results = compute_dose_ssim(reference, evaluated, config) %#ok<INUSD>
%% COMPUTE_DOSE_SSIM - Structural Similarity Index for 3D dose distributions
%
%   ssim_results = compute_dose_ssim(reference, evaluated, config)
%
%   Computes per-slice (axial) SSIM and a dose-weighted 3D aggregate.
%   Both arrays are normalized to a shared dynamic range before comparison.
%
%   INPUTS:
%       reference - 3D dose array (Gy)
%       evaluated - 3D dose array (Gy), same size as reference
%       config    - Config struct (reserved for future options)
%
%   OUTPUT:
%       ssim_results - Struct with:
%           .ssim_3d           - Dose-weighted aggregate SSIM
%           .ssim_per_slice    - Vector of per-axial-slice SSIM
%           .ssim_mean_slice   - Simple mean of per-slice SSIM
%           .dynamic_range     - [min, max] used for normalization
%           .num_slices        - Number of axial slices

    if ~isequal(size(reference), size(evaluated))
        error('compute_dose_ssim:SizeMismatch', ...
            'Reference size %s does not match evaluated size %s', ...
            mat2str(size(reference)), mat2str(size(evaluated)));
    end
    
    % Shared dynamic range for normalization
    combined_max = max(max(reference(:)), max(evaluated(:)));
    combined_min = min(min(reference(:)), min(evaluated(:)));
    dynamic_range = [combined_min, combined_max];
    
    % Normalize to [0, 1]
    if combined_max > combined_min
        ref_norm  = (reference - combined_min) / (combined_max - combined_min);
        eval_norm = (evaluated - combined_min) / (combined_max - combined_min);
    else
        ref_norm  = zeros(size(reference));
        eval_norm = zeros(size(evaluated));
    end
    
    % Per-slice SSIM (3rd dimension = axial slices)
    num_slices = size(reference, 3);
    ssim_per_slice = zeros(num_slices, 1);
    has_builtin_ssim = (exist('ssim', 'file') == 2);
    
    for k = 1:num_slices
        ref_slice  = ref_norm(:, :, k);
        eval_slice = eval_norm(:, :, k);
        
        % Both empty -> perfect match
        if max(ref_slice(:)) < 1e-10 && max(eval_slice(:)) < 1e-10
            ssim_per_slice(k) = 1.0;
            continue;
        end
        
        if has_builtin_ssim
            ssim_per_slice(k) = ssim(eval_slice, ref_slice, 'DynamicRange', 1.0);
        else
            ssim_per_slice(k) = compute_ssim_manual(ref_slice, eval_slice);
        end
    end
    
    % Dose-weighted 3D SSIM (weight each slice by its max dose)
    slice_weights = zeros(num_slices, 1);
    for k = 1:num_slices
        slice_weights(k) = max(max(reference(:, :, k)));
    end
    if sum(slice_weights) > 0
        slice_weights_norm = slice_weights / sum(slice_weights);
        ssim_3d = sum(ssim_per_slice .* slice_weights_norm);
    else
        ssim_3d = mean(ssim_per_slice);
    end
    
    % Pack output
    ssim_results.ssim_3d         = ssim_3d;
    ssim_results.ssim_per_slice  = ssim_per_slice;
    ssim_results.ssim_mean_slice = mean(ssim_per_slice);
    ssim_results.dynamic_range   = dynamic_range;
    ssim_results.num_slices      = num_slices;
end


function val = compute_ssim_manual(ref, eval_img)
%% COMPUTE_SSIM_MANUAL - Fallback SSIM for a 2D image pair
%   Standard SSIM with 11x11 Gaussian window, K1=0.01, K2=0.03, L=1.0

    K1 = 0.01;  K2 = 0.03;  L = 1.0;
    C1 = (K1 * L)^2;
    C2 = (K2 * L)^2;
    
    % Gaussian window
    win_size = 11; sigma = 1.5;
    [x, y] = meshgrid(-(win_size-1)/2:(win_size-1)/2);
    g = exp(-(x.^2 + y.^2) / (2 * sigma^2));
    g = g / sum(g(:));
    
    mu_ref  = conv2(ref, g, 'valid');
    mu_eval = conv2(eval_img, g, 'valid');
    
    sigma_ref_sq   = conv2(ref.^2, g, 'valid')         - mu_ref.^2;
    sigma_eval_sq  = conv2(eval_img.^2, g, 'valid')     - mu_eval.^2;
    sigma_ref_eval = conv2(ref .* eval_img, g, 'valid')  - mu_ref .* mu_eval;
    
    numerator   = (2 * mu_ref .* mu_eval + C1) .* (2 * sigma_ref_eval + C2);
    denominator = (mu_ref.^2 + mu_eval.^2 + C1) .* (sigma_ref_sq + sigma_eval_sq + C2);
    
    ssim_map = numerator ./ denominator;
    val = mean(ssim_map(:));
end


%% =========================================================================
%  SAVE RESULTS
%  (Signature matches master pipeline stub at line 526:
%   save_gamma_results(gamma1, gamma2, patient_id, session, config)  — 5 args)
%  =========================================================================

function save_gamma_results(gamma_ethos_vs_rs, gamma_rs_vs_recon, patient_id, session, config)
%% SAVE_GAMMA_RESULTS - Save gamma analysis results to disk
%
%   save_gamma_results(gamma_ethos_vs_rs, gamma_rs_vs_recon, patient_id, session, config)
%
%   Matches the master pipeline stub signature.

    analysis_dir = get_analysis_directory(patient_id, session, config);
    
    save(fullfile(analysis_dir, 'gamma_ethos_vs_rs.mat'), 'gamma_ethos_vs_rs', '-v7.3');
    fprintf('    Saved: gamma_ethos_vs_rs.mat\n');
    
    save(fullfile(analysis_dir, 'gamma_rs_vs_recon.mat'), 'gamma_rs_vs_recon', '-v7.3');
    fprintf('    Saved: gamma_rs_vs_recon.mat\n');
end


%% =========================================================================
%  VISUALIZATION
%  =========================================================================

function plot_analysis_results(reference, evaluated, gamma_result, ssim_result, spacing, comparison_name, output_dir, config) %#ok<INUSD>
%% PLOT_ANALYSIS_RESULTS - Generate dose comparison and gamma figures
%
%   Creates two figures per comparison:
%     1. Dose comparison: reference, evaluated, and difference side-by-side
%     2. Gamma map: gamma overlay, pass/fail map, and histogram
%
%   INPUTS:
%       reference       - 3D reference dose (Gy)
%       evaluated       - 3D evaluated dose (Gy)
%       gamma_result    - Struct from compute_gamma_analysis
%       ssim_result     - Struct from compute_dose_ssim (may be [])
%       spacing         - [dx, dy, dz] in mm (unused, reserved)
%       comparison_name - String for titles/filenames (e.g. 'ETHOS_vs_RS')
%       output_dir      - Directory for saving figures
%       config          - Config struct

    plot_slices = select_plot_slices(reference, config);
    num_plot = length(plot_slices);
    title_name = strrep(comparison_name, '_', ' ');
    
    dose_max = max(max(reference(:)), max(evaluated(:)));
    dose_clim = [0, dose_max];
    
    % ===== FIGURE 1: Dose Comparison =====
    fig1 = figure('Position', [100 100 400*num_plot 1000], 'Visible', 'off');
    
    for s = 1:num_plot
        k = plot_slices(s);
        ref_slice  = reference(:, :, k);
        eval_slice = evaluated(:, :, k);
        diff_slice = eval_slice - ref_slice;
        
        % Row 1: Reference
        subplot(3, num_plot, s);
        imagesc(ref_slice, dose_clim);
        colormap(gca, 'jet'); colorbar;
        title(sprintf('Reference (z=%d)', k));
        axis image; axis off;
        
        % Row 2: Evaluated
        subplot(3, num_plot, num_plot + s);
        imagesc(eval_slice, dose_clim);
        colormap(gca, 'jet'); colorbar;
        title(sprintf('Evaluated (z=%d)', k));
        axis image; axis off;
        
        % Row 3: Difference
        subplot(3, num_plot, 2*num_plot + s);
        diff_lim = max(abs(diff_slice(:)));
        if diff_lim < 1e-10, diff_lim = 1; end
        imagesc(diff_slice, [-diff_lim, diff_lim]);
        colormap(gca, rdbu_colormap()); colorbar;
        title(sprintf('Difference (z=%d)', k));
        axis image; axis off;
    end
    
    sgtitle(sprintf('Dose Comparison: %s', title_name), 'FontSize', 14, 'FontWeight', 'bold');
    
    fig_path = fullfile(output_dir, sprintf('dose_comparison_%s.png', lower(comparison_name)));
    exportgraphics(fig1, fig_path, 'Resolution', 150);
    close(fig1);
    fprintf('      Saved: %s\n', fig_path);
    
    % ===== FIGURE 2: Gamma Analysis =====
    fig2 = figure('Position', [100 100 400*num_plot+400 800], 'Visible', 'off');
    
    gamma_map = gamma_result.gamma_map;
    gamma_clim = [0, 2];
    
    for s = 1:num_plot
        k = plot_slices(s);
        gamma_slice = gamma_map(:, :, k);
        
        % Row 1: Gamma map
        subplot(2, num_plot + 1, s);
        imagesc(gamma_slice, gamma_clim);
        colormap(gca, gamma_colormap()); colorbar;
        title(sprintf('Gamma (z=%d)', k));
        axis image; axis off;
        
        % Row 2: Pass/fail overlay
        subplot(2, num_plot + 1, num_plot + 1 + s);
        pass_fail = nan(size(gamma_slice));
        valid = ~isnan(gamma_slice);
        pass_fail(valid & gamma_slice <= 1) = 1;  % Pass = green
        pass_fail(valid & gamma_slice > 1)  = 0;  % Fail = red
        imagesc(pass_fail, [0, 1]);
        colormap(gca, [1 0 0; 0 0.7 0]);
        colorbar('Ticks', [0.25, 0.75], 'TickLabels', {'Fail', 'Pass'});
        title(sprintf('Pass/Fail (z=%d)', k));
        axis image; axis off;
    end
    
    % Histogram in last column
    subplot(2, num_plot + 1, [num_plot + 1, 2*(num_plot + 1)]);
    valid_gammas = gamma_map(~isnan(gamma_map));
    histogram(valid_gammas, 0:0.05:3, 'FaceColor', [0.3 0.5 0.8], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.8);
    hold on; xline(1, 'r--', 'LineWidth', 2); hold off;
    xlabel('Gamma Index'); ylabel('Voxel Count');
    title(sprintf('Pass: %.1f%% | Mean: %.2f', ...
        gamma_result.pass_rate, gamma_result.mean_gamma));
    xlim([0, 3]); grid on;
    
    % Supertitle with SSIM if available
    if ~isempty(ssim_result) && isstruct(ssim_result)
        sgtitle(sprintf('Gamma: %s (%.0f%%/%.0fmm) | SSIM=%.4f', ...
            title_name, gamma_result.parameters.dose_pct, ...
            gamma_result.parameters.dist_mm, ssim_result.ssim_3d), ...
            'FontSize', 14, 'FontWeight', 'bold');
    else
        sgtitle(sprintf('Gamma: %s (%.0f%%/%.0fmm)', ...
            title_name, gamma_result.parameters.dose_pct, ...
            gamma_result.parameters.dist_mm), ...
            'FontSize', 14, 'FontWeight', 'bold');
    end
    
    fig_path = fullfile(output_dir, sprintf('gamma_%s.png', lower(comparison_name)));
    exportgraphics(fig2, fig_path, 'Resolution', 150);
    close(fig2);
    fprintf('      Saved: %s\n', fig_path);
end


function slices = select_plot_slices(dose, config)
    % Select representative axial slices for plotting.
    % If config.analysis_plot_slices is numeric, use directly.
    % If 'auto', pick 3 slices at 25th/50th/75th percentile of dose-bearing range.
    
    if isfield(config, 'analysis_plot_slices') && isnumeric(config.analysis_plot_slices)
        slices = config.analysis_plot_slices;
        slices = slices(slices >= 1 & slices <= size(dose, 3));
        if isempty(slices)
            slices = round(size(dose, 3) / 2);
        end
        return;
    end
    
    % Auto mode
    num_slices = size(dose, 3);
    slice_max = zeros(num_slices, 1);
    for k = 1:num_slices
        slice_max(k) = max(max(dose(:, :, k)));
    end
    
    active_slices = find(slice_max >= 0.1 * max(slice_max));
    
    if isempty(active_slices)
        slices = round(num_slices / 2);
    elseif length(active_slices) <= 3
        slices = active_slices(:)';
    else
        n = length(active_slices);
        slices = unique([...
            active_slices(max(1, round(0.25 * n))), ...
            active_slices(max(1, round(0.50 * n))), ...
            active_slices(max(1, round(0.75 * n)))]);
    end
end


function cmap = rdbu_colormap()
    % Red-white-blue diverging colormap for dose difference
    n = 256; half = n / 2;
    r = [linspace(0.1, 1, half)'; linspace(1, 0.8, half)'];
    g = [linspace(0.2, 1, half)'; linspace(1, 0.1, half)'];
    b = [linspace(0.7, 1, half)'; linspace(1, 0.1, half)'];
    cmap = [r, g, b];
end


function cmap = gamma_colormap()
    % Green -> yellow -> red colormap for gamma display
    n = 256; half = round(n / 2); remainder = n - half;
    r = [linspace(0, 1, half)';    ones(remainder, 1)];
    g = [ones(half, 1) * 0.8;     linspace(0.8, 0, remainder)'];
    b = [linspace(0.2, 0, half)';  zeros(remainder, 1)];
    cmap = [r, g, b];
end
