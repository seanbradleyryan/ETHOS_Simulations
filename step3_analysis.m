function results = step3_analysis(patient_id, session, config)
%% STEP3_ANALYSIS - Gamma analysis, SSIM, and visualization
%
%   results = step3_analysis(patient_id, session, config)
%
%   PURPOSE:
%   Quantify agreement between dose distributions using gamma analysis
%   (3%/3mm criteria) and Structural Similarity Index (SSIM). Performs two
%   comparisons:
%     1. ETHOS truth dose vs Raystation recalculation
%     2. Raystation dose vs Photoacoustic reconstruction
%
%   INPUTS:
%       patient_id  - String, patient identifier (e.g., '1194203')
%       session     - String, session name (e.g., 'Session_1')
%       config      - Struct with configuration parameters:
%           .working_dir           - Base directory path
%           .treatment_site        - Subfolder name (default: 'Pancreas')
%           .gamma_dose_pct        - Dose difference % (default: 3.0)
%           .gamma_dist_mm         - DTA mm (default: 3.0)
%           .gamma_dose_cutoff_pct - Low-dose cutoff % of max (default: 10.0)
%           .analysis_plot_results - Enable plotting (default: false)
%           .analysis_plot_slices  - Slice indices or 'auto' (default: 'auto')
%           .gruneisen_method      - For simulation directory path
%
%   OUTPUTS:
%       results - Struct containing:
%           .ethos_vs_rs.gamma     - Gamma analysis struct
%           .ethos_vs_rs.ssim      - SSIM analysis struct
%           .rs_vs_recon.gamma     - Gamma analysis struct
%           .rs_vs_recon.ssim      - SSIM analysis struct
%           .metadata              - Analysis metadata
%
%   FILES CREATED (in AnalysisResults/[PatientID]/[Session]/):
%       - gamma_ethos_vs_rs.mat
%       - gamma_rs_vs_recon.mat
%       - analysis_results.mat (combined)
%       - figures/ (if analysis_plot_results = true)
%           - gamma_ethos_vs_rs.png
%           - gamma_rs_vs_recon.png
%           - dose_comparison_ethos_rs.png
%           - dose_comparison_rs_recon.png
%
%   DEPENDENCIES:
%       - CalcGamma.m (Mark Geurts gamma calculation)
%       - load_processed_data.m
%       - Image Processing Toolbox (for ssim)
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
%   See also: CalcGamma, load_processed_data, compute_gamma_analysis, compute_dose_ssim

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

if ~isfield(config, 'treatment_site'),        config.treatment_site = 'Pancreas'; end
if ~isfield(config, 'gamma_dose_pct'),         config.gamma_dose_pct = 3.0; end
if ~isfield(config, 'gamma_dist_mm'),          config.gamma_dist_mm = 3.0; end
if ~isfield(config, 'gamma_dose_cutoff_pct'),  config.gamma_dose_cutoff_pct = 10.0; end
if ~isfield(config, 'analysis_plot_results'),   config.analysis_plot_results = false; end
if ~isfield(config, 'analysis_plot_slices'),    config.analysis_plot_slices = 'auto'; end
if ~isfield(config, 'gruneisen_method'),        config.gruneisen_method = 'threshold_2'; end

%% ======================== PRINT HEADER ========================

fprintf('\n=========================================================\n');
fprintf('  [STEP 3] Gamma Analysis and SSIM\n');
fprintf('  Patient: %s | Session: %s\n', patient_id, session);
fprintf('  Gamma criteria: %.1f%% / %.1f mm (cutoff: %.1f%% of max)\n', ...
    config.gamma_dose_pct, config.gamma_dist_mm, config.gamma_dose_cutoff_pct);
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

% --- Load ETHOS truth dose ---
fprintf('    Loading ETHOS truth dose...\n');
[ethos_dose, ethos_info] = load_ethos_truth_dose(patient_id, session, config);
fprintf('      Size: [%d x %d x %d], Max: %.4f Gy\n', ...
    size(ethos_dose, 1), size(ethos_dose, 2), size(ethos_dose, 3), max(ethos_dose(:)));

% --- Load Raystation processed data ---
fprintf('    Loading Raystation total dose...\n');
[~, ~, total_rs_dose, metadata] = load_processed_data(patient_id, session, config);
spacing = metadata.spacing;  % [dx, dy, dz] in mm
fprintf('      Size: [%d x %d x %d], Max: %.4f Gy\n', ...
    size(total_rs_dose, 1), size(total_rs_dose, 2), size(total_rs_dose, 3), max(total_rs_dose(:)));
fprintf('      Spacing: [%.2f, %.2f, %.2f] mm\n', spacing(1), spacing(2), spacing(3));

% --- Load reconstructed dose ---
fprintf('    Loading reconstructed dose...\n');
[total_recon, ~] = load_total_reconstruction(patient_id, session, config);
fprintf('      Size: [%d x %d x %d], Max: %.4f Gy\n', ...
    size(total_recon, 1), size(total_recon, 2), size(total_recon, 3), max(total_recon(:)));

%% ======================== DIMENSION MATCHING ========================

fprintf('\n  Checking dimension alignment...\n');

rs_size = size(total_rs_dose);
ethos_size = size(ethos_dose);
recon_size = size(total_recon);

% Check ETHOS vs RS dimensions
if ~isequal(ethos_size, rs_size)
    fprintf('    [WARNING] ETHOS dose size %s differs from RS dose size %s\n', ...
        mat2str(ethos_size), mat2str(rs_size));
    fprintf('    Resampling ETHOS dose to RS grid...\n');
    ethos_dose = resample_dose_to_grid(ethos_dose, ethos_info, rs_size, metadata);
    fprintf('      Resampled ETHOS dose size: %s, Max: %.4f Gy\n', ...
        mat2str(size(ethos_dose)), max(ethos_dose(:)));
end

% Check Recon vs RS dimensions
if ~isequal(recon_size, rs_size)
    fprintf('    [WARNING] Recon dose size %s differs from RS dose size %s\n', ...
        mat2str(recon_size), mat2str(rs_size));
    fprintf('    Resampling recon dose to RS grid...\n');
    % Recon should already be on the same grid from step 2, but handle if not
    total_recon = resample_3d_array(total_recon, rs_size);
    fprintf('      Resampled recon dose size: %s, Max: %.4f Gy\n', ...
        mat2str(size(total_recon)), max(total_recon(:)));
end

fprintf('    All distributions aligned to grid: [%d x %d x %d]\n', rs_size(1), rs_size(2), rs_size(3));

%% ======================== ANALYSIS 1: ETHOS vs RAYSTATION ========================

fprintf('\n  Computing ETHOS vs Raystation analysis...\n');

fprintf('    Running gamma analysis...\n');
results.ethos_vs_rs.gamma = compute_gamma_analysis(ethos_dose, total_rs_dose, spacing, config);
fprintf('      Pass rate: %.1f%% (%d/%d voxels)\n', ...
    results.ethos_vs_rs.gamma.pass_rate, ...
    results.ethos_vs_rs.gamma.num_passed, ...
    results.ethos_vs_rs.gamma.num_evaluated);
fprintf('      Mean gamma: %.3f, Max gamma: %.3f\n', ...
    results.ethos_vs_rs.gamma.mean_gamma, ...
    results.ethos_vs_rs.gamma.max_gamma);

fprintf('    Running SSIM analysis...\n');
results.ethos_vs_rs.ssim = compute_dose_ssim(ethos_dose, total_rs_dose, config);
fprintf('      3D SSIM: %.4f\n', results.ethos_vs_rs.ssim.ssim_3d);
fprintf('      Mean slice SSIM: %.4f\n', results.ethos_vs_rs.ssim.ssim_mean_slice);

%% ======================== ANALYSIS 2: RS vs RECONSTRUCTED ========================

fprintf('\n  Computing Raystation vs Reconstructed analysis...\n');

fprintf('    Running gamma analysis...\n');
results.rs_vs_recon.gamma = compute_gamma_analysis(total_rs_dose, total_recon, spacing, config);
fprintf('      Pass rate: %.1f%% (%d/%d voxels)\n', ...
    results.rs_vs_recon.gamma.pass_rate, ...
    results.rs_vs_recon.gamma.num_passed, ...
    results.rs_vs_recon.gamma.num_evaluated);
fprintf('      Mean gamma: %.3f, Max gamma: %.3f\n', ...
    results.rs_vs_recon.gamma.mean_gamma, ...
    results.rs_vs_recon.gamma.max_gamma);

fprintf('    Running SSIM analysis...\n');
results.rs_vs_recon.ssim = compute_dose_ssim(total_rs_dose, total_recon, config);
fprintf('      3D SSIM: %.4f\n', results.rs_vs_recon.ssim.ssim_3d);
fprintf('      Mean slice SSIM: %.4f\n', results.rs_vs_recon.ssim.ssim_mean_slice);

%% ======================== OPTIONAL PLOTTING ========================

if config.analysis_plot_results
    fprintf('\n  Generating visualization figures...\n');
    
    fig_dir = fullfile(analysis_dir, 'figures');
    if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
    
    plot_analysis_results(ethos_dose, total_rs_dose, ...
        results.ethos_vs_rs.gamma, results.ethos_vs_rs.ssim, ...
        spacing, 'ETHOS_vs_RS', fig_dir, config);
    
    plot_analysis_results(total_rs_dose, total_recon, ...
        results.rs_vs_recon.gamma, results.rs_vs_recon.ssim, ...
        spacing, 'RS_vs_Recon', fig_dir, config);
    
    fprintf('    Figures saved to: %s\n', fig_dir);
end

%% ======================== STORE METADATA ========================

results.metadata.patient_id = patient_id;
results.metadata.session = session;
results.metadata.timestamp = datetime('now');
results.metadata.config = config;
results.metadata.spacing_mm = spacing;
results.metadata.grid_size = rs_size;
results.metadata.ethos_dose_max_Gy = max(ethos_dose(:));
results.metadata.rs_dose_max_Gy = max(total_rs_dose(:));
results.metadata.recon_dose_max_Gy = max(total_recon(:));

%% ======================== SAVE RESULTS ========================

fprintf('\n  Saving analysis results...\n');

% Save individual gamma results (for backward compatibility with pipeline spec)
gamma_ethos_vs_rs = results.ethos_vs_rs.gamma;
save(fullfile(analysis_dir, 'gamma_ethos_vs_rs.mat'), 'gamma_ethos_vs_rs', '-v7.3');
fprintf('    Saved: gamma_ethos_vs_rs.mat\n');

gamma_rs_vs_recon = results.rs_vs_recon.gamma;
save(fullfile(analysis_dir, 'gamma_rs_vs_recon.mat'), 'gamma_rs_vs_recon', '-v7.3');
fprintf('    Saved: gamma_rs_vs_recon.mat\n');

% Save combined results
save(fullfile(analysis_dir, 'analysis_results.mat'), 'results', '-v7.3');
fprintf('    Saved: analysis_results.mat\n');

%% ======================== SUMMARY ========================

elapsed = toc(analysis_timer);
fprintf('\n=========================================================\n');
fprintf('  [STEP 3] Analysis Complete (%.1f sec)\n', elapsed);
fprintf('  ETHOS vs RS:    Gamma pass = %.1f%%  |  SSIM = %.4f\n', ...
    results.ethos_vs_rs.gamma.pass_rate, results.ethos_vs_rs.ssim.ssim_3d);
fprintf('  RS vs Recon:    Gamma pass = %.1f%%  |  SSIM = %.4f\n', ...
    results.rs_vs_recon.gamma.pass_rate, results.rs_vs_recon.ssim.ssim_3d);
fprintf('  Output directory: %s\n', analysis_dir);
fprintf('=========================================================\n\n');

end


%% ======================= HELPER FUNCTIONS ================================

function analysis_dir = get_analysis_directory(patient_id, session, config)
    % Return path to analysis output directory
    analysis_dir = fullfile(config.working_dir, 'AnalysisResults', patient_id, session);
end


function [ethos_dose, info] = load_ethos_truth_dose(patient_id, session, config)
    % Load ETHOS RTDOSE (truth) from sct directory and return dose + DICOM info
    
    sct_dir = fullfile(config.working_dir, 'EthosExports', patient_id, ...
        config.treatment_site, session, 'sct');
    
    rd_files = dir(fullfile(sct_dir, 'RD*.dcm'));
    if isempty(rd_files)
        error('step3_analysis:FileNotFound', ...
            'No RTDOSE file (RD*.dcm) found in: %s', sct_dir);
    end
    
    % Use first RTDOSE file
    rd_path = fullfile(sct_dir, rd_files(1).name);
    info = dicominfo(rd_path);
    ethos_dose = double(squeeze(dicomread(rd_path)));
    
    % Apply DoseGridScaling to convert to Gy
    if isfield(info, 'DoseGridScaling')
        ethos_dose = ethos_dose * info.DoseGridScaling;
    end
    
    % Extract geometry info for potential resampling
    if isfield(info, 'ImagePositionPatient')
        info.origin = info.ImagePositionPatient(:)';
    end
    if isfield(info, 'PixelSpacing')
        dx = info.PixelSpacing(1);
        dy = info.PixelSpacing(2);
        % Z spacing from GridFrameOffsetVector
        if isfield(info, 'GridFrameOffsetVector') && length(info.GridFrameOffsetVector) > 1
            dz = abs(info.GridFrameOffsetVector(2) - info.GridFrameOffsetVector(1));
        else
            dz = dx;  % fallback
        end
        info.spacing = [dx, dy, dz];
    end
end


function [total_recon, recon_metadata] = load_total_reconstruction(patient_id, session, config)
    % Load total reconstructed dose from simulation results
    
    sim_dir = fullfile(config.working_dir, 'SimulationResults', patient_id, ...
        session, config.gruneisen_method);
    
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


function resampled = resample_dose_to_grid(dose, dose_info, target_size, target_metadata)
    % Resample a dose distribution to match a target grid using trilinear interpolation
    %
    % Uses DICOM geometry info from dose_info and target_metadata to build
    % coordinate grids and interpolate.
    
    src_size = size(dose);
    
    % Build source coordinate vectors (mm)
    if isfield(dose_info, 'origin') && isfield(dose_info, 'spacing')
        src_origin = dose_info.origin;
        src_spacing = dose_info.spacing;
    else
        % Fallback: assume same origin, scale spacing
        src_origin = target_metadata.origin;
        src_spacing = target_metadata.spacing .* (target_size ./ src_size);
    end
    
    src_x = src_origin(1) + (0:src_size(1)-1) * src_spacing(1);
    src_y = src_origin(2) + (0:src_size(2)-1) * src_spacing(2);
    src_z = src_origin(3) + (0:src_size(3)-1) * src_spacing(3);
    
    % Build target coordinate vectors (mm)
    tgt_origin = target_metadata.origin;
    tgt_spacing = target_metadata.spacing;
    
    tgt_x = tgt_origin(1) + (0:target_size(1)-1) * tgt_spacing(1);
    tgt_y = tgt_origin(2) + (0:target_size(2)-1) * tgt_spacing(2);
    tgt_z = tgt_origin(3) + (0:target_size(3)-1) * tgt_spacing(3);
    
    % Create meshgrids
    [SrcX, SrcY, SrcZ] = ndgrid(src_x, src_y, src_z);
    [TgtX, TgtY, TgtZ] = ndgrid(tgt_x, tgt_y, tgt_z);
    
    % Interpolate
    F = griddedInterpolant(SrcX, SrcY, SrcZ, double(dose), 'linear', 'none');
    resampled = F(TgtX, TgtY, TgtZ);
    
    % Replace NaN (outside source domain) with 0
    resampled(isnan(resampled)) = 0;
end


function resampled = resample_3d_array(vol, target_size)
    % Simple resampling of 3D array to target size using imresize3
    % Used when coordinate metadata is unavailable
    
    if exist('imresize3', 'file')
        resampled = imresize3(double(vol), target_size, 'linear');
    else
        % Manual trilinear interpolation fallback
        src_size = size(vol);
        [Xi, Yi, Zi] = ndgrid(...
            linspace(1, src_size(1), target_size(1)), ...
            linspace(1, src_size(2), target_size(2)), ...
            linspace(1, src_size(3), target_size(3)));
        resampled = interp3(double(vol), Yi, Xi, Zi, 'linear', 0);
    end
end


%% ======================= CORE ANALYSIS FUNCTIONS =========================

function gamma = compute_gamma_analysis(reference, evaluated, spacing, config)
%% COMPUTE_GAMMA_ANALYSIS - Gamma analysis between two 3D dose distributions
%
%   Uses CalcGamma.m (Mark Geurts) with the reference/target struct format.
%   Applies a low-dose cutoff mask to exclude clinically irrelevant voxels.
%
%   INPUTS:
%       reference - 3D dose array (Gy), the ground truth
%       evaluated - 3D dose array (Gy), same size as reference
%       spacing   - [dx, dy, dz] voxel spacing in mm
%       config    - Config struct with gamma parameters
%
%   OUTPUT:
%       gamma - Struct with pass_rate, mean_gamma, max_gamma, gamma_map, etc.

    % Validate inputs
    if ~isequal(size(reference), size(evaluated))
        error('compute_gamma_analysis:SizeMismatch', ...
            'Reference size %s does not match evaluated size %s', ...
            mat2str(size(reference)), mat2str(size(evaluated)));
    end
    
    dose_pct = config.gamma_dose_pct;
    dist_mm = config.gamma_dist_mm;
    cutoff_pct = config.gamma_dose_cutoff_pct;
    
    % Build low-dose cutoff mask: ignore voxels below cutoff_pct of reference max
    max_ref_dose = max(reference(:));
    dose_threshold = max_ref_dose * cutoff_pct / 100;
    mask = (reference >= dose_threshold) | (evaluated >= dose_threshold);
    
    fprintf('      Dose threshold: %.4f Gy (%.1f%% of %.4f Gy)\n', ...
        dose_threshold, cutoff_pct, max_ref_dose);
    fprintf('      Voxels above threshold: %d / %d (%.1f%%)\n', ...
        sum(mask(:)), numel(mask), 100 * sum(mask(:)) / numel(mask));
    
    % Build CalcGamma input structures
    % CalcGamma expects: struct with .start, .width, .data
    grid_size = size(reference);
    
    ref_struct.start = [0, 0, 0];  % Relative coordinates are fine
    ref_struct.width = spacing;    % [dx, dy, dz] in mm
    ref_struct.data  = reference;
    
    tar_struct.start = [0, 0, 0];
    tar_struct.width = spacing;
    tar_struct.data  = evaluated;
    
    % Run CalcGamma (global gamma, restricted search for speed)
    fprintf('      Running CalcGamma (%.0f%%/%.0fmm, restricted search)...\n', dose_pct, dist_mm);
    gamma_timer = tic;
    
    gamma_map_raw = CalcGamma(ref_struct, tar_struct, dose_pct, dist_mm, ...
        'local', 0, 'restrict', 1, 'res', 20, 'limit', 2);
    
    gamma_elapsed = toc(gamma_timer);
    fprintf('      CalcGamma completed in %.1f sec\n', gamma_elapsed);
    
    % Apply mask
    gamma_map_masked = gamma_map_raw;
    gamma_map_masked(~mask) = NaN;
    
    % Compute statistics on masked region only
    valid_gammas = gamma_map_masked(mask);
    
    gamma.gamma_map     = gamma_map_masked;
    gamma.pass_rate     = 100 * sum(valid_gammas <= 1) / length(valid_gammas);
    gamma.mean_gamma    = mean(valid_gammas);
    gamma.max_gamma     = max(valid_gammas);
    gamma.num_evaluated = length(valid_gammas);
    gamma.num_passed    = sum(valid_gammas <= 1);
    gamma.computation_time_sec = gamma_elapsed;
    
    gamma.parameters.dose_pct   = dose_pct;
    gamma.parameters.dist_mm    = dist_mm;
    gamma.parameters.cutoff_pct = cutoff_pct;
    gamma.parameters.dose_threshold_Gy = dose_threshold;
    gamma.parameters.max_ref_dose_Gy = max_ref_dose;
end


function ssim_results = compute_dose_ssim(reference, evaluated, config)
%% COMPUTE_DOSE_SSIM - SSIM between two 3D dose distributions
%
%   Computes both an overall 3D SSIM and per-slice (axial) SSIM values.
%   Normalizes both distributions to a common dynamic range before comparison.
%
%   INPUTS:
%       reference - 3D dose array (Gy)
%       evaluated - 3D dose array (Gy), same size as reference
%       config    - Config struct (currently unused, reserved for future options)
%
%   OUTPUT:
%       ssim_results - Struct with ssim_3d, ssim_per_slice, ssim_mean_slice, etc.

    if ~isequal(size(reference), size(evaluated))
        error('compute_dose_ssim:SizeMismatch', ...
            'Reference size %s does not match evaluated size %s', ...
            mat2str(size(reference)), mat2str(size(evaluated)));
    end
    
    % Determine shared dynamic range for normalization
    combined_max = max(max(reference(:)), max(evaluated(:)));
    combined_min = min(min(reference(:)), min(evaluated(:)));
    dynamic_range = [combined_min, combined_max];
    
    % Normalize to [0, 1] for SSIM computation
    if combined_max > combined_min
        ref_norm = (reference - combined_min) / (combined_max - combined_min);
        eval_norm = (evaluated - combined_min) / (combined_max - combined_min);
    else
        ref_norm = zeros(size(reference));
        eval_norm = zeros(size(evaluated));
    end
    
    % --- Per-slice SSIM (axial slices = 3rd dimension) ---
    num_slices = size(reference, 3);
    ssim_per_slice = zeros(num_slices, 1);
    
    for k = 1:num_slices
        ref_slice = ref_norm(:, :, k);
        eval_slice = eval_norm(:, :, k);
        
        % Skip empty slices
        if max(ref_slice(:)) < 1e-10 && max(eval_slice(:)) < 1e-10
            ssim_per_slice(k) = 1.0;  % Both empty = perfect match
            continue;
        end
        
        % Use MATLAB's built-in ssim if available, else manual computation
        if exist('ssim', 'file')
            ssim_per_slice(k) = ssim(eval_slice, ref_slice, 'DynamicRange', 1.0);
        else
            ssim_per_slice(k) = compute_ssim_manual(ref_slice, eval_slice);
        end
    end
    
    ssim_results.ssim_per_slice = ssim_per_slice;
    ssim_results.ssim_mean_slice = mean(ssim_per_slice);
    
    % --- Overall 3D SSIM ---
    % Computed as mean of per-slice SSIM weighted by slice dose content
    slice_weights = zeros(num_slices, 1);
    for k = 1:num_slices
        slice_weights(k) = max(max(reference(:,:,k)));
    end
    
    if sum(slice_weights) > 0
        slice_weights = slice_weights / sum(slice_weights);
        ssim_results.ssim_3d = sum(ssim_per_slice .* slice_weights);
    else
        ssim_results.ssim_3d = mean(ssim_per_slice);
    end
    
    ssim_results.dynamic_range = dynamic_range;
    ssim_results.num_slices = num_slices;
end


function val = compute_ssim_manual(ref, eval_img)
    % Manual SSIM computation for a 2D image pair (fallback if Image Processing
    % Toolbox ssim function is not available)
    %
    % Uses standard SSIM formula with default constants:
    %   C1 = (K1*L)^2, C2 = (K2*L)^2 where K1=0.01, K2=0.03, L=1.0
    
    K1 = 0.01; K2 = 0.03; L = 1.0;
    C1 = (K1 * L)^2;
    C2 = (K2 * L)^2;
    
    % Gaussian window (11x11, sigma=1.5)
    win_size = 11;
    sigma = 1.5;
    [x, y] = meshgrid(-(win_size-1)/2:(win_size-1)/2);
    g = exp(-(x.^2 + y.^2) / (2 * sigma^2));
    g = g / sum(g(:));
    
    % Local means
    mu_ref = conv2(ref, g, 'valid');
    mu_eval = conv2(eval_img, g, 'valid');
    
    mu_ref_sq = mu_ref .^ 2;
    mu_eval_sq = mu_eval .^ 2;
    mu_ref_eval = mu_ref .* mu_eval;
    
    % Local variances and covariance
    sigma_ref_sq = conv2(ref.^2, g, 'valid') - mu_ref_sq;
    sigma_eval_sq = conv2(eval_img.^2, g, 'valid') - mu_eval_sq;
    sigma_ref_eval = conv2(ref .* eval_img, g, 'valid') - mu_ref_eval;
    
    % SSIM map
    numerator = (2 * mu_ref_eval + C1) .* (2 * sigma_ref_eval + C2);
    denominator = (mu_ref_sq + mu_eval_sq + C1) .* (sigma_ref_sq + sigma_eval_sq + C2);
    ssim_map = numerator ./ denominator;
    
    val = mean(ssim_map(:));
end


%% ======================= VISUALIZATION FUNCTION ==========================

function plot_analysis_results(reference, evaluated, gamma, ssim_results, spacing, comparison_name, output_dir, config)
%% PLOT_ANALYSIS_RESULTS - Generate comparison and gamma visualization figures
%
%   Creates two figures per comparison:
%     1. Dose comparison: side-by-side axial slices + dose difference
%     2. Gamma map: axial slices with gamma overlay + histogram
%
%   INPUTS:
%       reference       - 3D reference dose (Gy)
%       evaluated       - 3D evaluated dose (Gy)
%       gamma           - Gamma result struct from compute_gamma_analysis
%       ssim_results    - SSIM result struct from compute_dose_ssim
%       spacing         - [dx, dy, dz] in mm
%       comparison_name - String for figure titles/filenames (e.g., 'ETHOS_vs_RS')
%       output_dir      - Directory to save figures
%       config          - Config struct

    % Determine which slices to plot
    plot_slices = select_plot_slices(reference, config);
    num_plot = length(plot_slices);
    
    % Pretty name for titles
    title_name = strrep(comparison_name, '_', ' ');
    
    % Common dose colormap range
    dose_max = max(max(reference(:)), max(evaluated(:)));
    dose_clim = [0, dose_max];
    
    % ===== FIGURE 1: Dose Comparison =====
    fig1 = figure('Position', [100, 100, 400*num_plot, 1000], 'Visible', 'off');
    
    for s = 1:num_plot
        k = plot_slices(s);
        
        ref_slice = reference(:, :, k);
        eval_slice = evaluated(:, :, k);
        diff_slice = evaluated(:, :, k) - reference(:, :, k);
        
        % Reference dose
        subplot(3, num_plot, s);
        imagesc(ref_slice, dose_clim);
        colormap(gca, 'jet');
        colorbar;
        title(sprintf('Reference (slice %d)', k));
        axis image; axis off;
        
        % Evaluated dose
        subplot(3, num_plot, num_plot + s);
        imagesc(eval_slice, dose_clim);
        colormap(gca, 'jet');
        colorbar;
        title(sprintf('Evaluated (slice %d)', k));
        axis image; axis off;
        
        % Dose difference
        subplot(3, num_plot, 2*num_plot + s);
        diff_lim = max(abs(diff_slice(:)));
        if diff_lim < 1e-10, diff_lim = 1; end
        imagesc(diff_slice, [-diff_lim, diff_lim]);
        colormap(gca, rdbu_colormap());
        colorbar;
        title(sprintf('Difference (slice %d)', k));
        axis image; axis off;
    end
    
    sgtitle(sprintf('Dose Comparison: %s', title_name), 'FontSize', 14, 'FontWeight', 'bold');
    
    % Save
    dose_fig_path = fullfile(output_dir, sprintf('dose_comparison_%s.png', lower(comparison_name)));
    exportgraphics(fig1, dose_fig_path, 'Resolution', 150);
    close(fig1);
    fprintf('      Saved: %s\n', dose_fig_path);
    
    % ===== FIGURE 2: Gamma Analysis =====
    fig2 = figure('Position', [100, 100, 400*num_plot + 400, 800], 'Visible', 'off');
    
    gamma_map = gamma.gamma_map;
    gamma_clim = [0, 2];
    
    for s = 1:num_plot
        k = plot_slices(s);
        
        gamma_slice = gamma_map(:, :, k);
        ref_slice = reference(:, :, k);
        
        % Gamma map
        subplot(2, num_plot + 1, s);
        imagesc(gamma_slice, gamma_clim);
        colormap(gca, gamma_colormap());
        colorbar;
        title(sprintf('Gamma (slice %d)', k));
        axis image; axis off;
        
        % Pass/fail map
        subplot(2, num_plot + 1, num_plot + 1 + s);
        pass_fail = nan(size(gamma_slice));
        valid = ~isnan(gamma_slice);
        pass_fail(valid & gamma_slice <= 1) = 1;  % Pass
        pass_fail(valid & gamma_slice > 1) = 0;    % Fail
        imagesc(pass_fail, [0, 1]);
        colormap(gca, [1 0 0; 0 0.7 0]);  % Red=fail, Green=pass
        colorbar('Ticks', [0.25, 0.75], 'TickLabels', {'Fail', 'Pass'});
        title(sprintf('Pass/Fail (slice %d)', k));
        axis image; axis off;
    end
    
    % Gamma histogram in last column
    subplot(2, num_plot + 1, [num_plot + 1, 2*(num_plot + 1)]);
    valid_gammas = gamma_map(~isnan(gamma_map));
    histogram(valid_gammas, 0:0.05:3, 'FaceColor', [0.3 0.5 0.8], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.8);
    hold on;
    xline(1, 'r--', 'LineWidth', 2);
    hold off;
    xlabel('Gamma Index');
    ylabel('Voxel Count');
    title(sprintf('Gamma Histogram\nPass: %.1f%% | Mean: %.2f', ...
        gamma.pass_rate, gamma.mean_gamma));
    xlim([0, 3]);
    grid on;
    
    % Add SSIM info to title
    sgtitle(sprintf('Gamma Analysis: %s (%.0f%%/%.0fmm) | SSIM: %.4f', ...
        title_name, gamma.parameters.dose_pct, gamma.parameters.dist_mm, ...
        ssim_results.ssim_3d), 'FontSize', 14, 'FontWeight', 'bold');
    
    % Save
    gamma_fig_path = fullfile(output_dir, sprintf('gamma_%s.png', lower(comparison_name)));
    exportgraphics(fig2, gamma_fig_path, 'Resolution', 150);
    close(fig2);
    fprintf('      Saved: %s\n', gamma_fig_path);
end


function slices = select_plot_slices(dose, config)
    % Select representative axial slices for plotting
    %
    % If config.analysis_plot_slices is numeric, use directly.
    % If 'auto', pick 3 slices: 25th, 50th, 75th percentile of dose-containing slices.
    
    if isfield(config, 'analysis_plot_slices') && isnumeric(config.analysis_plot_slices)
        slices = config.analysis_plot_slices;
        % Clamp to valid range
        slices = slices(slices >= 1 & slices <= size(dose, 3));
        if isempty(slices)
            slices = round(size(dose, 3) / 2);
        end
        return;
    end
    
    % Auto mode: find slices with significant dose
    num_slices = size(dose, 3);
    slice_max = zeros(num_slices, 1);
    for k = 1:num_slices
        slice_max(k) = max(max(dose(:, :, k)));
    end
    
    dose_threshold = 0.1 * max(slice_max);  % 10% of peak
    active_slices = find(slice_max >= dose_threshold);
    
    if isempty(active_slices)
        slices = round(num_slices / 2);
    elseif length(active_slices) <= 3
        slices = active_slices(:)';
    else
        % Pick 25th, 50th, 75th percentile of active range
        n = length(active_slices);
        idx_25 = active_slices(max(1, round(0.25 * n)));
        idx_50 = active_slices(max(1, round(0.50 * n)));
        idx_75 = active_slices(max(1, round(0.75 * n)));
        slices = unique([idx_25, idx_50, idx_75]);
    end
end


function cmap = rdbu_colormap()
    % Red-white-blue diverging colormap for dose difference display
    n = 256;
    half = n / 2;
    
    % Blue to white
    r1 = linspace(0.1, 1, half)';
    g1 = linspace(0.2, 1, half)';
    b1 = linspace(0.7, 1, half)';
    
    % White to red
    r2 = linspace(1, 0.8, half)';
    g2 = linspace(1, 0.1, half)';
    b2 = linspace(1, 0.1, half)';
    
    cmap = [r1, g1, b1; r2, g2, b2];
end


function cmap = gamma_colormap()
    % Custom colormap for gamma display: green (pass) -> yellow -> red (fail)
    n = 256;
    
    % Green to yellow (gamma 0 to 1)
    half = round(n / 2);
    r1 = linspace(0, 1, half)';
    g1 = ones(half, 1) * 0.8;
    b1 = linspace(0.2, 0, half)';
    
    % Yellow to red (gamma 1 to 2)
    remainder = n - half;
    r2 = ones(remainder, 1);
    g2 = linspace(0.8, 0, remainder)';
    b2 = zeros(remainder, 1);
    
    cmap = [r1, g1, b1; r2, g2, b2];
end
