%% =========================================================================
%  STEP2_KWAVE_SIMULATION.m
%  k-Wave Photoacoustic Dose Reconstruction for ETHOS Pipeline
%  =========================================================================
%
%  PURPOSE:
%    Convert radiation dose to acoustic initial pressure using tissue-specific
%    Gruneisen coefficients, simulate photoacoustic wave propagation with 
%    k-Wave, and reconstruct dose via time reversal.
%
%  PHYSICS:
%    Initial pressure: p₀(r) = D(r) × Γ(r) × ρ(r)
%    Reconstruction:   D_recon = p₀_recon / (Γ × ρ)
%
%  FUNCTIONS PROVIDED:
%    - step2_run_all_simulations: Master function to process all fields
%    - create_acoustic_medium: Build k-Wave medium from CT data
%    - run_single_field_simulation: Forward + time-reversal for one field
%    - assign_tissue_properties: Map HU to acoustic properties
%
%  AUTHOR: ETHOS Pipeline Team
%  DATE: February 2026
%  VERSION: 2.0 (Modularized for master pipeline integration)
%
%  See also: kWaveGrid, kspaceFirstOrder3D, step15_process_doses
%  =========================================================================

%% ========================= MAIN ENTRY POINT ==============================

function [total_recon, field_recons, sim_log] = step2_run_all_simulations(patient_id, session, config)
%STEP2_RUN_ALL_SIMULATIONS Run k-Wave simulations for all fields
%
%   [total_recon, field_recons, sim_log] = step2_run_all_simulations(patient_id, session, config)
%
%   INPUTS:
%       patient_id  - String, patient identifier
%       session     - String, session name
%       config      - Configuration struct with:
%           .working_dir         - Base directory path
%           .gruneisen_method    - 'uniform' | 'threshold_1' | 'threshold_2'
%           .tissue_tables       - Lookup tables (from define_tissue_tables)
%           .dose_per_pulse_cGy  - Dose per LINAC pulse (default: 0.16)
%           .pml_size            - PML thickness in voxels (default: 10)
%           .cfl_number          - CFL stability number (default: 0.3)
%           .use_gpu             - Use GPU acceleration (default: true)
%           .use_parallel        - Use parfor for fields (default: false)
%           .plot_sim            - Show k-Wave plots (default: false)
%
%   OUTPUTS:
%       total_recon  - 3D array of total reconstructed dose (sum of all fields)
%       field_recons - Cell array of individual field reconstructions
%       sim_log      - Struct with timing and diagnostic information

    %% Input validation and defaults
    config = validate_and_set_defaults(config);
    
    fprintf('\n=========================================================\n');
    fprintf('  Step 2: k-Wave Photoacoustic Simulation\n');
    fprintf('=========================================================\n');
    fprintf('  Patient: %s, Session: %s\n', patient_id, session);
    fprintf('  Gruneisen method: %s\n', config.gruneisen_method);
    fprintf('  GPU enabled: %s\n', mat2str(config.use_gpu));
    fprintf('=========================================================\n\n');
    
    %% Initialize simulation log
    sim_log = struct();
    sim_log.patient_id = patient_id;
    sim_log.session = session;
    sim_log.config = config;
    sim_log.start_time = datetime('now');
    
    %% Load processed data
    fprintf('[Step 2.1] Loading processed field doses and CT...\n');
    
    processed_dir = fullfile(config.working_dir, 'RayStationFiles', ...
        patient_id, session, 'processed');
    
    if ~exist(processed_dir, 'dir')
        error('step2:DataNotFound', ...
            'Processed directory not found: %s\nRun step15_process_doses first.', processed_dir);
    end
    
    % Load field doses from individual files
    field_dose_files = dir(fullfile(processed_dir, 'field_dose_*.mat'));
    num_total_fields = length(field_dose_files);
    
    if num_total_fields == 0
        error('step2:NoFieldDoses', ...
            'No field dose files found in: %s', processed_dir);
    end
    
    fprintf('  Found %d field dose files\n', num_total_fields);
    
    % Load SCT resampled
    sct_file = fullfile(processed_dir, 'sct_resampled.mat');
    if ~exist(sct_file, 'file')
        error('step2:NoSCT', 'SCT resampled file not found: %s', sct_file);
    end
    
    sct_data = load(sct_file);
    if isfield(sct_data, 'sct_resampled')
        sct_resampled = sct_data.sct_resampled;
    else
        error('step2:InvalidSCT', 'sct_resampled.mat does not contain sct_resampled struct');
    end
    
    fprintf('  SCT dimensions: [%d × %d × %d]\n', ...
        sct_resampled.dimensions(1), sct_resampled.dimensions(2), sct_resampled.dimensions(3));
    
    %% Create acoustic medium
    fprintf('\n[Step 2.2] Creating acoustic medium from CT...\n');
    
    medium = create_acoustic_medium(sct_resampled, config);
    
    fprintf('  Density range: [%.0f, %.0f] kg/m³\n', ...
        min(medium.density(:)), max(medium.density(:)));
    fprintf('  Sound speed range: [%.0f, %.0f] m/s\n', ...
        min(medium.sound_speed(:)), max(medium.sound_speed(:)));
    
    sim_log.medium_stats.density_range = [min(medium.density(:)), max(medium.density(:))];
    sim_log.medium_stats.sound_speed_range = [min(medium.sound_speed(:)), max(medium.sound_speed(:))];
    
    %% Create output directory
    sim_dir = fullfile(config.working_dir, 'SimulationResults', ...
        patient_id, session, config.gruneisen_method);
    
    if ~exist(sim_dir, 'dir')
        mkdir(sim_dir);
        fprintf('  Created output directory: %s\n', sim_dir);
    end
    
    %% Process each field
    fprintf('\n[Step 2.3] Running k-Wave simulations for %d fields...\n', num_total_fields);
    
    grid_dims = sct_resampled.dimensions;
    total_recon = zeros(grid_dims);
    field_recons = cell(num_total_fields, 1);
    field_times = zeros(num_total_fields, 1);
    
    total_sim_start = tic;
    
    if config.use_parallel && num_total_fields > 1
        % Parallel processing with parfor
        fprintf('  Using parallel processing...\n');
        
        % Pre-load all field doses (memory permitting)
        field_doses = cell(num_total_fields, 1);
        for i = 1:num_total_fields
            data = load(fullfile(processed_dir, field_dose_files(i).name));
            field_doses{i} = data.field_dose;
        end
        
        recon_results = cell(num_total_fields, 1);
        parfor_times = zeros(num_total_fields, 1);
        
        parfor f_idx = 1:num_total_fields
            field_start = tic;
            
            field_dose = field_doses{f_idx};
            
            fprintf('    Field %d: gantry=%.1f°, processing...\n', ...
                f_idx, field_dose.gantry_angle);
            
            try
                recon_results{f_idx} = run_single_field_simulation(...
                    field_dose, sct_resampled, medium, config);
                parfor_times(f_idx) = toc(field_start);
                fprintf('    Field %d: complete (%.1f sec)\n', f_idx, parfor_times(f_idx));
            catch ME
                fprintf('    Field %d: FAILED - %s\n', f_idx, ME.message);
                recon_results{f_idx} = zeros(grid_dims);
                parfor_times(f_idx) = toc(field_start);
            end
        end
        
        % Aggregate results
        for f_idx = 1:num_total_fields
            field_recons{f_idx} = recon_results{f_idx};
            total_recon = total_recon + recon_results{f_idx};
            field_times(f_idx) = parfor_times(f_idx);
            
            % Save individual field reconstruction
            recon_dose = recon_results{f_idx};
            save(fullfile(sim_dir, sprintf('field_recon_%03d.mat', f_idx)), ...
                'recon_dose', '-v7.3');
        end
        
    else
        % Serial processing
        fprintf('  Using serial processing...\n');
        
        for f_idx = 1:num_total_fields
            field_start = tic;
            
            % Load this field's dose
            data = load(fullfile(processed_dir, field_dose_files(f_idx).name));
            field_dose = data.field_dose;
            
            fprintf('\n  --- Field %d/%d ---\n', f_idx, num_total_fields);
            fprintf('    Gantry angle: %.1f°\n', field_dose.gantry_angle);
            fprintf('    Max dose: %.4f Gy\n', max(field_dose.dose_Gy(:)));
            
            try
                recon_dose = run_single_field_simulation(...
                    field_dose, sct_resampled, medium, config);
                
                field_recons{f_idx} = recon_dose;
                total_recon = total_recon + recon_dose;
                
                % Save individual field reconstruction
                save(fullfile(sim_dir, sprintf('field_recon_%03d.mat', f_idx)), ...
                    'recon_dose', '-v7.3');
                
                field_times(f_idx) = toc(field_start);
                fprintf('    Reconstructed max: %.4f Gy (%.1f sec)\n', ...
                    max(recon_dose(:)), field_times(f_idx));
                
            catch ME
                fprintf('    ERROR: %s\n', ME.message);
                field_recons{f_idx} = zeros(grid_dims);
                field_times(f_idx) = toc(field_start);
                
                sim_log.errors{end+1} = struct('field', f_idx, 'message', ME.message);
            end
        end
    end
    
    total_sim_time = toc(total_sim_start);
    
    %% Save total reconstruction
    fprintf('\n[Step 2.4] Saving results...\n');
    
    metadata = struct();
    metadata.patient_id = patient_id;
    metadata.session = session;
    metadata.gruneisen_method = config.gruneisen_method;
    metadata.dimensions = grid_dims;
    metadata.spacing = sct_resampled.spacing;
    metadata.origin = sct_resampled.origin;
    metadata.num_fields = num_total_fields;
    metadata.timestamp = datetime('now');
    
    save(fullfile(sim_dir, 'total_recon_dose.mat'), ...
        'total_recon', 'metadata', '-v7.3');
    
    %% Finalize log
    sim_log.end_time = datetime('now');
    sim_log.total_time_sec = total_sim_time;
    sim_log.field_times_sec = field_times;
    sim_log.total_recon_max_Gy = max(total_recon(:));
    sim_log.num_fields_processed = num_total_fields;
    
    save(fullfile(sim_dir, 'simulation_log.mat'), 'sim_log', '-v7.3');
    
    %% Summary
    fprintf('\n=========================================================\n');
    fprintf('  Step 2 Complete\n');
    fprintf('=========================================================\n');
    fprintf('  Fields processed: %d\n', num_total_fields);
    fprintf('  Total simulation time: %.1f sec (%.1f sec/field avg)\n', ...
        total_sim_time, total_sim_time / num_total_fields);
    fprintf('  Max reconstructed dose: %.4f Gy\n', max(total_recon(:)));
    fprintf('  Output directory: %s\n', sim_dir);
    fprintf('=========================================================\n\n');
end


%% ========================= ACOUSTIC MEDIUM ===============================

function medium = create_acoustic_medium(sct_resampled, config)
%CREATE_ACOUSTIC_MEDIUM Build k-Wave medium structure from CT data
%
%   medium = create_acoustic_medium(sct_resampled, config)
%
%   INPUTS:
%       sct_resampled - Struct with:
%           .cubeHU      - 3D HU array
%           .dimensions  - [nx, ny, nz]
%           .spacing     - [dx, dy, dz] in mm
%       config - Configuration struct with:
%           .gruneisen_method - 'uniform' | 'threshold_1' | 'threshold_2'
%           .tissue_tables    - Lookup tables
%
%   OUTPUTS:
%       medium - k-Wave medium struct with:
%           .density      - 3D density array (kg/m³)
%           .sound_speed  - 3D sound speed array (m/s)
%           .alpha_coeff  - 3D absorption coefficient array
%           .alpha_power  - Scalar absorption power law exponent
%           .gruneisen    - 3D Gruneisen coefficient array

    % Extract HU cube
    if isfield(sct_resampled, 'cubeHU')
        ctCube = sct_resampled.cubeHU;
    elseif isfield(sct_resampled, 'cube')
        ctCube = sct_resampled.cube;
    else
        error('create_acoustic_medium:NoHU', ...
            'sct_resampled must contain cubeHU or cube field');
    end
    
    grid_dims = size(ctCube);
    
    fprintf('  Assigning tissue properties using method: %s\n', config.gruneisen_method);
    fprintf('  HU range in CT: [%.0f, %.0f]\n', min(ctCube(:)), max(ctCube(:)));
    
    % Assign tissue properties based on method
    [density, sound_speed, alpha_coeff, alpha_power, gruneisen] = ...
        assign_tissue_properties(ctCube, config.gruneisen_method, config.tissue_tables);
    
    % Ensure minimum values for numerical stability
    density = max(density, 1);           % Minimum 1 kg/m³
    sound_speed = max(sound_speed, 100); % Minimum 100 m/s
    alpha_coeff = max(alpha_coeff, 0);   % Non-negative absorption
    
    % Build medium struct
    medium = struct();
    medium.density = density;
    medium.sound_speed = sound_speed;
    medium.alpha_coeff = alpha_coeff;
    medium.alpha_power = alpha_power;  % Scalar for k-Wave
    medium.gruneisen = gruneisen;
    
    % Store tissue statistics for debugging
    medium.stats = struct();
    medium.stats.density_mean = mean(density(:));
    medium.stats.sound_speed_mean = mean(sound_speed(:));
    medium.stats.gruneisen_mean = mean(gruneisen(:));
    medium.stats.gruneisen_nonzero_frac = sum(gruneisen(:) > 0) / numel(gruneisen);
end


function [density, sound_speed, alpha_coeff, alpha_power, gruneisen] = ...
    assign_tissue_properties(ctCube, method, tables)
%ASSIGN_TISSUE_PROPERTIES Map HU values to acoustic properties
%
%   [density, sound_speed, alpha_coeff, alpha_power, gruneisen] = ...
%       assign_tissue_properties(ctCube, method, tables)
%
%   INPUTS:
%       ctCube  - 3D array of Hounsfield Units
%       method  - 'uniform' | 'threshold_1' | 'threshold_2'
%       tables  - Struct with tissue property lookup tables
%
%   OUTPUTS:
%       density     - 3D density array (kg/m³)
%       sound_speed - 3D sound speed array (m/s)
%       alpha_coeff - 3D absorption coefficient array
%       alpha_power - Scalar absorption power law exponent
%       gruneisen   - 3D Gruneisen coefficient array

    grid_dims = size(ctCube);
    
    switch lower(method)
        case 'uniform'
            % Single uniform properties for entire volume
            if isfield(tables, 'uniform')
                t = tables.uniform;
            else
                % Default uniform values
                t = struct();
                t.density = 1000;
                t.sound_speed = 1540;
                t.alpha_coeff = 0.5;
                t.alpha_power = 1.1;
                t.gruneisen = 1.0;
            end
            
            density = t.density * ones(grid_dims);
            sound_speed = t.sound_speed * ones(grid_dims);
            alpha_coeff = t.alpha_coeff * ones(grid_dims);
            alpha_power = t.alpha_power;
            gruneisen = t.gruneisen * ones(grid_dims);
            
        case 'threshold_1'
            % Detailed tissue segmentation (9 tissues)
            [density, sound_speed, alpha_coeff, alpha_power, gruneisen] = ...
                apply_threshold_segmentation(ctCube, tables.threshold_1);
            
        case 'threshold_2'
            % Simplified 4-tissue model
            [density, sound_speed, alpha_coeff, alpha_power, gruneisen] = ...
                apply_threshold_segmentation(ctCube, tables.threshold_2);
            
        otherwise
            error('assign_tissue_properties:InvalidMethod', ...
                'Unknown Gruneisen method: %s. Use: uniform, threshold_1, threshold_2', method);
    end
end


function [density, sound_speed, alpha_coeff, alpha_power, gruneisen] = ...
    apply_threshold_segmentation(ctCube, table)
%APPLY_THRESHOLD_SEGMENTATION Segment CT and assign properties by HU thresholds
%
%   Uses the HU boundaries in the table to segment the CT and assign
%   tissue-specific acoustic properties to each region.

    grid_dims = size(ctCube);
    
    % Initialize output arrays
    density = zeros(grid_dims);
    sound_speed = zeros(grid_dims);
    alpha_coeff = zeros(grid_dims);
    gruneisen = zeros(grid_dims);
    alpha_power_values = [];
    
    % Get number of tissue types
    num_tissues = length(table.tissue_names);
    hu_bounds = table.hu_boundaries;
    
    fprintf('  Segmenting %d tissue types:\n', num_tissues);
    
    % Process each tissue type
    for t = 1:num_tissues
        hu_min = hu_bounds(t);
        hu_max = hu_bounds(t + 1);
        
        % Create mask for this tissue type
        if t == 1
            mask = ctCube < hu_max;
        elseif t == num_tissues
            mask = ctCube >= hu_min;
        else
            mask = (ctCube >= hu_min) & (ctCube < hu_max);
        end
        
        num_voxels = sum(mask(:));
        pct = 100 * num_voxels / numel(ctCube);
        
        if num_voxels > 0
            fprintf('    %s: %d voxels (%.1f%%)\n', table.tissue_names{t}, num_voxels, pct);
            
            % Assign properties
            density(mask) = table.density(t);
            sound_speed(mask) = table.sound_speed(t);
            alpha_coeff(mask) = table.alpha_coeff(t);
            gruneisen(mask) = table.gruneisen(t);
            alpha_power_values = [alpha_power_values; repmat(table.alpha_power(t), num_voxels, 1)];
        end
    end
    
    % Handle any unassigned voxels (shouldn't happen with proper boundaries)
    unassigned = (density == 0);
    if any(unassigned(:))
        num_unassigned = sum(unassigned(:));
        fprintf('    WARNING: %d unassigned voxels (setting to water)\n', num_unassigned);
        density(unassigned) = 1000;
        sound_speed(unassigned) = 1480;
        alpha_coeff(unassigned) = 0.002;
        gruneisen(unassigned) = 0.11;
    end
    
    % k-Wave requires scalar alpha_power - use weighted mean
    if ~isempty(alpha_power_values)
        alpha_power = mean(alpha_power_values);
    else
        alpha_power = 1.1;  % Default
    end
    
    fprintf('    Mean alpha_power: %.2f\n', alpha_power);
end


%% ========================= SINGLE FIELD SIMULATION =======================

function recon_dose = run_single_field_simulation(field_dose, sct_resampled, medium, config)
%RUN_SINGLE_FIELD_SIMULATION Run forward and time-reversal simulation for one field
%
%   recon_dose = run_single_field_simulation(field_dose, sct_resampled, medium, config)
%
%   INPUTS:
%       field_dose    - Struct with:
%           .dose_Gy      - 3D dose array (Gy)
%           .gantry_angle - Gantry angle (degrees)
%           .meterset     - Monitor units (optional)
%       sct_resampled - Struct with CT data and geometry
%       medium        - k-Wave medium struct from create_acoustic_medium
%       config        - Configuration struct
%
%   OUTPUTS:
%       recon_dose    - 3D reconstructed dose array (Gy)
%
%   ALGORITHM:
%       1. Calculate number of pulses from meterset
%       2. Convert dose to initial pressure: p0 = D/pulses × Γ × ρ
%       3. Setup k-Wave grid and time stepping
%       4. Place planar sensor inferior to dose region
%       5. Run forward simulation
%       6. Run time-reversal reconstruction
%       7. Convert pressure back to dose

    %% Extract parameters
    dose_Gy = field_dose.dose_Gy;
    grid_dims = size(dose_Gy);
    
    % Get spacing (convert mm to m for k-Wave)
    spacing_mm = sct_resampled.spacing;
    dx = spacing_mm(1) / 1000;  % m
    dy = spacing_mm(2) / 1000;  % m
    dz = spacing_mm(3) / 1000;  % m
    
    % Get meterset (default to 100 MU if not available)
    if isfield(field_dose, 'meterset') && ~isempty(field_dose.meterset)
        meterset = field_dose.meterset;
    else
        meterset = 100;  % Default assumption
    end
    
    %% Calculate number of pulses
    num_pulses = ceil(meterset / config.dose_per_pulse_cGy);
    num_pulses = max(num_pulses, 1);  % At least 1 pulse
    
    %% Convert dose to initial pressure per pulse
    % p0 = (D / num_pulses) × Γ × ρ
    dose_per_pulse = dose_Gy / num_pulses;
    
    % Get tissue-specific Gruneisen and density
    gruneisen = medium.gruneisen;
    density = medium.density;
    
    % Initial pressure
    p0 = dose_per_pulse .* gruneisen .* density;
    
    fprintf('    Pulses: %d, Max p0: %.2e Pa\n', num_pulses, max(p0(:)));
    
    %% Check if dose is significant
    if max(p0(:)) < 1e-10
        fprintf('    WARNING: No significant initial pressure, returning zeros\n');
        recon_dose = zeros(grid_dims);
        return;
    end
    
    %% Setup k-Wave grid
    Nx = grid_dims(1);
    Ny = grid_dims(2);
    Nz = grid_dims(3);
    
    kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
    
    %% Calculate time stepping
    max_sound_speed = max(medium.sound_speed(:));
    min_spacing = min([dx, dy, dz]);
    
    % CFL-based time step
    dt = config.cfl_number * min_spacing / max_sound_speed;
    
    % Simulation time (allow wave to traverse grid ~2.5 times)
    grid_diagonal = sqrt((Nx*dx)^2 + (Ny*dy)^2 + (Nz*dz)^2);
    min_sound_speed = max(min(medium.sound_speed(:)), 100);  % Avoid div by zero
    sim_time = 2.5 * grid_diagonal / min_sound_speed;
    Nt = ceil(sim_time / dt);
    
    kgrid.dt = dt;
    kgrid.Nt = Nt;
    
    fprintf('    Grid: [%d×%d×%d], dt=%.2e s, Nt=%d\n', Nx, Ny, Nz, dt, Nt);
    
    %% Setup k-Wave medium (already have density and sound_speed)
    kwave_medium = struct();
    kwave_medium.density = medium.density;
    kwave_medium.sound_speed = medium.sound_speed;
    kwave_medium.alpha_coeff = medium.alpha_coeff;
    kwave_medium.alpha_power = medium.alpha_power;
    
    %% Place sensor based on dose distribution
    % Find dose extent
    dose_threshold = 0.01 * max(dose_Gy(:));  % 1% of max dose
    dose_mask = dose_Gy > dose_threshold;
    
    if ~any(dose_mask(:))
        fprintf('    WARNING: No significant dose found\n');
        recon_dose = zeros(grid_dims);
        return;
    end
    
    % Find extent in each dimension
    [xIdx, yIdx, zIdx] = ind2sub(grid_dims, find(dose_mask));
    
    min_x = max(1, min(xIdx) - 5);
    max_x = min(Nx, max(xIdx) + 5);
    min_y = min(yIdx);
    max_y = max(yIdx);
    min_z = max(1, min(zIdx) - 5);
    max_z = min(Nz, max(zIdx) + 5);
    
    % Place planar sensor inferior to dose (in Y direction)
    sensor_y = max(1, min_y - 10);  % 10 voxels below dose
    
    % Create sensor mask
    sensor = struct();
    sensor.mask = zeros(Nx, Ny, Nz);
    sensor.mask(min_x:max_x, sensor_y, min_z:max_z) = 1;
    
    num_sensor_points = sum(sensor.mask(:));
    fprintf('    Sensor: Y=%d, points=%d\n', sensor_y, num_sensor_points);
    
    %% Setup simulation options
    if config.use_gpu
        try
            gpuDevice;
            data_cast = 'gpuArray-single';
        catch
            fprintf('    GPU not available, using CPU\n');
            data_cast = 'single';
        end
    else
        data_cast = 'single';
    end
    
    input_args = {...
        'Smooth', false, ...
        'PMLInside', false, ...
        'PMLSize', config.pml_size, ...
        'DataCast', data_cast, ...
        'PlotSim', config.plot_sim};
    
    %% Forward simulation
    fprintf('    Running forward simulation...\n');
    
    source = struct();
    source.p0 = p0;
    
    try
        sensor_data = kspaceFirstOrder3D(kgrid, kwave_medium, source, sensor, input_args{:});
        fprintf('    Forward complete. Sensor data: [%d × %d]\n', ...
            size(sensor_data, 1), size(sensor_data, 2));
    catch ME
        fprintf('    Forward simulation FAILED: %s\n', ME.message);
        recon_dose = zeros(grid_dims);
        return;
    end
    
    %% Time-reversal reconstruction
    fprintf('    Running time-reversal reconstruction...\n');
    
    % Clear p0 from source, setup time reversal
    source = rmfield(source, 'p0');
    source.p_mask = sensor.mask;
    source.p = fliplr(sensor_data);  % Time-reversed data
    source.p_mode = 'dirichlet';
    
    % Record final pressure
    sensor_tr = struct();
    sensor_tr.record = {'p_final'};
    
    try
        tr_result = kspaceFirstOrder3D(kgrid, kwave_medium, source, sensor_tr, input_args{:});
        
        % Extract reconstructed pressure
        if isstruct(tr_result) && isfield(tr_result, 'p_final')
            p0_recon = tr_result.p_final;
        else
            p0_recon = tr_result;
        end
        
        % Ensure correct dimensions
        p0_recon = reshape(p0_recon, grid_dims);
        
        % Apply positivity constraint
        p0_recon = max(p0_recon, 0);
        
        fprintf('    TR complete. Reconstructed p0: [%.2e, %.2e] Pa\n', ...
            min(p0_recon(:)), max(p0_recon(:)));
        
    catch ME
        fprintf('    Time-reversal FAILED: %s\n', ME.message);
        recon_dose = zeros(grid_dims);
        return;
    end
    
    %% Convert pressure back to dose
    % D_recon = p0_recon / (Γ × ρ) × num_pulses
    
    conversion_factor = gruneisen .* density;
    
    % Avoid division by zero
    conversion_factor(conversion_factor < 1e-10) = 1;
    
    recon_dose_per_pulse = p0_recon ./ conversion_factor;
    recon_dose = recon_dose_per_pulse * num_pulses;
    
    fprintf('    Dose reconstruction complete: max=%.4f Gy\n', max(recon_dose(:)));
end


%% ========================= HELPER FUNCTIONS ==============================

function config = validate_and_set_defaults(config)
%VALIDATE_AND_SET_DEFAULTS Validate config struct and set default values

    % Required fields
    if ~isfield(config, 'working_dir')
        error('validate_config:MissingField', 'config.working_dir is required');
    end
    
    % Set defaults
    if ~isfield(config, 'gruneisen_method')
        config.gruneisen_method = 'threshold_2';
    end
    
    if ~isfield(config, 'dose_per_pulse_cGy')
        config.dose_per_pulse_cGy = 0.16;
    end
    
    if ~isfield(config, 'pml_size')
        config.pml_size = 10;
    end
    
    if ~isfield(config, 'cfl_number')
        config.cfl_number = 0.3;
    end
    
    if ~isfield(config, 'use_gpu')
        config.use_gpu = true;
    end
    
    if ~isfield(config, 'use_parallel')
        config.use_parallel = false;
    end
    
    if ~isfield(config, 'plot_sim')
        config.plot_sim = false;
    end
    
    % Ensure tissue tables are defined
    if ~isfield(config, 'tissue_tables')
        config.tissue_tables = define_tissue_tables();
    end
end


function tables = define_tissue_tables()
%DEFINE_TISSUE_TABLES Create tissue property lookup tables
%
%   tables = define_tissue_tables()
%
%   Returns struct with lookup tables for:
%     - uniform: Single property values for entire volume
%     - threshold_1: Detailed 9-tissue segmentation
%     - threshold_2: Simplified 4-tissue model

    tables = struct();
    
    %% UNIFORM: Single tissue properties
    tables.uniform = struct();
    tables.uniform.density = 1000;        % kg/m³
    tables.uniform.sound_speed = 1540;    % m/s
    tables.uniform.alpha_coeff = 0.5;     % dB/MHz^y/cm
    tables.uniform.alpha_power = 1.1;     
    tables.uniform.gruneisen = 1.0;
    
    %% THRESHOLD_1: Detailed tissue segmentation (9 tissues)
    tables.threshold_1 = struct();
    tables.threshold_1.hu_boundaries = [-1024, -900, -500, -200, -50, 13, 50, 80, 300, 3000, Inf];
    tables.threshold_1.tissue_names  = {'Air', 'Lung', 'Fat', 'Water', 'Blood', 'Muscle', 'SoftTissue', 'Bone', 'Metal'};
    tables.threshold_1.density       = [1.2,   400,   920,   1000,  1060,   1050,    1040,        1900,  7800];
    tables.threshold_1.sound_speed   = [343,   600,   1450,  1480,  1575,   1580,    1540,        3200,  5900];
    tables.threshold_1.alpha_coeff   = [0,     0.5,   0.48,  0.0022, 0.2,   0.5,     0.5,         10,    0];
    tables.threshold_1.alpha_power   = [1.0,   1.5,   1.5,   2.0,   1.3,    1.0,     1.1,         1.0,   1.0];
    tables.threshold_1.gruneisen     = [0,     0.5,   0.7,   0.11,  0.15,   0.2,     1.0,         0,     0];
    
    %% THRESHOLD_2: Simplified 4-tissue model
    % HU Range     | ρ (kg/m³) | c (m/s) | α      | Γ
    % < -200       | 1000      | 1480    | 0.0022 | 0.11
    % -200 to -50  | 920       | 1450    | 0.48   | 0.7
    % -50 to 100   | 1040      | 1540    | 0.5    | 1.0
    % ≥ 100        | 1900      | 3200    | 10     | 0
    
    tables.threshold_2 = struct();
    tables.threshold_2.hu_boundaries = [-1024, -200, -50, 100, Inf];
    tables.threshold_2.tissue_names  = {'Water', 'Fat', 'SoftTissue', 'Bone'};
    tables.threshold_2.density       = [1000,    920,  1040,          1900];
    tables.threshold_2.sound_speed   = [1480,    1450, 1540,          3200];
    tables.threshold_2.alpha_coeff   = [0.0022,  0.48, 0.5,           10];
    tables.threshold_2.alpha_power   = [2.0,     1.5,  1.1,           1.0];
    tables.threshold_2.gruneisen     = [0.11,    0.7,  1.0,           0];
end


%% ========================= UTILITY FUNCTIONS =============================

function save_field_reconstruction(recon_dose, field_idx, patient_id, session, config)
%SAVE_FIELD_RECONSTRUCTION Save individual field reconstruction to file
    sim_dir = fullfile(config.working_dir, 'SimulationResults', ...
        patient_id, session, config.gruneisen_method);
    
    if ~exist(sim_dir, 'dir')
        mkdir(sim_dir);
    end
    
    filename = sprintf('field_recon_%03d.mat', field_idx);
    save(fullfile(sim_dir, filename), 'recon_dose', '-v7.3');
end


function [total_recon, metadata] = load_total_reconstruction(patient_id, session, config)
%LOAD_TOTAL_RECONSTRUCTION Load previously computed total reconstruction
    sim_dir = fullfile(config.working_dir, 'SimulationResults', ...
        patient_id, session, config.gruneisen_method);
    
    filepath = fullfile(sim_dir, 'total_recon_dose.mat');
    
    if ~exist(filepath, 'file')
        error('load_total_reconstruction:FileNotFound', ...
            'Total reconstruction file not found: %s', filepath);
    end
    
    data = load(filepath);
    total_recon = data.total_recon;
    metadata = data.metadata;
end


function recon_dose = load_field_reconstruction(field_idx, patient_id, session, config)
%LOAD_FIELD_RECONSTRUCTION Load individual field reconstruction from file
    sim_dir = fullfile(config.working_dir, 'SimulationResults', ...
        patient_id, session, config.gruneisen_method);
    
    filename = sprintf('field_recon_%03d.mat', field_idx);
    filepath = fullfile(sim_dir, filename);
    
    if ~exist(filepath, 'file')
        error('load_field_reconstruction:FileNotFound', ...
            'Field reconstruction file not found: %s', filepath);
    end
    
    data = load(filepath);
    recon_dose = data.recon_dose;
end
