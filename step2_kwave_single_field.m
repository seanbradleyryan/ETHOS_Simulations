function [recon_dose, sim_info] = step2_kwave_single_field(field_dose, sct_resampled, config)
%STEP2_KWAVE_SINGLE_FIELD Run k-Wave simulation for a single radiation field
%
%   [recon_dose, sim_info] = step2_kwave_single_field(field_dose, sct_resampled, config)
%
%   This is a simplified single-function interface for running one field
%   simulation. It creates the acoustic medium internally and runs the
%   forward + time-reversal simulation.
%
%   INPUTS:
%       field_dose    - Struct with:
%           .dose_Gy      - 3D dose array (Gy)
%           .gantry_angle - Gantry angle (degrees)
%           .meterset     - Monitor units (optional, default: 100)
%           .spacing      - [dx, dy, dz] in mm (optional if in sct_resampled)
%
%       sct_resampled - Struct with:
%           .cubeHU      - 3D HU array
%           .dimensions  - [nx, ny, nz]
%           .spacing     - [dx, dy, dz] in mm
%           .origin      - [x, y, z] in mm (optional)
%
%       config        - Configuration struct with (all optional with defaults):
%           .gruneisen_method    - 'uniform' | 'threshold_1' | 'threshold_2' (default: 'threshold_2')
%           .tissue_tables       - Lookup tables (auto-generated if missing)
%           .dose_per_pulse_cGy  - Dose per pulse (default: 0.16)
%           .pml_size            - PML thickness in voxels (default: 10)
%           .cfl_number          - CFL stability (default: 0.3)
%           .use_gpu             - Use GPU (default: true)
%           .plot_sim            - Show plots (default: false)
%
%   OUTPUTS:
%       recon_dose - 3D reconstructed dose array (Gy), same size as input
%       sim_info   - Struct with simulation diagnostics
%
%   EXAMPLE:
%       % Load data
%       field = load('field_dose_001.mat'); field_dose = field.field_dose;
%       sct = load('sct_resampled.mat'); sct_resampled = sct.sct_resampled;
%       
%       % Configure
%       config.gruneisen_method = 'threshold_2';
%       config.use_gpu = true;
%       
%       % Run simulation
%       [recon_dose, info] = step2_kwave_single_field(field_dose, sct_resampled, config);
%
%   See also: step2_run_all_simulations, create_acoustic_medium

    %% Validate and set defaults
    if nargin < 3
        config = struct();
    end
    config = set_defaults(config);
    
    %% Initialize output
    sim_info = struct();
    sim_info.start_time = tic;
    
    %% Create acoustic medium
    medium = create_medium_internal(sct_resampled, config);
    sim_info.medium_created = true;
    
    %% Run simulation
    [recon_dose, forward_info, tr_info] = run_simulation_internal(...
        field_dose, sct_resampled, medium, config);
    
    %% Finalize
    sim_info.elapsed_sec = toc(sim_info.start_time);
    sim_info.forward = forward_info;
    sim_info.time_reversal = tr_info;
    sim_info.max_recon_dose = max(recon_dose(:));
    sim_info.max_orig_dose = max(field_dose.dose_Gy(:));
end


%% ========================================================================
%  INTERNAL FUNCTIONS
%% ========================================================================

function config = set_defaults(config)
    defaults = struct(...
        'gruneisen_method', 'threshold_2', ...
        'dose_per_pulse_cGy', 0.16, ...
        'pml_size', 10, ...
        'cfl_number', 0.3, ...
        'use_gpu', true, ...
        'plot_sim', false);
    
    fields = fieldnames(defaults);
    for i = 1:length(fields)
        if ~isfield(config, fields{i})
            config.(fields{i}) = defaults.(fields{i});
        end
    end
    
    if ~isfield(config, 'tissue_tables')
        config.tissue_tables = get_tissue_tables();
    end
end


function tables = get_tissue_tables()
    tables = struct();
    
    % Uniform
    tables.uniform.density = 1000;
    tables.uniform.sound_speed = 1540;
    tables.uniform.alpha_coeff = 0.5;
    tables.uniform.alpha_power = 1.1;
    tables.uniform.gruneisen = 1.0;
    
    % Threshold 1 (detailed)
    tables.threshold_1.hu_boundaries = [-1024, -900, -500, -200, -50, 13, 50, 80, 300, 3000, Inf];
    tables.threshold_1.tissue_names  = {'Air', 'Lung', 'Fat', 'Water', 'Blood', 'Muscle', 'SoftTissue', 'Bone', 'Metal'};
    tables.threshold_1.density       = [1.2,   400,   920,   1000,  1060,   1050,    1040,        1900,  7800];
    tables.threshold_1.sound_speed   = [343,   600,   1450,  1480,  1575,   1580,    1540,        3200,  5900];
    tables.threshold_1.alpha_coeff   = [0,     0.5,   0.48,  0.0022, 0.2,   0.5,     0.5,         10,    0];
    tables.threshold_1.alpha_power   = [1.0,   1.5,   1.5,   2.0,   1.3,    1.0,     1.1,         1.0,   1.0];
    tables.threshold_1.gruneisen     = [0,     0.5,   0.7,   0.11,  0.15,   0.2,     1.0,         0,     0];
    
    % Threshold 2 (simplified)
    tables.threshold_2.hu_boundaries = [-1024, -200, -50, 100, Inf];
    tables.threshold_2.tissue_names  = {'Water', 'Fat', 'SoftTissue', 'Bone'};
    tables.threshold_2.density       = [1000,    920,  1040,          1900];
    tables.threshold_2.sound_speed   = [1480,    1450, 1540,          3200];
    tables.threshold_2.alpha_coeff   = [0.0022,  0.48, 0.5,           10];
    tables.threshold_2.alpha_power   = [2.0,     1.5,  1.1,           1.0];
    tables.threshold_2.gruneisen     = [0.11,    0.7,  1.0,           0];
end


function medium = create_medium_internal(sct_resampled, config)
    % Extract HU cube
    if isfield(sct_resampled, 'cubeHU')
        ctCube = sct_resampled.cubeHU;
    else
        error('sct_resampled must contain cubeHU field');
    end
    
    grid_dims = size(ctCube);
    method = lower(config.gruneisen_method);
    tables = config.tissue_tables;
    
    if strcmp(method, 'uniform')
        t = tables.uniform;
        medium.density = t.density * ones(grid_dims);
        medium.sound_speed = t.sound_speed * ones(grid_dims);
        medium.alpha_coeff = t.alpha_coeff * ones(grid_dims);
        medium.alpha_power = t.alpha_power;
        medium.gruneisen = t.gruneisen * ones(grid_dims);
    else
        if strcmp(method, 'threshold_1')
            table = tables.threshold_1;
        else
            table = tables.threshold_2;
        end
        
        medium.density = zeros(grid_dims);
        medium.sound_speed = zeros(grid_dims);
        medium.alpha_coeff = zeros(grid_dims);
        medium.gruneisen = zeros(grid_dims);
        alpha_powers = [];
        
        num_tissues = length(table.tissue_names);
        hu_bounds = table.hu_boundaries;
        
        for t = 1:num_tissues
            if t == 1
                mask = ctCube < hu_bounds(2);
            elseif t == num_tissues
                mask = ctCube >= hu_bounds(t);
            else
                mask = (ctCube >= hu_bounds(t)) & (ctCube < hu_bounds(t+1));
            end
            
            n = sum(mask(:));
            if n > 0
                medium.density(mask) = table.density(t);
                medium.sound_speed(mask) = table.sound_speed(t);
                medium.alpha_coeff(mask) = table.alpha_coeff(t);
                medium.gruneisen(mask) = table.gruneisen(t);
                alpha_powers = [alpha_powers; repmat(table.alpha_power(t), n, 1)];
            end
        end
        
        % Handle unassigned (shouldn't happen)
        unassigned = medium.density == 0;
        if any(unassigned(:))
            medium.density(unassigned) = 1000;
            medium.sound_speed(unassigned) = 1480;
            medium.alpha_coeff(unassigned) = 0.002;
            medium.gruneisen(unassigned) = 0.11;
        end
        
        medium.alpha_power = mean(alpha_powers);
    end
    
    % Ensure stability
    medium.density = max(medium.density, 1);
    medium.sound_speed = max(medium.sound_speed, 100);
end


function [recon_dose, fwd_info, tr_info] = run_simulation_internal(field_dose, sct_resampled, medium, config)
    dose_Gy = field_dose.dose_Gy;
    grid_dims = size(dose_Gy);
    
    % Initialize info structs
    fwd_info = struct();
    tr_info = struct();
    
    % Get spacing
    spacing_mm = sct_resampled.spacing;
    dx = spacing_mm(1) / 1000;
    dy = spacing_mm(2) / 1000;
    dz = spacing_mm(3) / 1000;
    
    % Calculate pulses
    meterset = 100;
    if isfield(field_dose, 'meterset') && ~isempty(field_dose.meterset)
        meterset = field_dose.meterset;
    end
    num_pulses = max(1, ceil(meterset / config.dose_per_pulse_cGy));
    
    % Initial pressure
    dose_per_pulse = dose_Gy / num_pulses;
    p0 = dose_per_pulse .* medium.gruneisen .* medium.density;
    
    fwd_info.num_pulses = num_pulses;
    fwd_info.max_p0 = max(p0(:));
    
    % Check for significant signal
    if max(p0(:)) < 1e-10
        recon_dose = zeros(grid_dims);
        return;
    end
    
    % Setup k-Wave grid
    Nx = grid_dims(1); Ny = grid_dims(2); Nz = grid_dims(3);
    kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
    
    % Time stepping
    max_c = max(medium.sound_speed(:));
    min_c = max(min(medium.sound_speed(:)), 100);
    min_spacing = min([dx, dy, dz]);
    
    dt = config.cfl_number * min_spacing / max_c;
    grid_diag = sqrt((Nx*dx)^2 + (Ny*dy)^2 + (Nz*dz)^2);
    sim_time = 2.5 * grid_diag / min_c;
    Nt = ceil(sim_time / dt);
    
    kgrid.dt = dt;
    kgrid.Nt = Nt;
    
    fwd_info.dt = dt;
    fwd_info.Nt = Nt;
    
    % k-Wave medium
    kwave_medium.density = medium.density;
    kwave_medium.sound_speed = medium.sound_speed;
    kwave_medium.alpha_coeff = medium.alpha_coeff;
    kwave_medium.alpha_power = medium.alpha_power;
    
    % Sensor placement
    dose_thresh = 0.01 * max(dose_Gy(:));
    dose_mask = dose_Gy > dose_thresh;
    
    if ~any(dose_mask(:))
        recon_dose = zeros(grid_dims);
        return;
    end
    
    [xIdx, yIdx, zIdx] = ind2sub(grid_dims, find(dose_mask));
    min_x = max(1, min(xIdx)-5);  max_x = min(Nx, max(xIdx)+5);
    min_y = min(yIdx);
    min_z = max(1, min(zIdx)-5);  max_z = min(Nz, max(zIdx)+5);
    sensor_y = max(1, min_y - 10);
    
    sensor.mask = zeros(Nx, Ny, Nz);
    sensor.mask(min_x:max_x, sensor_y, min_z:max_z) = 1;
    fwd_info.sensor_points = sum(sensor.mask(:));
    
    % Data cast
    if config.use_gpu
        try
            gpuDevice;
            data_cast = 'gpuArray-single';
        catch
            data_cast = 'single';
        end
    else
        data_cast = 'single';
    end
    
    input_args = {'Smooth', false, 'PMLInside', false, 'PMLSize', config.pml_size, ...
                  'DataCast', data_cast, 'PlotSim', config.plot_sim};
    
    % Forward simulation
    source.p0 = p0;
    fwd_start = tic;
    
    try
        sensor_data = kspaceFirstOrder3D(kgrid, kwave_medium, source, sensor, input_args{:});
        fwd_info.elapsed_sec = toc(fwd_start);
        fwd_info.success = true;
    catch ME
        fwd_info.success = false;
        fwd_info.error = ME.message;
        recon_dose = zeros(grid_dims);
        return;
    end
    
    % Time reversal
    source = rmfield(source, 'p0');
    source.p_mask = sensor.mask;
    source.p = fliplr(sensor_data);
    source.p_mode = 'dirichlet';
    
    sensor_tr.record = {'p_final'};
    tr_start = tic;
    
    try
        tr_result = kspaceFirstOrder3D(kgrid, kwave_medium, source, sensor_tr, input_args{:});
        tr_info.elapsed_sec = toc(tr_start);
        tr_info.success = true;
        
        if isstruct(tr_result) && isfield(tr_result, 'p_final')
            p0_recon = tr_result.p_final;
        else
            p0_recon = tr_result;
        end
        
        p0_recon = reshape(p0_recon, grid_dims);
        p0_recon = max(p0_recon, 0);
        tr_info.max_p0_recon = max(p0_recon(:));
        
    catch ME
        tr_info.success = false;
        tr_info.error = ME.message;
        recon_dose = zeros(grid_dims);
        return;
    end
    
    % Convert back to dose
    conv_factor = medium.gruneisen .* medium.density;
    conv_factor(conv_factor < 1e-10) = 1;
    
    recon_dose_per_pulse = p0_recon ./ conv_factor;
    recon_dose = recon_dose_per_pulse * num_pulses;
end
