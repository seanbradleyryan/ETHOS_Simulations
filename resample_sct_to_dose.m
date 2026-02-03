%% resample_sct_to_dose.m
% Script to load Raystation dose file, load corresponding SCT, resample and
% crop SCT to match dose grid dimensions, and export as .mat file
%
% Author: Generated for ETHOS IMRT Processing Pipeline
% Date: 2026-02-03
%
% Purpose: Prepare SCT data that is spatially aligned and dimensionally
%          matched to the dose grid for subsequent k-Wave simulations

%% User Configuration
% Enable/disable plotting for verification
ENABLE_PLOTTING = true;  % Set to false to disable visualization

% Patient and session identifiers (loopable lists)
patient_ids = {'1194203'}; % Add more patient IDs as needed
sessions = {'Session_1'};  % Add more sessions as needed

%% Setup Paths and Loop Through Patients/Sessions
wd = '/mnt/weka/home/80030361/ETHOS_Simulations';

% Loop through each patient
for p_idx = 1:length(patient_ids)
    patient_id = patient_ids{p_idx};
    
    % Loop through each session
    for s_idx = 1:length(sessions)
        session = sessions{s_idx};
        
        fprintf('\n========================================\n');
        fprintf('Processing Patient: %s, Session: %s\n', patient_id, session);
        fprintf('========================================\n');
        
        raystation_dir = fullfile(wd, 'RaystationFiles', patient_id, session);
        sct_dir = fullfile(wd, 'EthosExports', patient_id, 'Pancreas', session, 'sct');

        % Verify directories exist
        if ~exist(raystation_dir, 'dir')
            warning('Raystation directory does not exist: %s. Skipping...', raystation_dir);
            continue;
        end
        if ~exist(sct_dir, 'dir')
            warning('SCT directory does not exist: %s. Skipping...', sct_dir);
            continue;
        end
        
        % Find RD*.dcm files in raystation directory
        rd_files = dir(fullfile(raystation_dir, 'RD*.dcm'));
        if isempty(rd_files)
            warning('No RD*.dcm files found in: %s. Skipping...', raystation_dir);
            continue;
        end
        
        % Use the first RD file found
        dose_filename = rd_files(1).name;
        fprintf('Found %d RD file(s). Using: %s\n', length(rd_files), dose_filename);

%% Load Raystation Dose
        fprintf('Loading Raystation dose file...\n');
        dose_filepath = fullfile(raystation_dir, dose_filename);

        % Read DICOM dose file
        dose_info = dicominfo(dose_filepath);
        dose_grid = double(squeeze(dicomread(dose_filepath)));

        % Apply dose grid scaling
        if isfield(dose_info, 'DoseGridScaling')
            dose_scaling = dose_info.DoseGridScaling;
            dose_grid = dose_grid * dose_scaling;
            fprintf('Applied dose grid scaling: %e\n', dose_scaling);
        end

        % Extract dose grid geometry
        dose_origin = dose_info.ImagePositionPatient; % [x, y, z] in mm
        dose_spacing = [dose_info.PixelSpacing(1); ...
                        dose_info.PixelSpacing(2); ...
                        dose_info.GridFrameOffsetVector(2) - dose_info.GridFrameOffsetVector(1)]; % [dx, dy, dz] in mm
        dose_dims = size(dose_grid); % [nx, ny, nz]

        fprintf('Dose grid dimensions: [%d, %d, %d]\n', dose_dims(1), dose_dims(2), dose_dims(3));
        fprintf('Dose grid spacing (mm): [%.3f, %.3f, %.3f]\n', dose_spacing(1), dose_spacing(2), dose_spacing(3));
        fprintf('Dose grid origin (mm): [%.3f, %.3f, %.3f]\n', dose_origin(1), dose_origin(2), dose_origin(3));

%% Load SCT Images
        fprintf('Loading SCT images...\n');

        % Get all DICOM files in SCT directory that contain 'CT' in the name
        all_files = dir(fullfile(sct_dir, '*.dcm'));
        sct_files = [];
        for i = 1:length(all_files)
            if contains(all_files(i).name, 'CT', 'IgnoreCase', true)
                sct_files = [sct_files; all_files(i)];
            end
        end
        
        if isempty(sct_files)
            warning('No CT DICOM files found in SCT directory: %s. Skipping...', sct_dir);
            continue;
        end
        
        fprintf('Found %d CT image files.\n', length(sct_files));

        % Sort files by Instance Number or Image Position
        sct_data = [];
        sct_positions = zeros(length(sct_files), 1);

        for i = 1:length(sct_files)
            info = dicominfo(fullfile(sct_dir, sct_files(i).name));
            sct_positions(i) = info.ImagePositionPatient(3); % z-position
        end

        [~, sort_idx] = sort(sct_positions);

        % Load sorted SCT slices
        fprintf('Reading %d SCT slices...\n', length(sct_files));
        for i = 1:length(sct_files)
            idx = sort_idx(i);
            img = dicomread(fullfile(sct_dir, sct_files(idx).name));
            if i == 1
                sct_data = zeros([size(img), length(sct_files)], 'int16');
            end
            sct_data(:, :, i) = img;
        end

        % Extract SCT geometry
        sct_info = dicominfo(fullfile(sct_dir, sct_files(sort_idx(1)).name));
        sct_origin = sct_info.ImagePositionPatient; % [x, y, z] in mm
        sct_spacing = [sct_info.PixelSpacing(1); ...
                       sct_info.PixelSpacing(2); ...
                       sct_info.SliceThickness]; % [dx, dy, dz] in mm
        sct_dims = size(sct_data); % [nx, ny, nz]

        fprintf('SCT dimensions: [%d, %d, %d]\n', sct_dims(1), sct_dims(2), sct_dims(3));
        fprintf('SCT spacing (mm): [%.3f, %.3f, %.3f]\n', sct_spacing(1), sct_spacing(2), sct_spacing(3));
        fprintf('SCT origin (mm): [%.3f, %.3f, %.3f]\n', sct_origin(1), sct_origin(2), sct_origin(3));

        % Convert to Hounsfield Units
        sct_hu = double(sct_data) * sct_info.RescaleSlope + sct_info.RescaleIntercept;

%% Resample and Crop SCT to Match Dose Grid
        fprintf('Resampling SCT to match dose grid...\n');

        % Create coordinate grids for dose
        [dose_x, dose_y, dose_z] = meshgrid(...
            dose_origin(1) + (0:dose_dims(2)-1) * dose_spacing(1), ...
            dose_origin(2) + (0:dose_dims(1)-1) * dose_spacing(2), ...
            dose_origin(3) + (0:dose_dims(3)-1) * dose_spacing(3));

        % Create coordinate grids for SCT
        [sct_x, sct_y, sct_z] = meshgrid(...
            sct_origin(1) + (0:sct_dims(2)-1) * sct_spacing(1), ...
            sct_origin(2) + (0:sct_dims(1)-1) * sct_spacing(2), ...
            sct_origin(3) + (0:sct_dims(3)-1) * sct_spacing(3));

        % Interpolate SCT to dose grid coordinates
        fprintf('Performing 3D interpolation (this may take a moment)...\n');
        sct_resampled = interp3(sct_x, sct_y, sct_z, sct_hu, ...
                                dose_x, dose_y, dose_z, 'linear', -1000);

        fprintf('Resampled SCT dimensions: [%d, %d, %d]\n', size(sct_resampled, 1), ...
                size(sct_resampled, 2), size(sct_resampled, 3));

%% Create Output Structure
        output_data = struct();
        output_data.sct_resampled = sct_resampled; % Hounsfield Units
        output_data.spacing = dose_spacing;        % [dx, dy, dz] in mm
        output_data.origin = dose_origin;          % [x, y, z] in mm
        output_data.dimensions = dose_dims;        % [nx, ny, nz]
        output_data.patient_id = patient_id;
        output_data.session = session;
        output_data.dose_filename = dose_filename;
        output_data.timestamp = datetime('now');

%% Save Output
        output_filename = sprintf('sct_dose_matched_%s_%s.mat', patient_id, session);
        output_filepath = fullfile(raystation_dir, output_filename);

        fprintf('Saving matched SCT data to: %s\n', output_filepath);
        save(output_filepath, 'output_data', '-v7.3');

        fprintf('Save complete!\n');

%% Optional Visualization
        if ENABLE_PLOTTING
            fprintf('Generating verification plots...\n');
            
            % Select middle slice in each dimension
            slice_x = round(dose_dims(1) / 2);
            slice_y = round(dose_dims(2) / 2);
            slice_z = round(dose_dims(3) / 2);
            
            figure('Name', 'SCT and Dose Grid Verification', 'Position', [100, 100, 1400, 800]);
            
            % Axial view (z-slice)
            subplot(2, 3, 1);
            imagesc(squeeze(sct_resampled(:, :, slice_z)));
            colorbar;
            colormap(gca, 'gray');
            title(sprintf('SCT Axial (z=%d)', slice_z));
            xlabel('X'); ylabel('Y');
            axis equal tight;
            
            subplot(2, 3, 4);
            imagesc(squeeze(dose_grid(:, :, slice_z)));
            colorbar;
            title(sprintf('Dose Axial (z=%d) - Reference', slice_z));
            xlabel('X'); ylabel('Y');
            axis equal tight;
            
            % Sagittal view (x-slice)
            subplot(2, 3, 2);
            imagesc(squeeze(sct_resampled(slice_x, :, :))');
            colorbar;
            colormap(gca, 'gray');
            title(sprintf('SCT Sagittal (x=%d)', slice_x));
            xlabel('Y'); ylabel('Z');
            axis equal tight;
            
            subplot(2, 3, 5);
            imagesc(squeeze(dose_grid(slice_x, :, :))');
            colorbar;
            title(sprintf('Dose Sagittal (x=%d) - Reference', slice_x));
            xlabel('Y'); ylabel('Z');
            axis equal tight;
            
            % Coronal view (y-slice)
            subplot(2, 3, 3);
            imagesc(squeeze(sct_resampled(:, slice_y, :))');
            colorbar;
            colormap(gca, 'gray');
            title(sprintf('SCT Coronal (y=%d)', slice_y));
            xlabel('X'); ylabel('Z');
            axis equal tight;
            
            subplot(2, 3, 6);
            imagesc(squeeze(dose_grid(:, slice_y, :))');
            colorbar;
            title(sprintf('Dose Coronal (y=%d) - Reference', slice_y));
            xlabel('X'); ylabel('Z');
            axis equal tight;
            
            % Add overall title
            sgtitle(sprintf('Patient: %s, Session: %s', patient_id, session), ...
                    'FontSize', 14, 'FontWeight', 'bold');
            
            fprintf('Verification plots displayed.\n');
        end

        fprintf('\n=== Processing Complete for %s / %s ===\n', patient_id, session);
        fprintf('Output file: %s\n', output_filename);
        fprintf('Output contains:\n');
        fprintf('  - sct_resampled: Resampled SCT in HU [%d x %d x %d]\n', dose_dims(1), dose_dims(2), dose_dims(3));
        fprintf('  - spacing: Voxel spacing [%.3f, %.3f, %.3f] mm\n', dose_spacing(1), dose_spacing(2), dose_spacing(3));
        fprintf('  - origin: Grid origin [%.3f, %.3f, %.3f] mm\n', dose_origin(1), dose_origin(2), dose_origin(3));
        fprintf('  - Additional metadata\n');
        
    end % End session loop
end % End patient loop

fprintf('\n========================================\n');
fprintf('All processing complete!\n');
fprintf('========================================\n');