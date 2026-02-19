%% VISUALIZE_SENSOR
%
%   Standalone script to visualize the ultrasound sensor placement on CT
%   images in the three anatomical planes (transverse, coronal, sagittal)
%   running through the center of the sensor. The sensor is overlaid in
%   bright red on grayscale CT.
%
%   PREREQUISITES:
%       - step15_process_doses must have been run successfully for the
%         target patient/session (creates processed/ directory with
%         sct_resampled.mat, field_dose_*.mat, and metadata.mat)
%       - determine_sensor_mask.m must be on the MATLAB path
%       - Image Processing Toolbox
%
%   USAGE:
%       1. Set patient_id, session, and config below
%       2. Optionally adjust field_index to select which beam field to use
%       3. Run the script
%
%   OUTPUTS:
%       Figure 1: Three-panel CT with sensor overlay (red) in transverse,
%                 coronal, and sagittal planes through the sensor center
%
%   COORDINATE SYSTEM:
%       Dose grid uses DICOM patient coordinates:
%           X = patient left-right (columns, dim 2)
%           Y = patient anterior-posterior (rows, dim 1). Lower index = more anterior.
%           Z = patient superior-inferior (slices, dim 3)
%       Array indexing: array(Y, X, Z)
%
%   AUTHOR: ETHOS Pipeline Team
%   DATE: February 2026
%   VERSION: 1.0
%
%   See also: determine_sensor_mask, step15_process_doses, load_processed_data

%% ======================== CONFIGURATION ========================

% --- Patient / Session ---
patient_id = '1194203';
session    = 'Session_1';

% --- Paths ---
config.working_dir    = '/mnt/weka/home/80030361/ETHOS_Simulations';
config.treatment_site = 'Pancreas';

% --- Field Selection ---
% Which field dose to use for sensor placement (1-based index into
% field_dose_*.mat files). Set to 1 to use the first field.
field_index = 1;

% --- Sensor Parameters (passed to determine_sensor_mask) ---
config.sensor_size_cm     = [10, 10];   % Physical sensor dims [X, Z] in cm
config.sensor_standoff_mm = 5;          % Gap between body surface and sensor (mm)
config.element_size_mm    = [];         % Element patch size (empty = no grouping)
config.jaw_margin_mm      = 10;         % Extra margin around jaw projection (mm)
config.sensor_placement   = 'anterior'; % Placement side
config.pml_size           = 10;         % PML thickness in voxels

% --- Display Parameters ---
hu_window = [-1000, 1000];  % HU display window [min, max]
sensor_color = [1, 0, 0];  % Bright red for sensor overlay
sensor_alpha = 0.7;        % Sensor overlay opacity (0 = transparent, 1 = opaque)

%% ======================== CONSTRUCT PATHS ========================

fprintf('\n=========================================================\n');
fprintf('  Sensor Placement Visualization\n');
fprintf('  Patient: %s, Session: %s\n', patient_id, session);
fprintf('=========================================================\n\n');

rs_dir = fullfile(config.working_dir, 'RayStationFiles', patient_id, session);
processed_dir = fullfile(rs_dir, 'processed');
sct_dir = fullfile(config.working_dir, 'EthosExports', patient_id, ...
    config.treatment_site, session, 'sct');

%% ======================== VERIFY FILES EXIST ========================

if ~isfolder(processed_dir)
    error('visualize_sensor:NotProcessed', ...
        ['Processed directory not found:\n  %s\n' ...
         'Run step15_process_doses first.'], processed_dir);
end

sct_file  = fullfile(processed_dir, 'sct_resampled.mat');
meta_file = fullfile(processed_dir, 'metadata.mat');

required_files = {sct_file, meta_file};
required_names = {'sct_resampled.mat', 'metadata.mat'};
for k = 1:length(required_files)
    if ~isfile(required_files{k})
        error('visualize_sensor:FileNotFound', ...
            '%s not found in:\n  %s\nRun step15_process_doses first.', ...
            required_names{k}, processed_dir);
    end
end

%% ======================== LOAD SCT RESAMPLED ========================

fprintf('[1/4] Loading SCT resampled data...\n');
loaded = load(sct_file);
sct_resampled = loaded.sct_resampled;
fprintf('  Grid: [%d x %d x %d]\n', size(sct_resampled.cubeHU));
fprintf('  Spacing (mm): [%.2f, %.2f, %.2f]\n', sct_resampled.spacing);
fprintf('  HU range: [%.0f, %.0f]\n', ...
    min(sct_resampled.cubeHU(:)), max(sct_resampled.cubeHU(:)));

% Ensure bodyMask exists (load from tissue_masks.mat if needed)
if ~isfield(sct_resampled, 'bodyMask')
    tissue_file = fullfile(processed_dir, 'tissue_masks.mat');
    if isfile(tissue_file)
        fprintf('  Loading tissue masks for bodyMask...\n');
        tissue_data = load(tissue_file);
        if isfield(tissue_data, 'body_mask')
            sct_resampled.bodyMask = tissue_data.body_mask;
        end
        if isfield(tissue_data, 'couch_mask')
            sct_resampled.couchMask = tissue_data.couch_mask;
        end
    end
    if ~isfield(sct_resampled, 'bodyMask')
        error('visualize_sensor:NoBodyMask', ...
            'bodyMask not found in sct_resampled or tissue_masks.mat.');
    end
end

% Ensure couchMask exists (default to empty if not)
if ~isfield(sct_resampled, 'couchMask')
    fprintf('  [WARNING] No couchMask found, using empty mask.\n');
    sct_resampled.couchMask = false(size(sct_resampled.bodyMask));
end

%% ======================== LOAD FIELD DOSE ========================

fprintf('[2/4] Loading field dose %d...\n', field_index);

field_file = fullfile(processed_dir, sprintf('field_dose_%03d.mat', field_index));
if ~isfile(field_file)
    error('visualize_sensor:FieldNotFound', ...
        'Field dose file not found: %s', field_file);
end

loaded = load(field_file);
field_dose = loaded.field_dose;
fprintf('  Gantry angle: %.1f deg\n', field_dose.gantry_angle);
fprintf('  Max dose: %.4f Gy\n', field_dose.max_dose_Gy);

%% ======================== BUILD BEAM METADATA ========================

fprintf('[3/4] Building beam metadata from RTPLAN...\n');

beam_metadata = build_beam_metadata_for_sensor(sct_dir, processed_dir);

if isempty(beam_metadata)
    fprintf('  [WARNING] No beam metadata available. Sensor placement will\n');
    fprintf('            proceed without beam exclusion zones.\n');
else
    fprintf('  Loaded metadata for %d beams\n', length(beam_metadata));
    for b = 1:length(beam_metadata)
        fprintf('    Beam %d: gantry %.1f deg\n', ...
            beam_metadata(b).beam_number, beam_metadata(b).gantry_angle);
    end
end

%% ======================== COMPUTE SENSOR MASK ========================

fprintf('[4/4] Computing sensor placement...\n');

[sensor_mask, sensor_info] = determine_sensor_mask( ...
    sct_resampled, field_dose, beam_metadata, config);

if ~sensor_info.placement_valid && sensor_info.num_sensor_voxels == 0
    error('visualize_sensor:NoPlacement', ...
        'Sensor placement failed. Check determine_sensor_mask warnings.');
end

fprintf('\n--- Sensor Placement Summary ---\n');
fprintf('  Y index (AP plane):    %d\n', sensor_info.sensor_y_index);
fprintf('  X range (LR voxels):   [%d, %d]\n', sensor_info.sensor_x_range);
fprintf('  Z range (SI voxels):   [%d, %d]\n', sensor_info.sensor_z_range);
fprintf('  Center (mm):           [%.1f, %.1f, %.1f]\n', sensor_info.sensor_center_mm);
fprintf('  Size (voxels):         [%d, %d]\n', sensor_info.sensor_size_voxels);
fprintf('  Sensor voxels:         %d\n', sensor_info.num_sensor_voxels);
fprintf('  Surface-to-dose (mm):  %.1f\n', sensor_info.surface_to_dose_distance_mm);
fprintf('  Placement valid:       %s\n', mat2str(sensor_info.placement_valid));

%% ======================== DETERMINE SLICE INDICES ========================
% Slice through the CENTER of the sensor in each anatomical plane

% Sensor center indices
center_y = sensor_info.sensor_y_index;                       % Single Y index (AP)
center_x = round(mean(sensor_info.sensor_x_range));          % Mid X (LR)
center_z = round(mean(sensor_info.sensor_z_range));          % Mid Z (SI)

fprintf('\n--- Visualization Slices ---\n');
fprintf('  Transverse (axial):  Z = %d  (slice through sensor center)\n', center_z);
fprintf('  Coronal:             Y = %d  (sensor AP plane)\n', center_y);
fprintf('  Sagittal:            X = %d  (slice through sensor center)\n', center_x);

%% ======================== EXTRACT CT AND SENSOR SLICES ========================

ct = sct_resampled.cubeHU;

% --- Transverse (axial): fixed Z slice ---
% Shows X (columns) vs Y (rows) at sensor center Z
trans_ct   = squeeze(ct(:, :, center_z));          % (Y, X)
trans_sens = squeeze(sensor_mask(:, :, center_z)); % (Y, X)

% --- Coronal: fixed Y slice ---
% Shows X (dim 2) vs Z (dim 3) at sensor Y plane
cor_ct   = squeeze(ct(center_y, :, :));            % (X, Z)
cor_sens = squeeze(sensor_mask(center_y, :, :));   % (X, Z)

% --- Sagittal: fixed X slice ---
% Shows Y (dim 1) vs Z (dim 3) at sensor center X
sag_ct   = squeeze(ct(:, center_x, :));            % (Y, Z)
sag_sens = squeeze(sensor_mask(:, center_x, :));   % (Y, Z)

%% ======================== BUILD RGB IMAGES ========================

trans_rgb = build_ct_sensor_overlay(trans_ct, trans_sens, hu_window, sensor_color, sensor_alpha);
cor_rgb   = build_ct_sensor_overlay(cor_ct',  cor_sens',  hu_window, sensor_color, sensor_alpha);
sag_rgb   = build_ct_sensor_overlay(sag_ct',  sag_sens',  hu_window, sensor_color, sensor_alpha);

%% ======================== COMPUTE PHYSICAL AXIS EXTENTS ========================

origin  = sct_resampled.origin;  % [x, y, z] mm
spacing = sct_resampled.spacing; % [dx, dy, dz] mm
dims    = size(ct);              % [Ny, Nx, Nz]

x_mm = origin(1) + (0:dims(2)-1) * spacing(1);  % columns
y_mm = origin(2) + (0:dims(1)-1) * spacing(2);  % rows
z_mm = origin(3) + (0:dims(3)-1) * spacing(3);  % slices

%% ======================== FIGURE: THREE-PANEL CT + SENSOR ========================

fig = figure('Name', 'Sensor Placement on CT', ...
    'Position', [40, 120, 1600, 550], 'Color', 'w');

% ----- Transverse (Axial) -----
ax1 = subplot(1, 3, 1);
imshow(trans_rgb, 'XData', x_mm, 'YData', y_mm);
axis on; axis equal tight;
hold on;
% Draw sensor bounding box
sx1 = origin(1) + (sensor_info.sensor_x_range(1) - 1) * spacing(1);
sx2 = origin(1) + (sensor_info.sensor_x_range(2) - 1) * spacing(1);
sy  = origin(2) + (sensor_info.sensor_y_index - 1) * spacing(2);
plot([sx1, sx2], [sy, sy], '-', 'Color', sensor_color, 'LineWidth', 2);
% Crosshair at sensor center
plot(origin(1) + (center_x - 1) * spacing(1), ...
     origin(2) + (center_y - 1) * spacing(2), ...
     '+', 'Color', [1 0.4 0.4], 'MarkerSize', 14, 'LineWidth', 1.5);
hold off;
xlabel('X — Left-Right (mm)');
ylabel('Y — Ant-Post (mm)');
title(sprintf('Transverse  (Z = %d, %.1f mm)', center_z, ...
    origin(3) + (center_z - 1) * spacing(3)), 'FontSize', 12);
set(ax1, 'YDir', 'normal');

% ----- Coronal -----
ax2 = subplot(1, 3, 2);
imshow(cor_rgb, 'XData', x_mm, 'YData', z_mm);
axis on; axis equal tight;
hold on;
% Draw sensor rectangle
sz1 = origin(3) + (sensor_info.sensor_z_range(1) - 1) * spacing(3);
sz2 = origin(3) + (sensor_info.sensor_z_range(2) - 1) * spacing(3);
rectangle('Position', [sx1, sz1, sx2 - sx1, sz2 - sz1], ...
    'EdgeColor', sensor_color, 'LineWidth', 2, 'LineStyle', '-');
% Crosshair
plot(origin(1) + (center_x - 1) * spacing(1), ...
     origin(3) + (center_z - 1) * spacing(3), ...
     '+', 'Color', [1 0.4 0.4], 'MarkerSize', 14, 'LineWidth', 1.5);
hold off;
xlabel('X — Left-Right (mm)');
ylabel('Z — Sup-Inf (mm)');
title(sprintf('Coronal  (Y = %d, %.1f mm)', center_y, ...
    origin(2) + (center_y - 1) * spacing(2)), 'FontSize', 12);
set(ax2, 'YDir', 'normal');

% ----- Sagittal -----
ax3 = subplot(1, 3, 3);
imshow(sag_rgb, 'XData', y_mm, 'YData', z_mm);
axis on; axis equal tight;
hold on;
% Draw sensor line (sensor is a single Y index spanning a Z range)
plot([sy, sy], [sz1, sz2], '-', 'Color', sensor_color, 'LineWidth', 2);
% Crosshair
plot(origin(2) + (center_y - 1) * spacing(2), ...
     origin(3) + (center_z - 1) * spacing(3), ...
     '+', 'Color', [1 0.4 0.4], 'MarkerSize', 14, 'LineWidth', 1.5);
hold off;
xlabel('Y — Ant-Post (mm)');
ylabel('Z — Sup-Inf (mm)');
title(sprintf('Sagittal  (X = %d, %.1f mm)', center_x, ...
    origin(1) + (center_x - 1) * spacing(1)), 'FontSize', 12);
set(ax3, 'YDir', 'normal');

% Super title
sgtitle(sprintf(['Sensor Placement — Patient %s, %s, Field %d ' ...
    '(Gantry %.1f°)\nSensor center: [%.1f, %.1f, %.1f] mm  |  ' ...
    'Size: %.0fx%.0f cm  |  Standoff: %.0f mm'], ...
    patient_id, session, field_index, field_dose.gantry_angle, ...
    sensor_info.sensor_center_mm(1), sensor_info.sensor_center_mm(2), ...
    sensor_info.sensor_center_mm(3), ...
    config.sensor_size_cm(1), config.sensor_size_cm(2), ...
    config.sensor_standoff_mm), ...
    'FontSize', 13, 'FontWeight', 'bold');

fprintf('\nVisualization complete.\n');

%% ========================================================================
%  LOCAL HELPER FUNCTIONS
%% ========================================================================

function rgb = build_ct_sensor_overlay(ct_slice, sensor_slice, hu_win, s_color, s_alpha)
%BUILD_CT_SENSOR_OVERLAY Create RGB image of CT with sensor highlighted in red
%
%   ct_slice     - 2D HU array
%   sensor_slice - 2D logical array (true = sensor voxel)
%   hu_win       - [min, max] HU display window
%   s_color      - [r, g, b] sensor color (0-1)
%   s_alpha      - Sensor overlay opacity (0-1)
%
%   Returns:
%       rgb - [M x N x 3] double image in [0, 1]

    % Normalize CT to grayscale [0, 1]
    ct_norm = (double(ct_slice) - hu_win(1)) / (hu_win(2) - hu_win(1));
    ct_norm = max(0, min(1, ct_norm));

    % Create grayscale RGB
    rgb = repmat(ct_norm, [1, 1, 3]);

    % Overlay sensor voxels
    if any(sensor_slice(:))
        mask = logical(sensor_slice);
        for c = 1:3
            ch = rgb(:, :, c);
            ch(mask) = (1 - s_alpha) * ch(mask) + s_alpha * s_color(c);
            rgb(:, :, c) = ch;
        end
    end
end


function beam_metadata = build_beam_metadata_for_sensor(sct_dir, processed_dir)
%BUILD_BEAM_METADATA_FOR_SENSOR Construct beam_metadata struct array
%
%   Attempts to load RTPLAN from sct_dir and extract beam_number,
%   gantry_angle, isocenter, jaw_x, and jaw_y for each beam.
%   Falls back to building minimal metadata from loaded field doses
%   if RTPLAN is unavailable.
%
%   INPUTS:
%       sct_dir       - Path to SCT directory containing RTPLAN*.dcm
%       processed_dir - Path to processed/ directory with field_dose_*.mat
%
%   OUTPUT:
%       beam_metadata - Struct array with fields needed by
%                       determine_sensor_mask, or empty [] on failure

    beam_metadata = [];

    % --- Try loading from RTPLAN ---
    rp_files = dir(fullfile(sct_dir, 'RTPLAN*.dcm'));
    if isempty(rp_files)
        rp_files = dir(fullfile(sct_dir, 'RP*.dcm'));
    end

    if ~isempty(rp_files)
        % Prefer adjusted MLC plan if available
        adjusted_idx = find(contains({rp_files.name}, 'adjusted_mlc'), 1);
        if ~isempty(adjusted_idx)
            rp_file = fullfile(sct_dir, rp_files(adjusted_idx).name);
        else
            rp_file = fullfile(sct_dir, rp_files(1).name);
        end

        try
            rtplan = dicominfo(rp_file);
            fprintf('  Loaded RTPLAN: %s\n', rp_files(1).name);

            if isfield(rtplan, 'BeamSequence')
                beam_fields = fieldnames(rtplan.BeamSequence);
                num_beams = length(beam_fields);

                for i = 1:num_beams
                    beam = rtplan.BeamSequence.(beam_fields{i});

                    % Beam number
                    if isfield(beam, 'BeamNumber')
                        beam_metadata(i).beam_number = beam.BeamNumber;
                    else
                        beam_metadata(i).beam_number = i;
                    end

                    % Gantry angle and isocenter from first control point
                    beam_metadata(i).gantry_angle = 0;
                    beam_metadata(i).isocenter = [];
                    beam_metadata(i).jaw_x = [];
                    beam_metadata(i).jaw_y = [];

                    if isfield(beam, 'ControlPointSequence')
                        cp_fields = fieldnames(beam.ControlPointSequence);
                        if ~isempty(cp_fields)
                            cp1 = beam.ControlPointSequence.(cp_fields{1});

                            if isfield(cp1, 'GantryAngle')
                                beam_metadata(i).gantry_angle = cp1.GantryAngle;
                            end

                            % Isocenter
                            if isfield(cp1, 'IsocenterPosition')
                                beam_metadata(i).isocenter = cp1.IsocenterPosition(:)';
                            end

                            % Jaw positions from BeamLimitingDevicePositionSequence
                            if isfield(cp1, 'BeamLimitingDevicePositionSequence')
                                bld_fields = fieldnames(cp1.BeamLimitingDevicePositionSequence);
                                for j = 1:length(bld_fields)
                                    bld = cp1.BeamLimitingDevicePositionSequence.(bld_fields{j});
                                    if isfield(bld, 'RTBeamLimitingDeviceType') && ...
                                       isfield(bld, 'LeafJawPositions')
                                        dev_type = upper(bld.RTBeamLimitingDeviceType);
                                        positions = bld.LeafJawPositions;
                                        if contains(dev_type, 'X') && ~contains(dev_type, 'MLC')
                                            beam_metadata(i).jaw_x = positions(:)';
                                        elseif contains(dev_type, 'Y') && ~contains(dev_type, 'MLC')
                                            beam_metadata(i).jaw_y = positions(:)';
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        catch ME
            warning('visualize_sensor:RTPLANError', ...
                'Failed to parse RTPLAN: %s\nProceeding without beam metadata.', ME.message);
            beam_metadata = [];
        end
    end

    % --- Fallback: build minimal metadata from field dose files ---
    if isempty(beam_metadata)
        fprintf('  Building beam metadata from field dose files (no jaw/isocenter)...\n');
        field_files = dir(fullfile(processed_dir, 'field_dose_*.mat'));
        if ~isempty(field_files)
            for i = 1:length(field_files)
                try
                    loaded = load(fullfile(processed_dir, field_files(i).name));
                    fd = loaded.field_dose;
                    beam_metadata(i).beam_number  = fd.field_num;
                    beam_metadata(i).gantry_angle = fd.gantry_angle;
                    beam_metadata(i).isocenter    = [];
                    beam_metadata(i).jaw_x        = [];
                    beam_metadata(i).jaw_y        = [];
                catch
                    continue;
                end
            end
        end
    end
end
