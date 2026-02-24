function [sensor_data_avg, sensor_data_expanded] = apply_element_averaging(sensor_data, sensor_info)
%APPLY_ELEMENT_AVERAGING Average sensor data within element patches and expand back
%
%   [sensor_data_avg, sensor_data_expanded] = apply_element_averaging(sensor_data, sensor_info)
%
%   PURPOSE:
%   After kspaceFirstOrder3D returns sensor_data (one time series per sensor
%   voxel), this function averages rows that share the same element_map index
%   to simulate finite-element signal averaging. For time reversal, the
%   averaged data is expanded back to full sensor resolution by duplicating
%   the averaged signal to all voxels within that element.
%
%   INPUTS:
%       sensor_data - [num_sensor_voxels x Nt] array from kspaceFirstOrder3D.
%                     Rows correspond to true voxels in sensor_mask, ordered
%                     by linear index (column-major).
%       sensor_info - Struct from determine_sensor_mask with:
%           .element_map     - 2D array (local_nx x local_nz) mapping each
%                              sensor voxel to an element index. 0 = removed voxel.
%           .num_elements    - Number of unique elements.
%           .sensor_x_range  - [x_start, x_end] voxel indices.
%           .sensor_z_range  - [z_start, z_end] voxel indices.
%           .sensor_y_index  - Y index of the sensor plane.
%
%   OUTPUTS:
%       sensor_data_avg     - [num_elements x Nt] averaged sensor data.
%       sensor_data_expanded - [num_sensor_voxels x Nt] expanded data where
%                              each voxel carries its element's averaged signal.
%                              Use this for time reversal input.
%
%   NOTES:
%       - If element_map is empty or num_elements == 0, returns inputs unchanged.
%       - The mapping from sensor_data rows to element indices relies on the
%         sensor mask being ordered by linear index (MATLAB column-major).
%         The sensor mask is at a single Y index, so the linear ordering
%         within that Y-plane is column-major over (X, Z).
%
%   EXAMPLE:
%       [sensor_data_avg, sensor_data_exp] = apply_element_averaging(sensor_data, sensor_info);
%       % Use sensor_data_exp for time reversal
%       source_tr.p = fliplr(sensor_data_exp);
%
%   AUTHOR: ETHOS Pipeline Team
%   DATE: February 2026
%   VERSION: 1.0
%
%   See also: determine_sensor_mask, run_single_field_simulation

    % Passthrough if no element averaging
    if isempty(sensor_info.element_map) || sensor_info.num_elements == 0
        sensor_data_avg = sensor_data;
        sensor_data_expanded = sensor_data;
        return;
    end
    
    num_voxels = size(sensor_data, 1);
    Nt = size(sensor_data, 2);
    num_elements = sensor_info.num_elements;
    element_map = sensor_info.element_map;
    
    % Flatten element_map to a vector matching the sensor voxel ordering.
    % The sensor mask is at a single Y index, and sensor_data rows are ordered
    % by linear index within the 3D grid. For a planar sensor at one Y index,
    % the ordering is column-major over (Y=fixed, X, Z) which means X varies
    % fastest, then Z.
    %
    % element_map is (local_nx x local_nz), where local_nx = X extent.
    % Flatten column-major: element_map(:) gives X-first ordering, which
    % matches the linear index ordering of the sensor plane.
    
    elem_vec = element_map(:);  % [local_nx * local_nz x 1]
    
    % Remove entries where element_map == 0 (removed voxels during validation)
    valid_mask = elem_vec > 0;
    elem_vec_valid = elem_vec(valid_mask);
    
    % Sanity check: number of valid elements should match sensor voxels
    if length(elem_vec_valid) ~= num_voxels
        warning('apply_element_averaging:SizeMismatch', ...
            'Element map valid entries (%d) != sensor voxels (%d). Skipping averaging.', ...
            length(elem_vec_valid), num_voxels);
        sensor_data_avg = sensor_data;
        sensor_data_expanded = sensor_data;
        return;
    end
    
    % Average sensor data within each element
    sensor_data_avg = zeros(num_elements, Nt);
    element_counts = zeros(num_elements, 1);
    
    for e = 1:num_elements
        member_rows = find(elem_vec_valid == e);
        if ~isempty(member_rows)
            sensor_data_avg(e, :) = mean(sensor_data(member_rows, :), 1);
            element_counts(e) = length(member_rows);
        end
    end
    
    fprintf('        [Sensor] Element averaging: %d elements, %d-%d voxels/element\n', ...
        num_elements, min(element_counts(element_counts > 0)), max(element_counts));
    
    % Expand back: duplicate averaged signal to all voxels in each element
    sensor_data_expanded = zeros(num_voxels, Nt);
    for v = 1:num_voxels
        e = elem_vec_valid(v);
        if e > 0
            sensor_data_expanded(v, :) = sensor_data_avg(e, :);
        end
    end
end
