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
%           .voxel_element_idx - Nvox x 1 vector mapping each sensor voxel to
%                                its element index (0 = unassigned). Ordered by
%                                find(sensor_mask) linear index, which matches
%                                the sensor_data row order from kspaceFirstOrder3D.
%                                AUTHORITATIVE mapping; required for tilted sensors.
%           .element_map       - LEGACY 2D fallback (local_nx x local_nz). Used
%                                only when voxel_element_idx is missing.
%           .num_elements      - Number of unique elements.
%
%   OUTPUTS:
%       sensor_data_avg     - [num_elements x Nt] averaged sensor data.
%       sensor_data_expanded - [num_sensor_voxels x Nt] expanded data where
%                              each voxel carries its element's averaged signal.
%                              Use this for time reversal input.
%
%   NOTES:
%       - If both voxel_element_idx and element_map are empty (or num_elements
%         == 0), returns inputs unchanged.
%       - Prefers voxel_element_idx (works for tilted and flat sensors); falls
%         back to flattening element_map column-major (legacy single-Y sensors
%         only).
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
    has_vox_idx = isfield(sensor_info, 'voxel_element_idx') && ...
                  ~isempty(sensor_info.voxel_element_idx);
    has_element_map = isfield(sensor_info, 'element_map') && ...
                      ~isempty(sensor_info.element_map);
    if (~has_vox_idx && ~has_element_map) || sensor_info.num_elements == 0
        sensor_data_avg = sensor_data;
        sensor_data_expanded = sensor_data;
        return;
    end

    num_voxels = size(sensor_data, 1);
    Nt = size(sensor_data, 2);
    num_elements = sensor_info.num_elements;

    if has_vox_idx
        % Authoritative mapping: one entry per true sensor voxel, ordered by
        % find(sensor_mask) linear index — exactly the order kspaceFirstOrder3D
        % returns sensor_data rows in. Works for both flat and tilted sensors.
        elem_vec_valid = double(sensor_info.voxel_element_idx(:));
        elem_vec_valid = elem_vec_valid(elem_vec_valid > 0);
    else
        % LEGACY fallback for sensors built before the tilt refactor: flatten
        % element_map column-major over (X, Z). Only correct when the sensor
        % lies on a single Y slab.
        elem_vec = sensor_info.element_map(:);
        valid_mask = elem_vec > 0;
        elem_vec_valid = elem_vec(valid_mask);
    end
    
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
