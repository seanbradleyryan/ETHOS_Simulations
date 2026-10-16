function tables = define_tissue_tables()
%DEFINE_TISSUE_TABLES Tissue property lookup tables for HU thresholding.
%
%   tables = define_tissue_tables()
%
%   Returns the HU->acoustic-property lookup tables used by
%   create_acoustic_medium. Two threshold schemes are provided:
%     .threshold_1 - 9-tissue HU anchor table, linearly interpolated between
%                    anchors; water outside the body contour
%     .threshold_2 - 4-tissue simplified model (+ real in-body air row)
%
%   This is the shared, canonical copy (extracted verbatim from the local
%   define_tissue_tables previously duplicated inside run_standalone_simulation
%   and the other simulation drivers). get_default_config populates
%   CONFIG.tissue_tables from this file and adds a '.uniform' row built from the
%   CONFIG.uniform_* values.
%
%   See also: create_acoustic_medium, get_default_config

    % THRESHOLD_1: HU anchor points, linearly interpolated in between.
    % Each column is one anchor. A tissue with an HU range gets two anchors
    % (start and end of the range) with equal values, so it is flat across
    % its range. HU below the first / above the last anchor take the end
    % values. NaN = no sourced value: that anchor is skipped for that
    % property, so it is interpolated from its neighbours.
    %   Blood (+13..+50) and Muscle (+35..+55) overlap: Blood is flat 13-35,
    %   then the 35-50 interpolation blends Blood -> Muscle.
    %   Bone ends at 3000 and Metal is >3000: the 3000 -> 3001 anchors make
    %   that a step.
    tables.threshold_1 = struct();
    tables.threshold_1.hu_points    = [-1000, -700,  -600,  -120,  -90,   0,      13,    35,    50,     55,     100,         300,         700,    3000,   3001];
    tables.threshold_1.point_names  = {'Air', 'Lung','Lung','Fat', 'Fat', 'Water','Blood','Blood','Muscle','Muscle','SoftTissue','SoftTissue','Bone', 'Bone', 'Metal'};
    tables.threshold_1.density      = [1.2,   400,   400,   911,   911,   993,    1060,  1060,  1090,   1090,   1090,        1090,        1850,   1850,   4500];
    tables.threshold_1.sound_speed  = [343,   644,   644,   1440,  1440,  1524,   1584,  1584,  1588,   1588,   1610,        1610,        2800,   2800,   6070];
    tables.threshold_1.alpha_coeff  = [0,     NaN,   NaN,   0.6,   0.6,   0.002,  0.2,   0.2,   1.1,    1.1,    0.4,         0.4,         8,      8,      0];
    tables.threshold_1.gruneisen    = [0.38,  0.2,   0.2,   0.8,   0.8,   0.20,   0.23,  0.23,  0.21,   0.21,   0.21,        0.21,        0.79,   0.79,   1.2];
    tables.threshold_1.alpha_power  = 1.1;   % single exponent for every voxel
    % Voxels outside the body contour are given this HU (water) before the
    % lookup, so the sensor-coupling region is a water bath.
    tables.threshold_1.outside_body_hu = 0;

    tables.threshold_2 = struct();
    tables.threshold_2.hu_boundaries = [-1000, -200, -50, 100, Inf];
    tables.threshold_2.tissue_names  = {'Water','Fat','SoftTissue','Bone'};
    tables.threshold_2.density       = [1000,   920, 1040,         1900];
    tables.threshold_2.sound_speed   = [1480,  1450, 1540,         3200];
    tables.threshold_2.alpha_coeff   = [0.0022, 0.48, 0.5,         10];
    tables.threshold_2.alpha_power   = [2.0,    1.5,  1.1,         1.0];
    tables.threshold_2.gruneisen     = [0.11,   0.7,  1.0,         1.0];
    % Real air row (applied only to low-HU voxels INSIDE the body contour;
    % outside-body low-HU stays water for sensor coupling). gruneisen = 0 so
    % air generates no PA signal; these voxels are masked from the recon dose.
    tables.threshold_2.air_hu_threshold = -300;
    tables.threshold_2.air = struct( ...
        'density',     1.2, ...
        'sound_speed', 343, ...
        'alpha_coeff', 0, ...
        'alpha_power', 1.0, ...
        'gruneisen',   0);
end
