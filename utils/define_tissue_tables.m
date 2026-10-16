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
%     .continuous_1 - Yu et al. 2019 Table I: piecewise-linear HU -> density
%                    -> sound speed / alpha_V / Cp / alpha_0; water outside
%                    the body contour
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

    % CONTINUOUS_1: Yu et al., Med Phys 46(8):3649 (2019), Table I.
    % HU -> density rho [g/cm^3], then rho -> sound speed, thermal expansion
    % alpha_V, specific heat Cp and attenuation prefactor alpha_0.
    % Gamma = alpha_V * c^2 / Cp is computed in create_acoustic_medium.
    % Each conversion is a list of straight-line pieces, one row per piece:
    %   [from, to, slope, intercept]  ->  y = slope * x + intercept  for from <= x < to
    tables.continuous_1 = struct();
    tables.continuous_1.hu_to_rho = [ ...        % x = HU, y = rho [g/cm^3]
        -Inf,    -98,    1.03e-3,     1.0309
        -98,     -1,     0.9383e-3,   1.0229
        -1,      11,    -0.475e-3,    1.0215
        11,      102,    1.07e-3,     1.0046
        102,     385,    0.4703e-3,   1.0657
        385,     Inf,    0.5898e-3,   1.0198];
    tables.continuous_1.rho_to_c = [ ...         % x = rho [g/cm^3], y = sound speed [m/s]
        -Inf,    0.93,   1173.1,      341.6344
        0.93,    1.12,   1120,        391
        1.12,    Inf,    2043.4,     -737.8876];
    tables.continuous_1.rho_to_alphaV = [ ...    % x = rho [g/cm^3], y = alpha_V [1/K]
        -Inf,    0.93,  -24.78e-4,    32.27e-4
        0.93,    0.9713, 0,           9.2e-4
        0.9713,  0.9882,-317.4556e-4, 317.5e-4
        0.9882,  1.12,   0,           3.835e-4
        1.12,    1.908, -3.8198e-4,   8.1132e-4
        1.908,   Inf,    0,           0.825e-4];
    tables.continuous_1.rho_to_Cp = [ ...        % x = rho [g/cm^3], y = Cp [J/(kg K)]
        -Inf,    0.93,   1466.1,      1002
        0.93,    0.9713, 0,           2348.3
        0.9713,  0.9882, 73193.2,    -68744.3
        0.9882,  1.12,   0,           3585.3
        1.12,    1.908, -2779.3,      6616.2
        1.908,   Inf,    0,           1312.8];
    tables.continuous_1.rho_to_alpha0 = [ ...    % x = rho [g/cm^3], y = alpha_0 [dB/(MHz^y cm)]
        -Inf,    0.93,   0.4265,     -0.0029
        0.93,    0.9713, 0,           0.4
        0.9713,  0.9882, 11.8343,    -11.0947
        0.9882,  1.12,   0,           0.6
        1.12,    Inf,    4.7,        -4.2248];
    % HU below Yu's lowest fitted value gives negative density, so it is
    % raised to this value before the lookup.
    tables.continuous_1.min_hu      = -1000;
    tables.continuous_1.alpha_power = 1.1;   % single exponent for every voxel
    % Water bath outside the body contour (same values as threshold_2 water).
    tables.continuous_1.outside_body = struct( ...
        'density',     1000, ...
        'sound_speed', 1480, ...
        'alpha_coeff', 0.0022, ...
        'gruneisen',   0.11);
end
