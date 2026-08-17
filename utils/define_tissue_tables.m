function tables = define_tissue_tables()
%DEFINE_TISSUE_TABLES Tissue property lookup tables for HU thresholding.
%
%   tables = define_tissue_tables()
%
%   Returns the HU->acoustic-property lookup tables used by
%   create_acoustic_medium. Two threshold schemes are provided:
%     .threshold_1 - 9-tissue model
%     .threshold_2 - 4-tissue simplified model (+ real in-body air row)
%
%   This is the shared, canonical copy (extracted verbatim from the local
%   define_tissue_tables previously duplicated inside run_standalone_simulation
%   and the other simulation drivers). get_default_config populates
%   CONFIG.tissue_tables from this file and adds a '.uniform' row built from the
%   CONFIG.uniform_* values.
%
%   See also: create_acoustic_medium, get_default_config

    tables.threshold_1 = struct();
    tables.threshold_1.hu_boundaries = [-1000, -900, -500, -200, -50, 13, 50, 80, 300, 3000, Inf];
    tables.threshold_1.tissue_names  = {'Air','Lung','Fat','Water','Blood','Muscle','SoftTissue','Bone','Metal'};
    tables.threshold_1.density       = [1.2,  400,  920, 1000, 1060, 1050, 1040, 1900, 7800];
    tables.threshold_1.sound_speed   = [343,  600, 1450, 1480, 1575, 1580, 1540, 3200, 5900];
    tables.threshold_1.alpha_coeff   = [0,   0.5, 0.48, 0.0022, 0.2, 0.5, 0.5,  10,   0];
    tables.threshold_1.alpha_power   = [1.0, 1.5, 1.5,  2.0,   1.3, 1.0, 1.1,  1.0,  1.0];
    tables.threshold_1.gruneisen     = [0,   0.5, 0.7,  0.11,  0.15, 0.2, 1.0,  0,    0];

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
