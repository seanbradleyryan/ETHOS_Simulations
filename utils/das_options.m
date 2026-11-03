function opts = das_options(config)
%DAS_OPTIONS Delay-and-sum options (das_reconstruct opts) read from a CONFIG.
%
%   opts = das_options(config)
%   defaults = das_options(struct())
%
%   PURPOSE:
%   Single definition of the delay-and-sum settings and their defaults. Each
%   CONFIG field das_<name> becomes opts.<name>; a missing field takes the
%   default below. Used by run_single_field_simulation (the pipeline DAS),
%   noise_ensemble_error_bars (the noise floor), compute_sim_config_hash (which
%   leaves a setting out of the hash when it equals its default) and
%   das_reconstruct (defaults for callers that pass a partial opts).
%
%   INPUTS:
%       config - CONFIG struct (full or partial). Read fields:
%           .das_use_elements    (true)      element vs per-voxel channels
%           .das_envelope        (true)      Hilbert envelope of the result
%           .das_use_aperture    (true)      apply the directivity cone
%           .das_aperture_cos    (0.997)     cone threshold (cos of half-angle)
%           .das_depth_weight    (true)      depth/r^2 vs plain 1/r^2 weighting
%           .das_interp          ('nearest') 'nearest' | 'linear'
%           .das_use_operator    (true)      precomputed DAS matrix (speed only)
%           .das_operator_max_gb (8)         memory limit for building the matrix
%           .das_check_operator  (false)     also run the channel loop and print
%                                            the difference (verification)
%           .use_gpu             (false)     multiply the matrix on the GPU
%       The first six change the reconstruction and are in the config hash
%       when they differ from these defaults; the last four only change speed
%       or logging and are not hashed.
%
%   OUTPUTS:
%       opts - Struct with fields use_elements, envelope, use_aperture,
%              aperture_cos, depth_weight, interp, use_operator,
%              max_operator_gb, check_operator, use_gpu.
%
%   EXAMPLE:
%       opts = das_options(CONFIG);
%       recon = das_reconstruct(sensorData, sensor, sensor_info, medium, ...
%                               Nx, Ny, Nz, dx, dy, dz, dt, opts);
%
%   See also: das_reconstruct, compute_sim_config_hash, run_single_field_simulation

    opts = struct();
    opts.use_elements    = read_field(config, 'das_use_elements',    true);
    opts.envelope        = read_field(config, 'das_envelope',        true);
    opts.use_aperture    = read_field(config, 'das_use_aperture',    true);
    opts.aperture_cos    = read_field(config, 'das_aperture_cos',    0.997);
    opts.depth_weight    = read_field(config, 'das_depth_weight',    true);
    opts.interp          = read_field(config, 'das_interp',          'nearest');
    opts.use_operator    = read_field(config, 'das_use_operator',    true);
    opts.max_operator_gb = read_field(config, 'das_operator_max_gb', 8);
    opts.check_operator  = read_field(config, 'das_check_operator',  false);
    opts.use_gpu         = read_field(config, 'use_gpu',             false);
end


function val = read_field(config, name, default_val)
    % config.(name) when present and not empty, else default_val.
    if isfield(config, name) && ~isempty(config.(name))
        val = config.(name);
    else
        val = default_val;
    end
end
