function [hash_hex, canonical] = compute_sim_config_hash(config)
%COMPUTE_SIM_CONFIG_HASH  Short hex digest of the sim-affecting CONFIG fields.
%
%   [hash_hex, canonical] = compute_sim_config_hash(config)
%
%   Builds a canonical sub-struct of CONFIG containing only the fields that
%   change the contents of a reconstructed dose, JSON-encodes it with stable
%   key order, and returns the first 8 hex characters of the MD5 digest plus
%   the canonical struct itself.
%
%   The field list is an explicit allow-list (sim_config_fields below), so
%   adding an unrelated CONFIG field later does NOT invalidate existing
%   caches.  Missing fields are recorded as the sentinel '__unset__' so
%   absence participates in the hash deterministically.
%
%   This is the single source of truth for the config hash used by:
%     - pipeline_simulate.m  (writes outputs keyed on the hash)
%     - step3_analysis.m     (reads outputs keyed on the hash)
%
%   INPUTS:
%       config - CONFIG struct (full or partial)
%
%   OUTPUTS:
%       hash_hex  - 8-character hex string (e.g. 'a1b2c3d4')
%       canonical - Struct containing only the allow-listed fields, with
%                   '__unset__' filled in for any that were missing
%
%   See also: pipeline_simulate, step3_analysis

    allowed   = sim_config_fields();
    canonical = struct();
    for i = 1:numel(allowed)
        name = allowed{i};
        if isfield(config, name)
            canonical.(name) = config.(name);
        else
            canonical.(name) = '__unset__';
        end
    end

    % Backward compatibility: 'gaussian' was the only (implicit) pulse shape
    % before pulse_shape existed, so hash it as the unset sentinel. This keeps
    % every pre-existing recon cache valid; only non-gaussian shapes (e.g.
    % 'rectangular') change the hash and get their own output files.
    if isfield(canonical, 'pulse_shape') && (ischar(canonical.pulse_shape) || isstring(canonical.pulse_shape)) ...
            && strcmpi(canonical.pulse_shape, 'gaussian')
        canonical.pulse_shape = '__unset__';
    end

    % Backward compatibility: sensor_side only matters for the lateral placement
    % method, and 'right' is its default. Hash the default as the unset sentinel
    % so adding sensor_side does not invalidate any existing (anterior) caches;
    % only 'left' changes the hash and gets its own output files.
    if isfield(canonical, 'sensor_side') && (ischar(canonical.sensor_side) || isstring(canonical.sensor_side)) ...
            && strcmpi(canonical.sensor_side, 'right')
        canonical.sensor_side = '__unset__';
    end

    json     = jsonencode(canonical);
    md       = java.security.MessageDigest.getInstance('MD5');
    md.update(uint8(json));
    raw      = typecast(md.digest(), 'uint8');
    hex      = sprintf('%02x', raw);
    hash_hex = hex(1:8);
end


function fields = sim_config_fields()
    % Allow-listed CONFIG fields that affect numerical reconstruction output.
    % use_gpu, paths, analysis knobs, and parallel controls are
    % intentionally excluded.  Sorted for deterministic JSON key order.
    fields = sort({ ...
        'alpha_power', ...
        'blind_recon_ct3', ...
        'blind_recon_reference_label', ...
        'cfl_number', ...
        'conv_deconv_lambda', ...
        'conv_noise_level', ...
        'convergence_tol', ...
        'convolution_kernel', ...
        'correction_factor', ...
        'dose_per_pulse_cGy', ...
        'downscale_factor', ...
        'force_uniform_attenuation', ...
        'force_uniform_density', ...
        'force_uniform_gruneisen', ...
        'force_uniform_speed', ...
        'gruneisen_method', ...
        'normalize', ...
        'num_time_reversal_iter', ...
        'pml_size', ...
        'pulse_shape', ...
        'reconstruction_method', ...
        'sensor_placement_method', ...
        'sensor_side', ...
        'sensor_x_index', ...
        'sensor_y_index', ...
        'use_attenuation', ...
        'use_grid_padding', ...
        'use_pressure_scale_correction'});
end
