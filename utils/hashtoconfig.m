function config_block = hashtoconfig(hash, config)
%HASHTOCONFIG  Recover the canonical CONFIG block that produced a config hash.
%
%   config_block = hashtoconfig(hash, config)
%
%   The reconstruction hash (compute_sim_config_hash) is an MD5 digest and
%   cannot be inverted.  Instead, pipeline_simulate records every hash it runs
%   in a per-session config_registry.json, storing the canonical (allow-listed)
%   CONFIG sub-struct alongside each hash.  This util scans those registries
%   under <working_dir>/SimulationResults and returns the config block for the
%   requested hash.
%
%   INPUTS:
%       hash   - 8-character hex config hash (e.g. 'a1b2c3d4'). Case-insensitive.
%       config - CONFIG struct; only .working_dir is required.
%
%   OUTPUTS:
%       config_block - Struct of the allow-listed CONFIG fields that produced
%                      the hash (same content as compute_sim_config_hash's
%                      canonical output). Errors if the hash is not found.
%
%   EXAMPLE:
%       config.working_dir = '/mnt/weka/home/80030361/ETHOS_Simulations';
%       cfg = hashtoconfig('a1b2c3d4', config);
%
%   DEPENDENCIES: none (reads config_registry.json written by pipeline_simulate)
%
%   See also: compute_sim_config_hash, load_recon_dose_data, pipeline_simulate

    if ~ischar(hash) && ~isstring(hash)
        error('hashtoconfig:InvalidInput', 'hash must be a string.');
    end
    hash = lower(char(hash));

    if ~isstruct(config) || ~isfield(config, 'working_dir')
        error('hashtoconfig:MissingConfig', ...
            'config must be a struct with a .working_dir field.');
    end

    sim_root = fullfile(config.working_dir, 'SimulationResults');
    if ~exist(sim_root, 'dir')
        error('hashtoconfig:NoSimResults', ...
            'SimulationResults directory not found: %s', sim_root);
    end

    % Registries live at SimulationResults/[patient]/[session]/[method]/config_registry.json
    reg_files = dir(fullfile(sim_root, '**', 'config_registry.json'));
    if isempty(reg_files)
        error('hashtoconfig:NoRegistry', ...
            'No config_registry.json files found under %s', sim_root);
    end

    key = ['h_' hash];  % matches pipeline_simulate's registry_key()
    config_block = [];
    for i = 1:numel(reg_files)
        reg_path = fullfile(reg_files(i).folder, reg_files(i).name);
        registry = read_registry(reg_path);
        if isfield(registry, key) && isfield(registry.(key), 'config')
            config_block = registry.(key).config;
            fprintf('[hashtoconfig] hash %s found in %s\n', hash, reg_path);
            return;
        end
    end

    error('hashtoconfig:HashNotFound', ...
        'Config hash %s not found in any registry under %s', hash, sim_root);
end


function registry = read_registry(reg_path)
    % Read one config_registry.json into a struct; empty struct on any failure.
    registry = struct();
    fid = fopen(reg_path, 'r');
    if fid < 0, return; end
    raw = fread(fid, '*char')';
    fclose(fid);
    if isempty(strtrim(raw)), return; end
    try
        decoded = jsondecode(raw);
        if isstruct(decoded)
            registry = decoded;
        end
    catch
        % Skip unreadable/corrupt registry files silently; other registries
        % may still hold the hash.
    end
end
