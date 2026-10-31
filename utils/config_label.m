function label = config_label(config)
%CONFIG_LABEL Readable folder name saying how a CONFIG differs from the default.
%
%   label = config_label(config)
%
%   PURPOSE:
%   AnalysisResults/<patient>/<session>/<label>/ holds the analysis of one
%   simulation config. get_default_config is the canonical calculation, so the
%   label is 'canonical' when every simulation field matches it, and otherwise
%   lists each field that differs as <field>.<value>, joined by '__':
%       'canonical'
%       'reconstruction_method.das'
%       'downscale_factor.2__gruneisen_method.threshold_1'
%   Only the fields in the config hash are compared (compute_sim_config_hash),
%   so one label <-> one hash; paths, plotting and analysis knobs are ignored.
%
%   'unset' means the field was absent from that run's CONFIG. A hash keeps some
%   legacy values as absent (pulse_shape 'gaussian', sensor_side 'right'), so a
%   gaussian-pulse run is labelled 'pulse_shape.unset'.
%
%   INPUTS:
%       config - Full simulation CONFIG, or the canonical block of a hash read
%                back from config_registry.json (hashtoconfig). Both give the
%                same label.
%
%   OUTPUTS:
%       label - Char folder name.
%
%   EXAMPLE:
%       CONFIG = get_default_config('reconstruction_method', 'das');
%       config_label(CONFIG)                       % 'reconstruction_method.das'
%       config_label(hashtoconfig('a1b2c3d4', CONFIG))
%
%   See also: get_default_config, compute_sim_config_hash, hashtoconfig,
%             get_simulation_dir

    [~, run]      = compute_sim_config_hash(config);
    [~, standard] = compute_sim_config_hash(get_default_config());

    names = union(fieldnames(run), fieldnames(standard));   % sorted
    parts = {};
    for i = 1:numel(names)
        name = names{i}
        if isfield(run, name) && isfield(standard, name) && isequal(run.(name), standard.(name))
            continue;
        end
        if isfield(run, name)
            value = run.(name)
        else
            value = '__unset__';
        end
        parts{end + 1} = [name, '.', value_to_text(value)]; %#ok<AGROW>
    end

    if isempty(parts)
        label = 'canonical';
    else
        label = strjoin(parts, '__');
    end
end


function txt = value_to_text(value)
    % Short, filesystem-safe text for one config value.
    if isequal(value, '__unset__')
        txt = 'unset';
    elseif ischar(value) || isstring(value)
        txt = char(value);
    elseif islogical(value) && isscalar(value)
        if value, txt = 'true'; else, txt = 'false'; end
    elseif isnumeric(value) && isempty(value)
        txt = 'none';
    elseif isnumeric(value) || islogical(value)
        txt = strjoin(arrayfun(@(v) sprintf('%g', v), double(value(:)'), ...
            'UniformOutput', false), '-');
    else
        txt = 'custom';
    end
    txt = regexprep(txt, '[^A-Za-z0-9_.+-]', '-');
end
