function append_gamma_log(config, gamma_results, dose_filename)
%APPEND_GAMMA_LOG Append a run's CONFIG + gamma pass rates to a running .mat log.
%   append_gamma_log(config, gamma_results, dose_filename)
%
%   No-op unless config.log_gamma is true and gamma_results has pass_rates.
%   Maintains a per-criterion best record (highest pass rate + the CONFIG that
%   produced it) in config.gamma_log_file (default gamma_log.mat under
%   config.working_dir). Shared copy of the standalone drivers' logging block.

    if ~(isfield(config, 'log_gamma') && config.log_gamma) || ...
            isempty(gamma_results) || ~isfield(gamma_results, 'pass_rates')
        return;
    end

    if isfield(config, 'gamma_log_file') && ~isempty(config.gamma_log_file)
        log_path = config.gamma_log_file;
    else
        log_path = 'gamma_log.mat';
    end
    if isempty(fileparts(log_path))
        log_path = fullfile(config.working_dir, log_path);
    elseif ~isfolder(fileparts(log_path))
        log_path = fullfile(config.working_dir, log_path);
    end

    entry = struct();
    entry.timestamp     = datestr(now, 'yyyy-mm-dd HH:MM:SS'); %#ok<TNOW1,DATST>
    entry.config        = config;
    entry.dose_filename = dose_filename;
    entry.criteria      = gamma_results.criteria(:, 3);
    entry.pass_rates    = gamma_results.pass_rates(:);

    n_crit = numel(entry.pass_rates);
    blank_best = repmat(struct('criterion', '', 'pass_rate', -Inf, ...
        'config', [], 'timestamp', ''), n_crit, 1);
    for gc = 1:n_crit, blank_best(gc).criterion = entry.criteria{gc}; end

    if isfile(log_path)
        L = load(log_path);
        if isfield(L, 'log_entries'), log_entries = L.log_entries; else, log_entries = struct([]); end
        if isfield(L, 'best'), best = L.best; else, best = blank_best; end
    else
        log_entries = struct([]);
        best = blank_best;
    end

    if isempty(log_entries), log_entries = entry; else, log_entries(end+1) = entry; end

    for gc = 1:n_crit
        pr = entry.pass_rates(gc);
        if ~isnan(pr) && pr > best(gc).pass_rate
            best(gc).criterion = entry.criteria{gc};
            best(gc).pass_rate = pr;
            best(gc).config    = config;
            best(gc).timestamp = entry.timestamp;
        end
    end

    save(log_path, 'log_entries', 'best', '-v7.3');
    fprintf('\n[GammaLog] Appended run to %s (%d total entries).\n', log_path, numel(log_entries));
end
