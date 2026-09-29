function done_beams = step25_metrics_watcher(patient_id, session, config, stop_file)
%STEP25_METRICS_WATCHER Run Step 2.5 on each beam as soon as Step 2 finishes it.
%
%   done_beams = step25_metrics_watcher(patient_id, session, config, stop_file)
%
%   PURPOSE:
%   Lets the CPU-bound Step 2.5 (per-segment gamma + SSIM) run AT THE SAME TIME
%   as the GPU-bound Step 2 (k-Wave), instead of after it. pipeline_simulate
%   starts this function as a background batch job with its own pool of CPU
%   workers, then runs the Step 2 GPU parfor as usual. The two sides never talk
%   directly: Step 2 writes recon files to disk, and this loop picks them up.
%
%   Every poll it looks for beams whose recon files ALL exist on disk (and have
%   not changed for 30 s, so none is still being saved) and that it has not
%   scored yet, and runs step25_segment_metrics on just those beams
%   (step25 loads a whole beam at once, so a beam is the smallest unit of work).
%   It stops once every beam is scored, or once pipeline_simulate creates
%   stop_file (after one final pass, so the last beams Step 2 wrote are caught).
%
%   Best effort: a beam that fails here is not retried. The final Step 2.5 call
%   in pipeline_simulate computes any segment that is still not folded in, and
%   it also computes the noise floor and writes the summary (both are switched
%   off here).
%
%   INPUTS:
%       patient_id - Char/string patient identifier (e.g. '1194203').
%       session    - Char/string session name (e.g. 'Session_1').
%       config     - The pipeline_simulate CONFIG struct. Uses .working_dir,
%                    .gruneisen_method, .metrics_beams ([] => all beams),
%                    .metrics_watcher_poll_sec (seconds between polls), plus
%                    everything step25_segment_metrics reads.
%       stop_file  - Path of the flag file pipeline_simulate creates when
%                    Step 2 is over. Unique per pipeline_simulate process so
%                    sibling instances do not stop each other's watchers.
%
%   OUTPUTS:
%       done_beams - Beam numbers this watcher handed to step25_segment_metrics.
%
%   ALGORITHM:
%       1. List every field (list_processed_field_doses) and its expected recon
%          file for the active config hash, grouped by beam.
%       2. Loop: note whether stop_file exists; find unscored beams whose recon
%          files all exist and are >30 s old; run step25_segment_metrics on
%          them; pause.
%       3. Exit when all beams are scored, or after the pass that saw stop_file.
%
%   EXAMPLE (as launched by pipeline_simulate):
%       job = batch(@step25_metrics_watcher, 0, ...
%           {'1194203', 'Session_1', CONFIG, stop_file}, 'Pool', 23);
%
%   DEPENDENCIES:
%       step25_segment_metrics, list_processed_field_doses,
%       compute_sim_config_hash, Parallel Computing Toolbox (batch pool).
%
%   See also: step25_segment_metrics, pipeline_simulate, batch

    if ~ischar(patient_id) && ~isstring(patient_id)
        error('step25_metrics_watcher:InvalidInput', ...
            'patient_id must be a string or character array.');
    end
    patient_id = char(patient_id);
    if ~ischar(session) && ~isstring(session)
        error('step25_metrics_watcher:InvalidInput', ...
            'session must be a string or character array.');
    end
    session   = char(session);
    stop_file = char(stop_file);

    % --- 1. Expected recon file of every field, and the beam it belongs to ---
    hash8   = compute_sim_config_hash(config);
    sim_dir = fullfile(config.working_dir, 'SimulationResults', patient_id, ...
        session, config.gruneisen_method);

    [field_index, ~] = list_processed_field_doses(patient_id, session, config);
    field_beams = [field_index.beam_index];
    recon_files = cell(1, numel(field_index));
    for k = 1:numel(field_index)
        [~, stem] = fileparts(field_index(k).source_mat_filename);
        recon_files{k} = fullfile(sim_dir, sprintf('%s_recon_%s.mat', stem, hash8));
    end

    beams = unique(field_beams(~isnan(field_beams)));
    if isfield(config, 'metrics_beams') && ~isempty(config.metrics_beams)
        beams = intersect(beams, config.metrics_beams(:)');
    end

    % Each pass only folds per-segment metrics into the recon files. The noise
    % floor and the summary are left to the final Step 2.5 call.
    config.metrics_noise_floor   = false;
    config.metrics_write_summary = false;

    fprintf('[WATCHER] Watching %d beam(s) for %s / %s (config %s).\n', ...
        numel(beams), patient_id, session, hash8);

    % --- 2. Poll loop ---
    done_beams = [];
    stop_seen  = false;
    while ~stop_seen
        % Check the flag BEFORE scanning, so the pass after the flag appears
        % still picks up the last recon files the GPUs wrote.
        stop_seen = isfile(stop_file);

        ready = [];
        for b = setdiff(beams, done_beams)
            if all(cellfun(@is_finished_file, recon_files(field_beams == b)))
                ready(end+1) = b; %#ok<AGROW>
            end
        end

        if ~isempty(ready)
            fprintf('[WATCHER] %s  scoring beam(s) %s\n', ...
                datestr(now, 'HH:MM:SS'), mat2str(ready));
            config.metrics_beams = ready;
            try
                step25_segment_metrics(patient_id, session, config);
            catch ME
                fprintf(['[WATCHER] Beam(s) %s failed (%s); left for the final ' ...
                    'Step 2.5 call.\n'], mat2str(ready), ME.message);
            end
            done_beams = [done_beams, ready]; %#ok<AGROW>
        end

        % --- 3. Exit conditions ---
        if numel(done_beams) == numel(beams)
            break;
        end
        if ~stop_seen
            pause(config.metrics_watcher_poll_sec);
        end
    end

    fprintf('[WATCHER] Done: scored %d of %d beam(s).\n', numel(done_beams), numel(beams));
end


function tf = is_finished_file(path)
%IS_FINISHED_FILE True if the file exists and has not been modified for 30 s,
%  so a recon file a GPU worker is still saving is never read (or appended to).
    info = dir(path);
    tf = ~isempty(info) && (now - info(1).datenum) * 24 * 3600 > 30;
end
