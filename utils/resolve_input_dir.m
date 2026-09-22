function scan_dir = resolve_input_dir(primary_dir, fallback_dir, patterns)
%RESOLVE_INPUT_DIR Choose an input directory: primary first, fallback second.
%
%   scan_dir = resolve_input_dir(primary_dir, fallback_dir, patterns)
%
%   PURPOSE:
%   Lets the compress pipeline scan RayStationFiles as the primary input
%   location while falling back to the native EthosExports copy when
%   RayStationFiles holds none of the requested files. Processed outputs are
%   written elsewhere (RayStationFiles/.../processed) regardless of which
%   directory the inputs came from — this only picks where to READ.
%
%   INPUTS:
%       primary_dir  - char, preferred directory (RayStationFiles/<id>/<session>)
%       fallback_dir - char, fallback directory (EthosExports/.../<session>/sct)
%       patterns     - cellstr (or char) of dir() glob patterns to look for,
%                      e.g. {'dose_*.npz'} or {'dose_*.mat','dose_*.dcm'}
%
%   OUTPUT:
%       scan_dir - primary_dir if it contains any file matching one of the
%                  patterns; else fallback_dir if IT contains a match; else
%                  primary_dir unchanged (so the caller reports its usual
%                  "not found" behavior against the expected location).
%
%   See also: pipeline_compress, step14_npz_to_mat, step15_process_doses

    if ischar(patterns) || isstring(patterns)
        patterns = {char(patterns)};
    end

    if dir_has_any(primary_dir, patterns)
        scan_dir = primary_dir;
        return;
    end

    if dir_has_any(fallback_dir, patterns)
        scan_dir = fallback_dir;
        fprintf('  [FALLBACK] No matching input in %s\n', primary_dir);
        fprintf('             Using native EthosExports copy: %s\n', fallback_dir);
        return;
    end

    scan_dir = primary_dir;
end


function tf = dir_has_any(d, patterns)
%DIR_HAS_ANY True if directory d exists and holds a file matching a pattern.
    tf = false;
    if isempty(d) || ~isfolder(d)
        return;
    end
    for k = 1:numel(patterns)
        if ~isempty(dir(fullfile(d, patterns{k})))
            tf = true;
            return;
        end
    end
end
