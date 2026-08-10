classdef Analysis
%ANALYSIS  Static gamma + SSIM primitives over two co-gridded dose volumes.
%
%   PURPOSE:
%   The metric core shared by every study. Operates on plain 3D arrays so it is
%   independent of where the doses came from. gamma() defaults to 3%/3mm; ssim()
%   scores over the SAME 10%-of-reference evaluation mask that gamma uses, so
%   the two metrics are computed over an identical voxel set and are directly
%   comparable per study.
%
%   Gamma matches the pipeline convention exactly: global gamma
%   (CalcGamma(...,'local',0)), distance search capped at 2x the DTA criterion,
%   'restrict',1 for speed. Pass rate = 100 * mean(gamma(mask) <= 1).
%
%   See also: CalcGamma, ssim, ethos.StudyRunner, ethos.SimResult

    methods (Static)
        function G = gamma(ref, tgt, spacing_mm, varargin)
        %GAMMA  Gamma index map + pass rate between reference and target volumes.
        %   G = ethos.Analysis.gamma(ref, tgt, spacing_mm, Name, Value, ...)
        %   Options:
        %     'DosePct'  (default 3)   - dose-difference criterion (%).
        %     'DistMm'   (default 3)   - distance-to-agreement criterion (mm).
        %     'Cutoff'   (default 0.10)- low-dose cutoff as a fraction of max(ref).
        %     'EvalMask' (default [])  - explicit logical mask; overrides Cutoff.
        %   Returns G with: gamma_map, pass_rate, eval_mask, criteria [pct dist],
        %   n_eval, cutoff.
            opts = ethos.Analysis.parseGammaOpts(varargin);
            ref  = double(ref);
            tgt  = double(tgt);

            eval_mask = ethos.Analysis.resolveMask(ref, opts.EvalMask, opts.Cutoff);

            ref_struct = struct('start', [0, 0, 0], 'width', spacing_mm(:)', 'data', ref);
            tgt_struct = struct('start', [0, 0, 0], 'width', spacing_mm(:)', 'data', tgt);

            gmap = CalcGamma(ref_struct, tgt_struct, opts.DosePct, opts.DistMm, ...
                'local', 0, 'limit', 2 * opts.DistMm, 'restrict', 1);

            % Report/plot only inside the eval region (NaN elsewhere).
            gamma_map = nan(size(gmap));
            gamma_map(eval_mask) = gmap(eval_mask);

            G = struct();
            G.gamma_map = gamma_map;
            G.eval_mask = eval_mask;
            G.criteria  = [opts.DosePct, opts.DistMm];
            G.cutoff    = opts.Cutoff;
            G.n_eval    = nnz(eval_mask);
            if G.n_eval > 0
                G.pass_rate = 100 * mean(gmap(eval_mask) <= 1);
            else
                G.pass_rate = NaN;
            end
        end

        function S = ssim(ref, tgt, varargin)
        %SSIM  Structural similarity scored over the gamma eval mask.
        %   S = ethos.Analysis.ssim(ref, tgt, Name, Value, ...)
        %   Options: 'Cutoff' (default 0.10) and 'EvalMask' (overrides Cutoff),
        %   matching ethos.Analysis.gamma so both metrics score the same region.
        %   Returns S with: ssim_index (mask-mean), ssim_global, ssim_map,
        %   eval_mask, cutoff.
            opts = ethos.Analysis.parseGammaOpts(varargin);
            ref  = double(ref);
            tgt  = double(tgt);

            eval_mask = ethos.Analysis.resolveMask(ref, opts.EvalMask, opts.Cutoff);

            [ssim_global, ssim_map] = ssim(tgt, ref);

            S = struct();
            S.ssim_map    = ssim_map;
            S.ssim_global = ssim_global;
            S.eval_mask   = eval_mask;
            S.cutoff      = opts.Cutoff;
            if any(eval_mask(:))
                S.ssim_index = mean(ssim_map(eval_mask));
            else
                S.ssim_index = NaN;
            end
        end

        function SW = gammaSweep(ref, tgt, spacing_mm, crit_vec, cutoff)
        %GAMMASWEEP  Pass rate for each n%/n mm criterion (mirrors gamma_sweep_pass_rates).
        %   Reference = ref, target = tgt, 10%-of-ref cutoff (override via cutoff).
            if nargin < 5 || isempty(cutoff), cutoff = 0.10; end
            ref = double(ref);
            tgt = double(tgt);
            ref_struct = struct('start', [0, 0, 0], 'width', spacing_mm(:)', 'data', ref);
            tgt_struct = struct('start', [0, 0, 0], 'width', spacing_mm(:)', 'data', tgt);
            eval_mask  = ref >= cutoff * max(ref(:));

            pr = nan(numel(crit_vec), 1);
            for k = 1:numel(crit_vec)
                c = crit_vec(k);
                try
                    gmap  = CalcGamma(ref_struct, tgt_struct, c, c, ...
                        'local', 0, 'limit', c * 2, 'restrict', 1);
                    pr(k) = 100 * mean(gmap(eval_mask) <= 1);
                catch
                    pr(k) = NaN;
                end
            end

            SW = struct('crit', crit_vec(:), 'pass_rates', pr, 'eval_mask', eval_mask);
        end
    end

    methods (Static, Access = private)
        function opts = parseGammaOpts(args)
            opts = struct('DosePct', 3, 'DistMm', 3, 'Cutoff', 0.10, 'EvalMask', []);
            if mod(numel(args), 2) ~= 0
                error('ethos:Analysis:BadOptions', 'Options must be Name,Value pairs.');
            end
            for k = 1:2:numel(args)
                name = lower(char(args{k}));
                val  = args{k + 1};
                switch name
                    case 'dosepct',  opts.DosePct  = val;
                    case 'distmm',   opts.DistMm   = val;
                    case 'cutoff',   opts.Cutoff   = val;
                    case 'evalmask', opts.EvalMask = logical(val);
                    case 'criteria'  % [pct dist] shorthand
                        opts.DosePct = val(1);
                        opts.DistMm  = val(2);
                    otherwise
                        error('ethos:Analysis:UnknownOption', 'Unknown option ''%s''.', name);
                end
            end
        end

        function mask = resolveMask(ref, explicitMask, cutoff)
            if ~isempty(explicitMask)
                mask = logical(explicitMask);
            else
                mask = ref >= cutoff * max(ref(:));
            end
        end
    end
end
