classdef StudyRunner < handle
%STUDYRUNNER  Source-agnostic study orchestration over SimResults.
%
%   PURPOSE:
%   Runs the studies (compare / sweep / per-dose gamma batch) and the metric
%   methods (gamma, ssim, gammaSweep) on top of ethos.Simulator. Every study
%   first turns SimCases into SimResults via the Simulator's load-or-simulate
%   resolve(), then analyzes the SimResults -- so it never branches on whether a
%   dose was simulated or loaded. gamma() defaults to cfg.gamma_default (3%/3mm);
%   ssim() scores over the same 10%-of-truth mask so the two are comparable.
%
%   USAGE:
%       sr = ethos.StudyRunner(ethos.SimConfig.gammaBatch().withOverrides( ...
%               'working_dir', WD, 'patient_id', '1194203', 'session', 'Session_1'));
%       r  = sr.resolve(caseCT1);
%       G  = sr.gamma(r);                 % recon vs RayStation truth, 3%/3mm
%       S  = sr.ssim(r);                  % SSIM over the gamma mask
%       C  = sr.compare(caseCT1, caseCT3);% adapted-vs-reference change detection
%       W  = sr.sweep('num_time_reversal_iter', 1:10, cases);
%       B  = sr.gammaBatch([15 112; 3 44]);
%
%   See also: ethos.Simulator, ethos.Analysis, ethos.SimPlotter, ethos.SimCase

    properties
        cfg           % ethos.SimConfig
        sim           % ethos.Simulator
    end

    methods
        function obj = StudyRunner(cfg)
            if nargin < 1 || isempty(cfg)
                cfg = ethos.SimConfig();
            end
            obj.cfg = cfg;
            obj.sim = ethos.Simulator(cfg);
        end

        function r = resolve(obj, simCase)
        %RESOLVE  Convenience pass-through to the Simulator's load-or-simulate.
            r = obj.sim.resolve(simCase);
        end

        % ---- Metric methods -------------------------------------------------

        function G = gamma(obj, result, varargin)
        %GAMMA  Gamma of recon vs truth, evaluated over the truth's 10% mask.
        %   Options: 'Reference' (default 'recon_dose'), 'Target' (default
        %   'rs_truth'), 'MaskFrom' (default = Target), 'Criteria' (default
        %   cfg.gamma_default), 'Cutoff' (default 0.10).
            p = obj.parseMetricOpts(varargin);
            refVol  = obj.getVol(result, p.Reference);
            tgtVol  = obj.getVol(result, p.Target);
            maskVol = obj.getVol(result, p.MaskFrom);
            mask    = maskVol >= p.Cutoff * max(maskVol(:));
            G = ethos.Analysis.gamma(refVol, tgtVol, result.spacing(), ...
                'DosePct', p.Criteria(1), 'DistMm', p.Criteria(2), 'EvalMask', mask);
        end

        function G = gammaChange(obj, resultA, resultB, varargin)
        %GAMMACHANGE  Gamma of reconA vs reconB (change detection), mask = A truth.
            p = obj.parseMetricOpts(varargin);
            refVol  = resultA.recon_dose;
            tgtVol  = resultB.recon_dose;
            maskVol = obj.getVol(resultA, p.MaskFrom);
            mask    = maskVol >= p.Cutoff * max(maskVol(:));
            G = ethos.Analysis.gamma(refVol, tgtVol, resultA.spacing(), ...
                'DosePct', p.Criteria(1), 'DistMm', p.Criteria(2), 'EvalMask', mask);
        end

        function S = ssim(obj, result, varargin)
        %SSIM  SSIM over the same 10%-of-truth mask gamma uses.
            p = obj.parseMetricOpts(varargin);
            reconVol = obj.getVol(result, p.Reference);
            truthVol = obj.getVol(result, p.Target);
            mask     = truthVol >= p.Cutoff * max(truthVol(:));
            S = ethos.Analysis.ssim(truthVol, reconVol, 'EvalMask', mask);
        end

        function SW = gammaSweep(obj, result, crit_vec, varargin)
        %GAMMASWEEP  Pass-rate vs n%/n mm curve for recon vs truth.
            if nargin < 3 || isempty(crit_vec), crit_vec = 1:0.5:5; end
            p = obj.parseMetricOpts(varargin);
            refVol = obj.getVol(result, p.Reference);
            tgtVol = obj.getVol(result, p.Target);
            SW = ethos.Analysis.gammaSweep(refVol, tgtVol, result.spacing(), ...
                crit_vec, p.Cutoff);
        end

        % ---- Studies --------------------------------------------------------

        function C = compare(obj, caseA, caseB, varargin)
        %COMPARE  Adapted-vs-reference change detection between two cases.
        %   Resolves both (load-or-sim), computes gammaChange + per-case gamma
        %   and ssim, and plots when cfg.plot_results is on.
            rA = obj.resolveSensor(caseA);
            rB = obj.resolveSensor(caseB);
            C = struct();
            C.resultA     = rA;
            C.resultB     = rB;
            C.gamma_change = obj.gammaChange(rA, rB, varargin{:});
            C.gamma_A     = obj.gamma(rA, varargin{:});
            C.ssim_A      = obj.ssim(rA, varargin{:});
            fprintf('[compare] %s vs %s: change-gamma pass %.1f%% (%g%%/%gmm)\n', ...
                caseA.label(), caseB.label(), C.gamma_change.pass_rate, ...
                C.gamma_change.criteria(1), C.gamma_change.criteria(2));
            if obj.cfg.plot_results
                obj.plotResult(rA, C.gamma_A);
                ethos.SimPlotter.changePanels(rA, rB, C.gamma_change, ...
                    obj.cfg.viz_smooth_sigma);
            end
        end

        function C = compareBlind(obj, fwdCase, refCase, varargin)
        %COMPAREBLIND  Reference full recon vs blind recon of the adapted anatomy.
        %   refCase (e.g. CT_1) is resolved normally; fwdCase (e.g. CT_3) is
        %   forward-propagated on its own geometry but time-reversed on refCase's
        %   geometry (Simulator.runBlind). Gamma-compares the two recons on the
        %   shared reference grid -- the live-IRAI change-detection task.
            rRef   = obj.resolve(refCase);
            rRef   = obj.ensureSensor(rRef, refCase);
            rBlind = obj.sim.runBlind(fwdCase, refCase);
            C = struct();
            C.reference_result = rRef;
            C.blind_result     = rBlind;
            C.gamma_change     = obj.gammaChange(rRef, rBlind, varargin{:});
            fprintf('[compareBlind] %s (ref) vs %s (blind): change-gamma %.1f%%\n', ...
                refCase.label(), fwdCase.label(), C.gamma_change.pass_rate);
            if obj.cfg.plot_results
                ethos.SimPlotter.changePanels(rRef, rBlind, C.gamma_change, ...
                    obj.cfg.viz_smooth_sigma);
                ethos.SimPlotter.sensorPlacement(rRef);
            end
        end

        function W = sweep(obj, param_name, values, cases)
        %SWEEP  Vary one CONFIG parameter; report mean gamma pass rate per value.
        %   Resolves `cases` under each override (Simulate by default so the
        %   swept parameter actually changes the recon) and averages the recon-
        %   vs-truth pass rate across cases.
            values = values(:)';
            nV = numel(values);
            pass_mean = nan(1, nV);
            pass_all  = nan(numel(cases), nV);
            for v = 1:nV
                cfgv = obj.cfg.withOverrides(param_name, values(v));
                srv  = ethos.StudyRunner(cfgv);
                for c = 1:numel(cases)
                    r = srv.resolve(cases(c));
                    G = srv.gamma(r);
                    pass_all(c, v) = G.pass_rate;
                end
                pass_mean(v) = mean(pass_all(:, v), 'omitnan');
                fprintf('[sweep] %s = %g -> mean pass %.1f%%\n', ...
                    param_name, values(v), pass_mean(v));
            end
            W = struct('param', param_name, 'values', values, ...
                'pass_mean', pass_mean, 'pass_all', pass_all, ...
                'criteria', obj.cfg.gamma_default);
            if obj.cfg.plot_results
                obj.plotSweep(W);
            end
        end

        function W = paramSweep(obj, simCase, param, values, varargin)
        %PARAMSWEEP  Vary one axis over a SINGLE field; keep every recon + metric.
        %   For each value the case is re-simulated (dose_source forced to
        %   Simulate so the swept parameter actually changes the recon), then
        %   scored. Returns per-value SimResults so any figure can be drawn
        %   afterwards; overlays pass rate + SSIM vs value when plotting.
        %
        %   Name/Value:
        %     'Apply'   fn @(cfg,val)->cfg - how to fold a value into the config.
        %               Default: cfg.withOverrides(param, val). Use this for axes
        %               that are not scalar fields, e.g.
        %               @(c,v) c.withTissue('threshold_2','air','sound_speed',v).
        %     'ScoreVs' 'truth' (default) | 'first' | a SimResult.
        %               'truth' scores each recon vs its RayStation truth; 'first'
        %               scores every run vs the FIRST value's recon (a
        %               self-convergence metric, e.g. dt/Nt sweeps); a SimResult
        %               scores against that fixed reference run.
        %     'Label'   x-axis label (default = param).
        %     'Plot'    logical (default cfg.plot_results).
            p = obj.parseSweepOpts(param, varargin);
            values = values(:)';
            nV = numel(values);
            runs = repmat(struct('value', [], 'result', ethos.SimResult(), ...
                'gamma', struct(), 'ssim', struct()), 1, nV);
            pass = nan(1, nV); ssimv = nan(1, nV);

            % Pass 1: simulate every value (keep the results).
            for v = 1:nV
                cfgv = p.Apply(obj.cfg, values(v));
                cfgv = cfgv.withOverrides('dose_source', ethos.DoseSource.Simulate, ...
                    'plot_results', false);
                runs(v).value  = values(v);
                runs(v).result = ethos.StudyRunner(cfgv).resolve(simCase);
            end

            % Resolve the scoring reference.
            refResult = [];
            if ischar(p.ScoreVs) || isstring(p.ScoreVs)
                if strcmpi(char(p.ScoreVs), 'first'), refResult = runs(1).result; end
            elseif isa(p.ScoreVs, 'ethos.SimResult')
                refResult = p.ScoreVs;
            end

            % Pass 2: score each run (vs truth, or vs the reference run).
            for v = 1:nV
                [G, S] = obj.scorePair(refResult, runs(v).result);
                runs(v).gamma = G; runs(v).ssim = S;
                pass(v) = G.pass_rate; ssimv(v) = S.ssim_index;
                fprintf('[paramSweep] %s = %g -> pass %.1f%%  ssim %.3f\n', ...
                    p.Label, values(v), pass(v), ssimv(v));
            end

            W = struct('param', param, 'label', p.Label, 'values', values, ...
                'pass', pass, 'ssim', ssimv, 'runs', runs, ...
                'criteria', obj.cfg.gamma_default);
            if p.Plot
                obj.plotParamSweep(W);
            end
        end

        function B = gammaBatch(obj, specs, varargin)
        %GAMMABATCH  Per-dose A1 (recon vs truth) + A2 (CT_1 vs CT_3) gamma.
        %   specs: Nx2 [beam segment] matrix, or a struct array with .beam/.segment.
        %   For each spec both CT labels are resolved (load-or-sim). A1 = recon_CT1
        %   vs its RayStation truth; A2 = recon_CT1 vs recon_CT3 (change).
            specs = obj.normalizeSpecs(specs);
            n = numel(specs);
            rows = repmat(struct('beam', [], 'segment', [], ...
                'A1_pass', NaN, 'A2_pass', NaN, 'A1_ssim', NaN), 1, n);
            for i = 1:n
                b = specs(i).beam; s = specs(i).segment;
                cCT1 = ethos.SimCase.list(obj.cfg.patient_id, obj.cfg.session, ...
                    obj.cfg, 'Beam', b, 'Segment', s, 'CTLabel', 'CT_1');
                cCT3 = ethos.SimCase.list(obj.cfg.patient_id, obj.cfg.session, ...
                    obj.cfg, 'Beam', b, 'Segment', s, 'CTLabel', 'CT_3');
                if isempty(cCT1)
                    warning('ethos:StudyRunner:NoField', ...
                        'No CT_1 field for beam %g seg %g; skipping.', b, s);
                    continue;
                end
                rCT1 = obj.resolveSensor(cCT1(1));
                A1 = obj.gamma(rCT1, varargin{:});
                S1 = obj.ssim(rCT1, varargin{:});
                rows(i).beam = b; rows(i).segment = s;
                rows(i).A1_pass = A1.pass_rate;
                rows(i).A1_ssim = S1.ssim_index;
                if ~isempty(cCT3)
                    rCT3 = obj.resolve(cCT3(1));
                    A2 = obj.gammaChange(rCT1, rCT3, varargin{:});
                    rows(i).A2_pass = A2.pass_rate;
                end
                fprintf('[gammaBatch] B%g_%g: A1 %.1f%%  A2 %.1f%%  ssim %.3f\n', ...
                    b, s, rows(i).A1_pass, rows(i).A2_pass, rows(i).A1_ssim);
                if obj.cfg.plot_results
                    obj.plotResult(rCT1, A1);
                end
            end
            B = struct('rows', rows, 'criteria', obj.cfg.gamma_default, ...
                'patient_id', obj.cfg.patient_id, 'session', obj.cfg.session);
        end

        function BS = beamSummary(obj, varargin)
        %BEAMSUMMARY  Per-beam gamma/SSIM over ALL segments, both CT images.
        %   Load-only companion to gammaBatch: for every CT_1 field with a CT_3
        %   twin it scores four comparisons per segment, then reduces to a
        %   per-beam mean +/- std. Missing/mismatched entries are skipped.
        %   Options: 'Metric' 'gamma'(default)|'ssim', 'PlanType' (default 'reference').
        %   Comparisons (ref vs tgt, mask from ref's 10%):
        %     t1t3 truth_CT1 vs truth_CT3 | t1r1 truth_CT1 vs recon_CT1
        %     t1r3 truth_CT1 vs recon_CT3 | t3r3 truth_CT3 vs recon_CT3
            p = struct('Metric', 'gamma', 'PlanType', 'reference');
            for k = 1:2:numel(varargin), p.(varargin{k}) = varargin{k + 1}; end

            c1 = ethos.SimCase.list(obj.cfg.patient_id, obj.cfg.session, ...
                obj.cfg, 'PlanType', p.PlanType, 'CTLabel', 'CT_1');
            beams = unique([c1.beam]);
            nB = numel(beams);
            keys = {'t1t3', 't1r1', 't1r3', 't3r3'};
            acc = struct();
            for kk = 1:numel(keys), acc.(keys{kk}) = cell(1, nB); end

            for i = 1:numel(c1)
                bi = find(beams == c1(i).beam, 1);
                c3 = ethos.SimCase.list(obj.cfg.patient_id, obj.cfg.session, ...
                    obj.cfg, 'Beam', c1(i).beam, 'Segment', c1(i).segment, ...
                    'PlanType', p.PlanType, 'CTLabel', 'CT_3');
                if isempty(c3), continue; end
                try
                    r1 = obj.resolve(c1(i));
                    r3 = obj.resolve(c3(1));
                    sp = r1.spacing();
                    acc.t1t3{bi}(end+1) = obj.metricPass(p.Metric, r1.rs_truth,   r3.rs_truth,   sp);
                    acc.t1r1{bi}(end+1) = obj.metricPass(p.Metric, r1.rs_truth,   r1.recon_dose, sp);
                    acc.t1r3{bi}(end+1) = obj.metricPass(p.Metric, r1.rs_truth,   r3.recon_dose, sp);
                    acc.t3r3{bi}(end+1) = obj.metricPass(p.Metric, r3.rs_truth,   r3.recon_dose, sp);
                catch ME
                    warning('ethos:StudyRunner:SegSkip', ...
                        'Beam %g seg %g skipped: %s', c1(i).beam, c1(i).segment, ME.message);
                end
            end

            BS = struct('beams', beams, 'metric', p.Metric, ...
                'criteria', obj.cfg.gamma_default);
            for kk = 1:numel(keys)
                key = keys{kk};
                BS.([key '_mean']) = cellfun(@(v) mean(v, 'omitnan'), acc.(key));
                BS.([key '_std'])  = cellfun(@(v) std(v,  'omitnan'), acc.(key));
            end
            if obj.cfg.plot_results
                obj.plotBeamSummary(BS);
            end
        end

        % ---- Sensor / plotting ---------------------------------------------

        function result = ensureSensor(obj, result, simCase)
        %ENSURESENSOR  Attach a sensor mask (+ density) to a loaded result.
        %   Loaded recons carry no sensor mask; rebuild it from the case geometry
        %   so sensor-placement plots work in analyze-only studies.
            if ~isempty(result.sensor_mask)
                return;
            end
            sc = simCase.materialize(obj.cfg);
            if ~isempty(sc.precomputed_sensor) && isfield(sc.precomputed_sensor, 'sensor_mask')
                result = result.withSensor(sc.precomputed_sensor.sensor_mask, ...
                    sc.precomputed_sensor.sensor_info);
            end
            if isempty(result.density) && isfield(sc.cbct_resampled, 'cubeDensity')
                result.density = sc.cbct_resampled.cubeDensity;
            end
        end

        function plotResult(obj, result, G)
        %PLOTRESULT  Dispatch the per-result figures gated by cfg.enable.
            if nargin < 3, G = []; end
            en = obj.cfg.enable;
            if isfield(en, 'reconVsTruth') && en.reconVsTruth
                ethos.SimPlotter.reconVsTruth(result, G, ...
                    {'RayStation truth', 'Recon'}, obj.cfg.viz_smooth_sigma);
            end
            if isfield(en, 'sensorPlacement') && en.sensorPlacement
                ethos.SimPlotter.sensorPlacement(result);
            end
            if isfield(en, 'convergence') && en.convergence && result.isSimulated()
                ethos.SimPlotter.convergence(result);
            end
        end

        function plotSweep(~, W)
        %PLOTSWEEP  Mean gamma pass rate vs the swept parameter value.
            figure('Name', 'Parameter Sweep', 'Color', 'w', 'NumberTitle', 'off');
            plot(W.values, W.pass_mean, '-o', 'LineWidth', 1.8, 'MarkerSize', 6, ...
                'MarkerFaceColor', [0.2 0.4 0.8]);
            grid on; box on; ylim([0 105]);
            yline(90, 'r--', '90%');
            xlabel(strrep(W.param, '_', '\_'));
            ylabel(sprintf('Mean gamma pass rate (%%)  @ %g%%/%g mm', ...
                W.criteria(1), W.criteria(2)));
            title(sprintf('Sweep: %s', strrep(W.param, '_', '\_')));
            drawnow;
        end

        function plotBeamSummary(~, BS)
        %PLOTBEAMSUMMARY  Two per-beam figures: change detection + own-CT fidelity.
            b = BS.beams(:)';
            ylab = 'Gamma pass (%)';
            if strcmpi(BS.metric, 'ssim'), ylab = 'SSIM (%)'; end

            figure('Name', 'Change Detection', 'Color', 'w', 'NumberTitle', 'off');
            hold on;
            errorbar(b, BS.t1r1_mean, BS.t1r1_std, 'o-', 'Color', [0.15 0.6 0.2], ...
                'DisplayName', 'recon CT\_1 vs truth CT\_1');
            errorbar(b, BS.t1r3_mean, BS.t1r3_std, 'o-', 'Color', [0.8 0.15 0.15], ...
                'DisplayName', 'recon CT\_3 vs truth CT\_1');
            errorbar(b, BS.t1t3_mean, BS.t1t3_std, 'o-', 'Color', [0.2 0.4 0.8], ...
                'DisplayName', 'truth CT\_1 vs truth CT\_3');
            yline(90, 'k--', '90%'); grid on; box on; ylim([0 105]);
            xlabel('Beam number'); ylabel(ylab); legend('Location', 'best');
            title('Change detection (mean \pm std over segments)');

            figure('Name', 'Reconstruction Fidelity', 'Color', 'w', 'NumberTitle', 'off');
            hold on;
            errorbar(b, BS.t1r1_mean, BS.t1r1_std, 'o-', 'Color', [0.15 0.6 0.2], ...
                'DisplayName', 'truth CT\_1 vs recon CT\_1');
            errorbar(b, BS.t3r3_mean, BS.t3r3_std, 'o-', 'Color', [0.8 0.15 0.15], ...
                'DisplayName', 'truth CT\_3 vs recon CT\_3');
            yline(90, 'k--', '90%'); grid on; box on; ylim([0 105]);
            xlabel('Beam number'); ylabel(ylab); legend('Location', 'best');
            title('Reconstruction fidelity (mean \pm std over segments)');
            drawnow;
        end

        function plotParamSweep(~, W)
        %PLOTPARAMSWEEP  Pass rate (left axis) + SSIM (right axis) vs swept value.
            figure('Name', 'Parameter Sweep', 'Color', 'w', 'NumberTitle', 'off');
            yyaxis left;
            plot(W.values, W.pass, '-o', 'LineWidth', 1.8, 'MarkerSize', 6);
            ylim([0 105]); ylabel(sprintf('Gamma pass (%%)  @ %g%%/%g mm', ...
                W.criteria(1), W.criteria(2)));
            yline(90, 'k--', '90%');
            yyaxis right;
            plot(W.values, W.ssim, '-s', 'LineWidth', 1.4, 'MarkerSize', 6);
            ylabel('SSIM (masked)');
            grid on; box on;
            xlabel(strrep(W.label, '_', '\_'));
            title(sprintf('Sweep: %s', strrep(W.label, '_', '\_')));
            drawnow;
        end

        function plotPassRateCurve(~, SW, label_str)
        %PLOTPASSRATECURVE  Pass-rate vs n%/n mm from a gammaSweep result.
            if nargin < 3, label_str = ''; end
            figure('Name', 'Gamma Pass Rate vs Criterion', 'Color', 'w', ...
                'NumberTitle', 'off');
            plot(SW.crit, SW.pass_rates, '-o', 'LineWidth', 1.8, 'MarkerSize', 6, ...
                'MarkerFaceColor', [0.2 0.4 0.8]);
            grid on; box on; ylim([0 105]);
            yline(90, 'r--', '90%');
            xlabel('Gamma criterion  n (%/mm)');
            ylabel('Pass rate (%)');
            title(sprintf('Gamma pass rate vs criterion  %s', ...
                strrep(label_str, '_', '\_')));
            drawnow;
        end
    end

    methods (Access = private)
        function r = resolveSensor(obj, simCase)
        %RESOLVESENSOR  Resolve then ensure a sensor mask for plotting.
            r = obj.resolve(simCase);
            r = obj.ensureSensor(r, simCase);
        end

        function p = parseMetricOpts(obj, args)
            p = struct('Reference', 'recon_dose', 'Target', 'rs_truth', ...
                'MaskFrom', '', 'Criteria', obj.cfg.gamma_default, 'Cutoff', 0.10);
            if mod(numel(args), 2) ~= 0
                error('ethos:StudyRunner:BadOptions', 'Options must be Name,Value pairs.');
            end
            for k = 1:2:numel(args)
                name = char(args{k});
                p.(name) = args{k + 1};
            end
            if isempty(p.MaskFrom), p.MaskFrom = p.Target; end
        end

        function v = getVol(~, result, name)
            switch lower(name)
                case 'recon_dose',  v = result.recon_dose;
                case 'rs_truth',    v = result.rs_truth;
                case 'ethos_truth', v = result.ethos_truth;
                otherwise
                    if isprop(result, name)
                        v = result.(name);
                    else
                        error('ethos:StudyRunner:BadVol', ...
                            'Unknown volume selector ''%s''.', name);
                    end
            end
            if isempty(v)
                error('ethos:StudyRunner:EmptyVol', ...
                    'Result has no ''%s'' volume for this analysis.', name);
            end
            v = double(v);
        end

        function p = parseSweepOpts(obj, param, args)
            p = struct('Apply', @(c, v) c.withOverrides(param, v), ...
                'Label', param, 'Plot', obj.cfg.plot_results, 'ScoreVs', 'truth');
            if mod(numel(args), 2) ~= 0
                error('ethos:StudyRunner:BadOptions', 'Options must be Name,Value pairs.');
            end
            for k = 1:2:numel(args)
                name = char(args{k});
                switch lower(name)
                    case 'apply',   p.Apply   = args{k + 1};
                    case 'label',   p.Label   = args{k + 1};
                    case 'plot',    p.Plot    = logical(args{k + 1});
                    case 'scorevs', p.ScoreVs = args{k + 1};
                    otherwise
                        error('ethos:StudyRunner:UnknownOption', ...
                            'Unknown paramSweep option ''%s''.', name);
                end
            end
        end

        function val = metricPass(obj, metric, refVol, tgtVol, spacing_mm)
        %METRICPASS  Scalar score (%) of tgt vs ref over ref's 10% mask.
            refVol = double(refVol); tgtVol = double(tgtVol);
            mask = refVol >= 0.10 * max(refVol(:));
            crit = obj.cfg.gamma_default;
            if strcmpi(metric, 'ssim')
                S = ethos.Analysis.ssim(refVol, tgtVol, 'EvalMask', mask);
                val = 100 * S.ssim_index;
            else
                G = ethos.Analysis.gamma(refVol, tgtVol, spacing_mm, ...
                    'DosePct', crit(1), 'DistMm', crit(2), 'EvalMask', mask);
                val = G.pass_rate;
            end
        end

        function [G, S] = scorePair(obj, refResult, r)
        %SCOREPAIR  Gamma+SSIM of r vs its truth (refResult empty) or vs refResult.
            if isempty(refResult)
                G = obj.gamma(r);
                S = obj.ssim(r);
            else
                maskVol = refResult.rs_truth;
                if isempty(maskVol), maskVol = refResult.recon_dose; end
                mask = maskVol >= 0.10 * max(maskVol(:));
                crit = obj.cfg.gamma_default;
                G = ethos.Analysis.gamma(refResult.recon_dose, r.recon_dose, ...
                    r.spacing(), 'DosePct', crit(1), 'DistMm', crit(2), 'EvalMask', mask);
                S = ethos.Analysis.ssim(refResult.recon_dose, r.recon_dose, 'EvalMask', mask);
            end
        end

        function specs = normalizeSpecs(~, specs)
            if isnumeric(specs)
                if size(specs, 2) ~= 2
                    error('ethos:StudyRunner:BadSpecs', ...
                        'Numeric specs must be an Nx2 [beam segment] matrix.');
                end
                s = repmat(struct('beam', [], 'segment', []), 1, size(specs, 1));
                for i = 1:size(specs, 1)
                    s(i).beam = specs(i, 1);
                    s(i).segment = specs(i, 2);
                end
                specs = s;
            end
        end
    end
end
