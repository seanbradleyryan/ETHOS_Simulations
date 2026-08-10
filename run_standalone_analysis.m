%% =========================================================================
%  RUN_STANDALONE_ANALYSIS.m
%  Analysis-only standalone pipeline (NO k-Wave), on the +ethos class layer.
%  Loads already-reconstructed doses from disk (DoseSource.Load) and runs the
%  gamma / SSIM analysis + plots. Three modes:
%    'single'     - recon vs its RayStation truth (gamma 3%/3 mm + panels).
%    'comparison' - CT_1 recon vs CT_3 recon (change detection).
%    'sweep'      - comparison PLUS the gamma pass-rate-vs-criterion chart.
%
%  Replaces the former ~1000-line script. The noise-ensemble error bars and the
%  p0 convergence plot need k-Wave and are intentionally not part of this
%  load-only path (see study_pass_rates_individual.m for the noise ensemble).
%  Remote/HIPAA execution only.
%  =========================================================================

clear; clc; close all;

% ---- Configuration ------------------------------------------------------
WORKING_DIR   = '/mnt/weka/home/80030361/ETHOS_Simulations';
PATIENT_ID    = '1194203';
SESSION       = 'Session_1';
ANALYSIS_MODE = 'single';        % 'single' | 'comparison' | 'sweep'
BEAM          = [];              % [] = first matched field
SEGMENT       = [];
CRITERION_N   = (1:0.5:5)';      % pass-rate sweep grid (sweep mode)

cfg = ethos.SimConfig.standalone().withOverrides( ...
    'working_dir', WORKING_DIR, 'patient_id', PATIENT_ID, 'session', SESSION, ...
    'gruneisen_method', 'threshold_2', ...
    'dose_source', 'load', ...        % analysis-only: never runs k-Wave
    'plot_results', true);

sr = ethos.StudyRunner(cfg);

% ---- Helper to pick a field (beam/segment optional) ---------------------
filt = {};
if ~isempty(BEAM),    filt = [filt, {'Beam', BEAM}];       end
if ~isempty(SEGMENT), filt = [filt, {'Segment', SEGMENT}]; end

switch lower(ANALYSIS_MODE)
    case 'single'
        c1 = ethos.SimCase.list(PATIENT_ID, SESSION, cfg, filt{:}, 'CTLabel', 'CT_1');
        assert(~isempty(c1), 'No CT_1 field matched the filter.');
        r  = sr.resolve(c1(1));
        r  = sr.ensureSensor(r, c1(1));
        G  = sr.gamma(r);
        S  = sr.ssim(r);
        fprintf('Gamma %.1f%%  SSIM %.3f  (%s)\n', G.pass_rate, S.ssim_index, c1(1).label());
        ethos.SimPlotter.reconVsTruth(r, G);
        ethos.SimPlotter.sensorPlacement(r);

    case {'comparison', 'sweep'}
        c1 = ethos.SimCase.list(PATIENT_ID, SESSION, cfg, filt{:}, 'CTLabel', 'CT_1');
        c3 = ethos.SimCase.list(PATIENT_ID, SESSION, cfg, filt{:}, 'CTLabel', 'CT_3');
        assert(~isempty(c1) && ~isempty(c3), 'Need both CT_1 and CT_3 fields.');
        rCT1 = sr.resolve(c1(1));   rCT1 = sr.ensureSensor(rCT1, c1(1));
        rCT3 = sr.resolve(c3(1));
        G = sr.gammaChange(rCT1, rCT3);
        fprintf('CT_1 vs CT_3 change gamma: %.1f%%\n', G.pass_rate);
        ethos.SimPlotter.changePanels(rCT1, rCT3, G, cfg.viz_smooth_sigma);
        ethos.SimPlotter.sensorPlacement(rCT1);

        if strcmpi(ANALYSIS_MODE, 'sweep')
            % Pass-rate-vs-criterion for the CT_1-vs-CT_3 recon comparison.
            SW = ethos.Analysis.gammaSweep(rCT1.recon_dose, rCT3.recon_dose, ...
                rCT1.spacing(), CRITERION_N);
            sr.plotPassRateCurve(SW, sprintf('%s: CT_1 vs CT_3', c1(1).label()));
        end

    otherwise
        error('ANALYSIS_MODE must be single | comparison | sweep.');
end
