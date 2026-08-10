%% =========================================================================
%  RUN_DT_COMPARISON.m
%  Time-step (dt) convergence sweep for a single field, on the +ethos layer.
%  dt is set by the CFL condition (dt = cfl * dx_min / c_max), so coarsening dt
%  by a multiplier is exactly scaling cfl_number. Each multiplier is re-simulated
%  and gamma-compared (3%/3 mm) against the multiplier-1 recon -- a self-
%  convergence metric (ScoreVs 'first'). Air stays at the physically-correct
%  343 m/s (the threshold_2 default).
%
%  Replaces the former ~1500-line self-contained sweep. StudyRunner.paramSweep
%  keeps every recon so any figure can be drawn afterward.
%  Remote/HIPAA execution only.
%  =========================================================================

clear; clc; close all;

WORKING_DIR = '/mnt/weka/home/80030361/ETHOS_Simulations';
PATIENT_ID  = '1194203';
SESSION     = 'Session_1';
BEAM        = 15;           % dose_1194203_Session_1_reference_CT_1_B15_112.mat
SEGMENT     = 112;
BASE_CFL        = 0.3;
DT_MULTIPLIERS  = [1 1.25 1.5 1.75 2 3 4];

cfg = ethos.SimConfig.standalone().withOverrides( ...
    'working_dir', WORKING_DIR, 'patient_id', PATIENT_ID, 'session', SESSION, ...
    'gruneisen_method', 'threshold_2', 'num_time_reversal_iter', 10, ...
    'normalize', true, 'plot_results', true);

sr = ethos.StudyRunner(cfg);
c  = ethos.SimCase.list(PATIENT_ID, SESSION, cfg, 'Beam', BEAM, ...
    'Segment', SEGMENT, 'PlanType', 'reference', 'CTLabel', 'CT_1');

% dt multiplier -> cfl_number = BASE_CFL * multiplier; score every run vs mult=1.
W = sr.paramSweep(c(1), 'cfl_number', DT_MULTIPLIERS, ...
    'Apply',   @(cf, m) cf.withOverrides('cfl_number', BASE_CFL * m), ...
    'ScoreVs', 'first', ...
    'Label',   'dt multiplier (x CFL-limited step)');

% W.runs(k).result / .gamma / .ssim hold every reconstruction and its metrics.
save(fullfile(WORKING_DIR, 'dt_comparison_results.mat'), 'W');
