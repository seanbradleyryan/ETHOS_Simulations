%% =========================================================================
%  RUN_AIR_SOUND_SPEED_COMPARISON.m
%  Compare reconstruction with in-body air sound speed = 343 m/s (physically
%  correct, slow) vs 1480 m/s (water speed, fast), otherwise identical, and
%  gamma-analyze the difference. On the +ethos layer the air speed is swept via
%  SimConfig.withTissue; each recon is gamma-compared (3%/3 mm) against the
%  343 m/s run (ScoreVs 'first'). Noise is disabled so any difference is
%  attributable to the air-speed change alone.
%
%  Replaces the former ~1400-line self-contained script.
%  Remote/HIPAA execution only.
%  =========================================================================

clear; clc; close all;

WORKING_DIR = '/mnt/weka/home/80030361/ETHOS_Simulations';
PATIENT_ID  = '1194203';
SESSION     = 'Session_1';
BEAM        = 15;           % dose_1194203_Session_1_reference_CT_1_B15_112.mat
SEGMENT     = 112;
AIR_SPEEDS  = [343, 1480];  % first entry is the physically-correct baseline

cfg = ethos.SimConfig.standalone().withOverrides( ...
    'working_dir', WORKING_DIR, 'patient_id', PATIENT_ID, 'session', SESSION, ...
    'gruneisen_method', 'threshold_2', 'num_time_reversal_iter', 10, ...
    'convolution_kernel', 0, ...     % disable the pulse/noise block for both runs
    'normalize', true, 'plot_results', true);

sr = ethos.StudyRunner(cfg);
c  = ethos.SimCase.list(PATIENT_ID, SESSION, cfg, 'Beam', BEAM, ...
    'Segment', SEGMENT, 'PlanType', 'reference', 'CTLabel', 'CT_1');

W = sr.paramSweep(c(1), 'air_sound_speed', AIR_SPEEDS, ...
    'Apply',   @(cf, v) cf.withTissue('threshold_2', 'air', 'sound_speed', v), ...
    'ScoreVs', 'first', ...
    'Label',   'in-body air sound speed (m/s)');

% Difference panel between the two air-speed recons.
if numel(W.runs) >= 2
    ethos.SimPlotter.changePanels(W.runs(1).result, W.runs(2).result, ...
        W.runs(2).gamma, cfg.viz_smooth_sigma);
end
save(fullfile(WORKING_DIR, 'air_sound_speed_results.mat'), 'W');
