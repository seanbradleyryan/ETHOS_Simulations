%% =========================================================================
%  RUN_STANDALONE_COMPARISON.m
%  Two-dose (CT_1 vs CT_3) photoacoustic comparison on the +ethos class layer.
%  The reference-CT dose is reconstructed normally; the other-CT dose is
%  forward-propagated on its own geometry but time-reversed BLINDLY on the
%  reference geometry (ethos.Simulator.runBlind wraps run_second_field_simulation).
%  The two reconstructions are gamma-compared at 3%/3 mm.
%
%  Replaces the former ~2350-line self-contained script (inline dual sims,
%  CONFIG block, and plotting helpers are now the classes).
%  Remote/HIPAA execution only -- do not run locally.
%  =========================================================================

clear; clc; close all;

% ---- What to compare ----------------------------------------------------
WORKING_DIR  = '/mnt/weka/home/80030361/ETHOS_Simulations';
PATIENT_ID   = '1194203';
SESSION      = 'Session_1';
BEAM         = 15;          % dose_1194203_Session_1_reference_CT_3_B15_112.mat
SEGMENT      = 112;
PLAN_TYPE    = 'reference';
REFERENCE_CT = 'CT_1';      % geometry we invert on (live-IRAI reference)
ADAPTED_CT   = 'CT_3';      % true anatomy forward-propagated, then blind-recon

% ---- Config: standalone preset + this script's knobs --------------------
cfg = ethos.SimConfig.standalone().withOverrides( ...
    'working_dir', WORKING_DIR, 'patient_id', PATIENT_ID, 'session', SESSION, ...
    'gruneisen_method', 'threshold_2', ...
    'correction_factor', 20, ...
    'conv_noise_level', 0.125, ...
    'num_time_reversal_iter', 10, ...
    'normalize', true, ...
    'blind_recon_ct3', true, ...
    'blind_recon_reference_label', REFERENCE_CT, ...
    'dose_source', 'simulate', ...
    'plot_results', true);

% ---- Locate both CT counterparts of this beam/segment -------------------
sr    = ethos.StudyRunner(cfg);
refC  = ethos.SimCase.list(PATIENT_ID, SESSION, cfg, 'Beam', BEAM, ...
    'Segment', SEGMENT, 'PlanType', PLAN_TYPE, 'CTLabel', REFERENCE_CT);
fwdC  = ethos.SimCase.list(PATIENT_ID, SESSION, cfg, 'Beam', BEAM, ...
    'Segment', SEGMENT, 'PlanType', PLAN_TYPE, 'CTLabel', ADAPTED_CT);

% ---- Reference full recon vs blind recon of the adapted anatomy ---------
C = sr.compareBlind(fwdC(1), refC(1));

fprintf('Change-detection gamma pass: %.1f%%  (%g%%/%g mm)\n', ...
    C.gamma_change.pass_rate, C.gamma_change.criteria(1), C.gamma_change.criteria(2));

% C.reference_result / C.blind_result are SimResults; draw any extra figures
% from them, e.g.:
%   ethos.SimPlotter.reconVsTruth(C.reference_result, sr.gamma(C.reference_result));
