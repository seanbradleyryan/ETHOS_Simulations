%% =========================================================================
%  RUN_STANDALONE_SIMULATION.m
%  Standalone single-field k-Wave forward + time-reversal, on the +ethos
%  class layer. Replaces the former ~2000-line self-contained script (its
%  inline forward/recon logic now routes through run_single_field_simulation
%  via ethos.Simulator; its CONFIG block, medium/sensor builders and plotting
%  helpers are the classes). Load-or-simulate is automatic (DoseSource).
%
%  Remote/HIPAA execution only -- do not run locally.
%  See CLAUDE-SIMULATION_CONTEXT.md ("+ethos Class Layer").
%  =========================================================================

clear; clc; close all;

% ---- What to run --------------------------------------------------------
WORKING_DIR = '/mnt/weka/home/80030361/ETHOS_Simulations';
PATIENT_ID  = '1194203';
SESSION     = 'Session_1';
BEAM        = 15;          % dose_1194203_Session_1_reference_CT_1_B15_112.mat
SEGMENT     = 112;
CT_LABEL    = 'CT_1';
PLAN_TYPE   = 'reference';

% ---- Config: standalone preset + the knobs this script used -------------
cfg = ethos.SimConfig.standalone().withOverrides( ...
    'working_dir', WORKING_DIR, 'patient_id', PATIENT_ID, 'session', SESSION, ...
    'gruneisen_method', 'threshold_2', ...
    'correction_factor', 20, ...        % ad-hoc scalar the old script ran with
    'conv_noise_level', 0.125, ...
    'reconstruction_method', 'tr', ...
    'num_time_reversal_iter', 10, ...
    'normalize', true, ...
    'dose_source', 'auto');             % load cached recon if present, else sim

fprintf('Config hash: %s\n', cfg.hash());

% ---- Resolve one field (load-or-simulate) -------------------------------
sr      = ethos.StudyRunner(cfg);
oneCase = ethos.SimCase.list(PATIENT_ID, SESSION, cfg, 'Beam', BEAM, ...
    'Segment', SEGMENT, 'PlanType', PLAN_TYPE, 'CTLabel', CT_LABEL);
result  = sr.resolve(oneCase(1));
result  = sr.ensureSensor(result, oneCase(1));

% ---- Metrics (gamma default 3%/3 mm; SSIM over the same mask) -----------
G = sr.gamma(result);
S = sr.ssim(result);
fprintf('Gamma pass rate: %.1f%%  (%g%%/%g mm)\n', ...
    G.pass_rate, G.criteria(1), G.criteria(2));
fprintf('SSIM (masked):   %.3f\n', S.ssim_index);

% ---- Figures ------------------------------------------------------------
ethos.SimPlotter.reconVsTruth(result, G, {'RayStation truth', 'Recon'}, cfg.viz_smooth_sigma);
ethos.SimPlotter.sensorPlacement(result);
if result.isSimulated()
    ethos.SimPlotter.convergence(result);
end

% ---- Optional: cache the recon to disk for later analyze-only runs -------
% sr.sim.save_recon = true; result = sr.sim.run(oneCase(1));
