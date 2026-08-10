%% =========================================================================
%  STUDY_PASS_RATES_ALLSEGMENTS.m
%  Per-BEAM photoacoustic gamma/SSIM summary, averaged over ALL segments, on the
%  +ethos class layer (load-only, no k-Wave). For every CT_1 field with a CT_3
%  twin, StudyRunner.beamSummary scores four per-segment comparisons and reduces
%  to a per-beam mean +/- std, then draws the change-detection and reconstruction-
%  fidelity figures.
%
%    truth_CT1 vs truth_CT3 | truth_CT1 vs recon_CT1
%    truth_CT1 vs recon_CT3 | truth_CT3 vs recon_CT3
%
%  Replaces the former ~1280-line script. Metric is gamma (3%/3 mm) or SSIM over
%  the same 10%-of-truth mask. Missing/mismatched segments are skipped.
%  Remote/HIPAA execution only.
%  =========================================================================

clear; clc; close all;

WORKING_DIR = '/mnt/weka/home/80030361/ETHOS_Simulations';
PATIENT_ID  = '1194203';
SESSION     = 'Session_1';
METRIC      = 'gamma';       % 'gamma' | 'ssim'
PLAN_TYPE   = 'reference';

cfg = ethos.SimConfig.standalone().withOverrides( ...
    'working_dir', WORKING_DIR, 'patient_id', PATIENT_ID, 'session', SESSION, ...
    'gruneisen_method', 'threshold_2', ...
    'dose_source', 'load', ...       % analysis-only
    'gamma_default', [3 3], ...
    'plot_results', true);

sr = ethos.StudyRunner(cfg);
BS = sr.beamSummary('Metric', METRIC, 'PlanType', PLAN_TYPE);

out_dir = fullfile(WORKING_DIR, 'AnalysisResults', PATIENT_ID, SESSION);
if ~isfolder(out_dir), mkdir(out_dir); end
save(fullfile(out_dir, sprintf('pass_rates_allsegments_%s.mat', METRIC)), 'BS');
