%% Setup Parameters

clear; clc; close all; 

patient_ids = {'1194203'}; 
sessions = {'Session_1','Session_2'}; 
imagepath = 'Pancreas/Session_1/';


dicomPath = '/mnt/weka/home/80030361/iraigeneral/EthosExports/1194203/Pancreas/Session_1/sct'; 
matRadPath = '/mnt/weka/home/80030361/MATLAB/Addons/matRad'; 
matRadPath2 = fullfile(matRadPath,'matRad'); 
outputDir = './quick_results/'; 

%% Setup

if ~exist('matRad_cfg','file')
    fprintf('Adding matRad to path ...\n'); 
    addpath(genpath(matRadPath)); 
    addpath(genpath(matRadPath2)); 
end

% Create output dir
if ~exist(outputDir, 'dir')
    mkdir(outputDir); 
end

%% Load dicom files

fprintf('\nStep 1: Load dicom files'); 
files = dir(fullfile(dicomPath, '*.dcm')); 
if isempty(files)
    error('No DICOM files found in: %s',dicomPath); 
end

% Categorize by modality

rtplanFile = ''; 
rtdoseFile = '';
rtstructFile = ''; 
fprintf('Scanning DICOM files...\n'); 
for i = 1:length(files)
    filepath = fullfile(files(i).folder,files(i).name); 
    try
        info = dicominfo(filepath); 
        switch info.Modality
            case 'RTPLAN'
                rtplanFile = filepath; 
                fprintf(' Found RTPLAN: %s\n',files(i).name); 
            case 'RTDOSE'
                rtdoseFile = filepath; 
                fprintf(' Found RTDOSE: %s\n',files(i).name); 
            case 'RTSTRUCT'
                rtstructFile = filepath; 
                fprintf('Found RTSTRUCT: %s\n', files(i).name); 
        end
    catch
    end
end

if isempty(rtplanFile) | isempty(rtdoseFile) | isempty(rtstructFile)
    error('At least one RT file not found')
end


pln = matRad_importDicomRTPlan(rtplanFile)
ct = matRad_importDicomCt(dicomPath); 


