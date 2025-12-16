%% Load CT stuff
clear all

id = {'1194203'}; 
session = {'Session_1'}; 

ctFolder = fullfile(pwd,'EthosExports',id,'Pancreas',session); 
dicoms = dicomCollection(ctFolder{1}); 
%% Dosegrid setup
rowIndex = strcmp(dicoms.Modality,"RTDOSE"); 
filePath = dicoms.Filenames{rowIndex};
info = dicominfo(filePath);
rawDose = dicomread(filePath);
doseGrid = double(squeeze(rawDose)); 
scaling = info.DoseGridScaling; 
doseGrid = doseGrid * scaling; 
disp(['Dose Units: ' info.DoseUnits]);
disp(['Voxel Spacing (x,y,z): ' num2str(info.PixelSpacing(2)) ' mm (x), ' ...
      num2str(info.PixelSpacing(1)) ' mm (y), ' ...
      num2str(abs(info.GridFrameOffsetVector(2) - info.GridFrameOffsetVector(1))) ' mm (z)']);

n_pulses = 1200; % @ Change me
doseGrid = doseGrid / n_pulses; 


%% Prepare Density grid

ctRowIndex = strcmp(dicoms.SeriesDescription, 'sct');
ctRow = find(ctRowIndex, 1);
ctFilenames = dicoms.Filenames{ctRow};  % Cell array of filenames for CT series
[ctVolume, ctSpatial] = dicomreadVolume(ctFilenames);  % Read CT volume
ctHU = squeeze(ctVolume);  % 3D array of HU (initial raw values)

% Get metadata from the first CT slice for coordinate setup
firstFile = ctFilenames{1};
firstInfo = dicominfo(firstFile);
ipp = firstInfo.ImagePositionPatient;  % Top-left corner [x, y, z]
pixelSpacing = firstInfo.PixelSpacing;  % [row_spacing (dy), col_spacing (dx)] in mm
sliceThickness = firstInfo.SliceThickness;
pixelSpacing = [pixelSpacing; sliceThickness];  % [dy; dx; dz]
numRows = firstInfo.Rows;  % ny (rows, y-dir)
numCols = firstInfo.Columns;  % nx (cols, x-dir)
numSlices = length(ctSpatial.PatientPositions(:,3));
% Apply rescale to get actual HU values
slope = double(firstInfo.RescaleSlope);
intercept = double(firstInfo.RescaleIntercept);
ctHU = double(ctHU);
ctHU = slope * ctHU + intercept;
% Rest of the script remains unchanged...
% Get ImageOrientationPatient to determine direction cosines
iop = firstInfo.ImageOrientationPatient;
col_cos = iop(1:3);  % Direction for increasing column (x-dir typically)
row_cos = iop(4:6);  % Direction for increasing row (y-dir typically)

% Compute signed spacings (assuming axis-aligned for simplicity)
dx = pixelSpacing(2) * col_cos(1);  % Signed dx for x-dir
dy = pixelSpacing(1) * row_cos(2);  % Signed dy for y-dir
% For z, use actual positions (assuming positive, but we'll check)

% Generate coordinate vectors with signed increments
ctX = double(ipp(1)) + double((0:numCols-1)) * double(dx);  % X coordinates (columns)
ctY = double(ipp(2)) + double((0:numRows-1)) * double(dy);  % Y coordinates (rows)
ctZ = ctSpatial.PatientPositions(:,3)';  % Z coordinates (slices); use actual positions, assumed increasing
ctZ = double(ctZ);

% Ensure coordinates are monotonically increasing; flip array dimensions if necessary
if ctX(2) < ctX(1)  % If decreasing
    ctX = flip(ctX);
    ctHU = flip(ctHU, 2);  % Flip along columns (dim2)
end
if ctY(2) < ctY(1)  % If decreasing
    ctY = flip(ctY);
    ctHU = flip(ctHU, 1);  % Flip along rows (dim1)
end
if ctZ(2) < ctZ(1)  % If decreasing (unlikely, but check)
    ctZ = flip(ctZ);
    ctHU = flip(ctHU, 3);  % Flip along slices (dim3)
end

% Assume dose grid metadata (from RTDOSE dicominfo)
doseInfo = info;
doseOrigin = doseInfo.ImagePositionPatient;  % [x y z]
doseSpacing = doseInfo.PixelSpacing;  % [row_spacing (dy), col_spacing (dx)]
doseSliceOffsets = doseInfo.GridFrameOffsetVector;  % Z offsets
doseDims = [doseInfo.Rows, doseInfo.Columns, numel(doseSliceOffsets)];  % [ny (rows), nx (cols), nz (slices)]

% Use same direction cosines for dose grid (since referenced to CT)
dose_dx = doseSpacing(2) * col_cos(1);  % Signed dx
dose_dy = doseSpacing(1) * row_cos(2);  % Signed dy

% Dose Z coordinates (assume increasing, but sort if needed)
doseZ = doseOrigin(3) + doseSliceOffsets;  % Z coords
if doseZ(2) < doseZ(1)
    doseZ = flip(doseZ);
    doseGrid = flip(doseGrid, 3);  % Flip dose values if needed (though not necessary for query)
end
doseZ = double(doseZ);

% Create vectors for dose coordinates
doseX_vec = double(doseOrigin(1)) + double((0:doseDims(2)-1)) * double(dose_dx);
doseY_vec = double(doseOrigin(2)) + double((0:doseDims(1)-1)) * double(dose_dy);
doseZ_vec = doseZ;

% Create query grids for dose voxels using ndgrid (order: y, x, z for consistency with interp3? No, ndgrid outputs vary accordingly)
[doseGridY, doseGridX, doseGridZ] = ndgrid(doseY_vec, doseX_vec, doseZ_vec);  % Note: swapped to match interp3 order (Y varies dim1, X dim2, Z dim3)
doseGridX = double(doseGridX); 
doseGridY = double(doseGridY); 
doseGridZ = double(doseGridZ); 

% Debugging: Print coordinate ranges to check for overlap
disp('CT Coordinate Ranges:');
disp(['X: ' num2str(min(ctX)) ' to ' num2str(max(ctX))]);
disp(['Y: ' num2str(min(ctY)) ' to ' num2str(max(ctY))]);
disp(['Z: ' num2str(min(ctZ)) ' to ' num2str(max(ctZ))]);

disp('Dose Grid Query Ranges:');
disp(['X: ' num2str(min(doseGridX(:))) ' to ' num2str(max(doseGridX(:)))]);
disp(['Y: ' num2str(min(doseGridY(:))) ' to ' num2str(max(doseGridY(:)))]);
disp(['Z: ' num2str(min(doseGridZ(:))) ' to ' num2str(max(doseGridZ(:)))]);

% If ranges don't overlap, that's why extrapolation occurs. Adjust origins or check DICOM metadata.

% Resample CT HU to dose grid (scaling up/down as needed)
ctHU = double(ctHU); 
resampledHU = interp3(ctX,ctY,ctZ,ctHU,doseGridX,doseGridY,doseGridZ,'linear',-1000); 

% density = 1 + resampledHU; 
% density(density < .001) = .001; 
% density = density * 1000; 

save('ctdensity.mat', 'resampledHU');
