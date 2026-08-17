%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Parameters %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pitch = 3.65; % in mm
spacing = 0.2; % in mm
nElementsPerSide = 32;
speedOfSound = 1540; % Speed of sound in solidified oil in meters per second (m/s
%speedOfSound = 1430; % Speed of sound in solidified oil in meters per second (m/s)
Dr = 200; % Image size in pixels
D = 800; % Number of pixels in depth (x)

% Define the lateral, elevation, and depth ranges (in mm)
lateralRange_mm = [-50, 50]; % Lateral range in mm (-5 cm to +5 cm)
elevationRange_mm = [-50, 50]; % Elevation range in mm (-5 cm to +5 cm)
depthRange_mm = [0, 200]; % Depth range in mm (0 cm to 20 cm)

% Flag to choose between empirical and calculated element spacing
useEmpiricalSpacing = true;

% Convert pitch and spacing from mm to meters for consistency
pitch_meters = pitch / 1000; % Convert from mm to meters
spacing_meters = spacing / 1000; % Convert from mm to meters
lateralRange = lateralRange_mm / 1000; % Convert lateral range to meters
elevationRange = elevationRange_mm / 1000; % Convert elevation range to meters
depthRange = depthRange_mm / 1000; % Convert depth range to meters

% Calculate element spacing based on the flag
if useEmpiricalSpacing
    % Empirical element spacing formula (used in the working code)
    elmspace = 0.85 * 1.54 / 0.35; % space between elements (empirical) in meters
    disp('Using empirical element spacing:');
else
    % Element spacing based on the provided pitch and spacing (converted to meters)
    elmspace = pitch_meters + spacing_meters; % calculated from ultrasound parameters in meters
    disp('Using calculated element spacing:');
end
disp(['Element Spacing = ', num2str(elmspace), ' meters']);

% Updated Sampling frequency
fs = 0.7812 * 4; % in MHz, still needs to be converted to Hz for time delay calculations
fs_Hz = fs * 1e6; % Convert from MHz to Hz

% Initialize element positions matrix (1024 x 2) in meters
elementPositions = zeros(nElementsPerSide^2, 2);

% Define element positions in a 32x32 grid (converted to meters)
for i = 1:nElementsPerSide
    for j = 1:nElementsPerSide
        % Calculate index
        index = (i-1) * nElementsPerSide + j;
        
        % Assign positions (centered around origin (0,0)) in meters
        elementPositions(index, 1) = (i - (nElementsPerSide + 1)/2) * elmspace;
        elementPositions(index, 2) = (j - (nElementsPerSide + 1)/2) * elmspace;
    end
end

% Display the element positions for verification (optional)
disp('Element Positions (in meters):');
disp(elementPositions);

% Initialize the grid for lateral (y), elevation (z), and depth (x) ranges
lateralRangeVals = linspace(lateralRange(1), lateralRange(2), Dr); % Lateral (y-axis) values in meters
elevationRangeVals = linspace(elevationRange(1), elevationRange(2), Dr); % Elevation (z-axis) values in meters
depthRangeVals = linspace(depthRange(1), depthRange(2), D); % Depth (x-axis) values in meters

% Create meshgrid for lateral and elevation focus points
[yGrid, zGrid] = meshgrid(lateralRangeVals, elevationRangeVals); % Grid for y and z

% Initialize the padded image matrix
im_pad = zeros(2 * Dr + 1, 2 * Dr + 1, D);
im_pad_dump = im_pad; % optional dump for debugging

% Set up a parallel pool with 64 cores
parpool('local', 64);

% Parallelized computation over array elements instead of focus points
parfor iii = 3:30
    disp(['Processing element: ', num2str(iii)]);
    im_pad = im_pad + pcalculation(iii, Dr, elmspace, Data1024, speedOfSound, fs_Hz, lateralRangeVals, elevationRangeVals, depthRangeVals);
end

% Perform envelope detection and scaling on the image
parfor i = 1:(2 * Dr + 1)
    C1(:, i, :) = envelopeDetection(squeeze(im_pad(:, i, :)));
end

% Close the parallel pool after computation
delete(gcp('nocreate'));

% Display and save the final scaled image
figure; imagesc(lateralRange_mm, elevationRange_mm, squeeze(C1(:, :, D/2))); % Scale to the middle depth slice
axis image;
save('output_file.mat', 'C1', '-v7.3');
disp("DAS finished");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% pcalculation Function %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function im_pad_temp = pcalculation(iii, Dr, elmspace, PAD, speedOfSound, fs_Hz, lateralRangeVals, elevationRangeVals, depthRangeVals)
    im_pad_temp = zeros(2 * Dr + 1, 2 * Dr + 1, length(depthRangeVals));
    for ii = 3:30
        % Calculate element positions (converted to meters)
        xi = Dr + 1 + (-16.5 + iii) * elmspace;
        yi = Dr + 1 + (-16.5 + ii) * elmspace;
        zi = 1;
        
        % Iterate over the entire 3D grid
        for x = 1:2 * Dr + 1
            for y = 1:2 * Dr + 1
                for z = 1:length(depthRangeVals)
                    r = floor(sqrt((x - xi)^2 + (y - yi)^2 + (z - zi)^2)); % Calculate distance in meters
                    ppi = r - 0;
                    
                    % Ensure the delay is within bounds and apply scaling (z/r.^2)
                    if ppi <= size(PAD, 1) - 1 && ppi > 20 && z / r > 0.997
                        im_pad_temp(x, y, z) = im_pad_temp(x, y, z) + PAD(ppi, iii, ii) * z / r.^2;
                    end
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Envelope Detection Function %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function envData = envelopeDetection(data)
    % Simple envelope detection using the Hilbert transform
    envData = abs(hilbert(data));
end
