close all
clear all
EthosPrepCT; 

densityfile = load("ctdensity.mat"); 
HU = densityfile.resampledHU; 

density = zeros(size(HU)); 
soundspeed = zeros(size(HU)); 
gruneisen = zeros(size(HU)); 
attenuation = zeros(size(HU)); 

%% Defines dx,dy,dz based on the info in the dicom file before padding by getting d from the largest dimension array
[orig_ny, orig_nx, orig_nz] = size(doseGrid);
spacing_y_mm = info.PixelSpacing(1);  % rows (y)
spacing_x_mm = info.PixelSpacing(2);  % columns (x)
spacing_z_mm = abs(info.GridFrameOffsetVector(2) - info.GridFrameOffsetVector(1));  % slices (z)
dim_lengths = [orig_ny, orig_nx, orig_nz];
spacings_mm = [spacing_y_mm, spacing_x_mm, spacing_z_mm];
% Find the dimension with the largest array length
[~, max_dim_idx] = max(dim_lengths);
chosen_spacing_mm = spacings_mm(max_dim_idx);
dx = spacings_mm(1);
dy = spacings_mm(2);
dz = spacings_mm(3); % Scale by the inverse of the amount streteched 

pml_size = 10; 

dosegriddims = size(doseGrid); 
Nx = dosegriddims(1); 
Ny = dosegriddims(2); 
Nz = dosegriddims(3); 
kgrid = kWaveGrid(Nx,dx,Ny,dy,Nz,dz); 

% kgrid.Nt = 600; 
% kgrid.dt = 1e-7; 
% 
%% Setup acoustic parameters

density(HU < -200 ) = 1000; % Water/Air 
density(HU >= -200 & HU < -50) = 1000; % Fat
density(HU >= -50 & HU < 100) = 1040; % Soft Tissue
density(HU >= 100) = 1900; % Bone

soundspeed(HU < -200 ) = 1500; % Water/Air 
soundspeed(HU >= -200 & HU < -50) = 1480; % Fat
soundspeed(HU >= -50 & HU < 100) = 1540; % Soft Tissue
soundspeed(HU >= 100) = 2000; % Bone

% attenuation(HU < -200 ) = .0022; % Water/Air 
% attenuation(HU >= -200 & HU < -50) = .5; % Fat
% attenuation(HU >= -50 & HU < 100) = 1; % Soft Tissue
% attenuation(HU >= 100) = 10; % Bone

gruneisen(HU < -200 ) = .11; % Water/Air 
gruneisen(HU >= -200 & HU < -50) = .8; % Fat
gruneisen(HU >= -50 & HU < 100) = .3; % Soft Tissue
gruneisen(HU >= 100) = .8; % Bone

medium.density = density; 
medium.sound_speed = soundspeed; 
medium.alpha_coeff = attenuation; 
medium.alpha_power = 1.85; 

%% Pressure and sensor setup

p0 = doseGrid .* gruneisen .* medium.density / 100; 
source.p0 = p0; 

sensor.mask = zeros(Nx,Ny,Nz); 
sensor.mask(10,50:100,60:110) = 1; 

% imagesc(squeeze(density(12,60:110,50:100))); 
% imagesc(squeeze(density(12,:,:)))
% imagesc(squeeze(doseGrid(maxx,:,:))); 
% imagesc(squeeze(doseGrid(maxx,50:100,60:110)))
% Place Sensor

input_args = {'Smooth', false, 'PMLInside', false, 'PMLSize', pml_size, 'DataCast', 'gpuArray-single', 'PlotSim', true};
    sensor_data = kspaceFirstOrder3D(kgrid,medium,source,sensor,input_args{:});

%% Time Reversal

source = rmfield(source,'p0'); 
source.p_mask = sensor.mask; 
source.p = fliplr(sensor_data); 
source.p_mode = 'dirichlet'; 
sensor.record = {'p_final'}; 
p0_recon = kspaceFirstOrder3D(kgrid,medium,source,sensor,input_args{:}); 
p0_recon.p_final = p0_recon.p_final .* (p0_recon.p_final > 0); 


n_iterations = 5; 

for loop = 2:n_iterations
                disp("Iteration stage 1 for loop # ")
            loop
            p0_recon = gather(p0_recon); 
            source = rmfield(source,'p');  
            source.p0 = p0_recon.p_final; 
            sensor = rmfield(sensor,'record'); 
            sensor_data2 = kspaceFirstOrder3D(kgrid,medium,source,sensor,input_args{:}); 
            data_difference = sensor_data - sensor_data2; 
            disp("Iteration stage 2")
            source.p_mask = sensor.mask; 
            source.p = fliplr(data_difference); 
            source = rmfield(source,'p0'); 
            source.p_mode = 'dirichlet'; 
            sensor.record = {'p_final'}; 
            p0_update = kspaceFirstOrder3D(kgrid,medium,source,sensor,input_args{:}); 
            p0_recon.p_final = p0_recon.p_final + p0_update.p_final; 
            p0_recon.p_final = p0_recon.p_final .* (p0_recon.p_final > 0); 
            % p0_recon.p_final = p0_recon.p_final * 2; 
            eval(['p0_' num2str(loop) ' = gather(p0_recon.p_final)']); 
end

recon_dose = p0_5 ./ (gruneisen .* density); 
recon_dose = gather(recon_dose); 


%% Plotting

[maxval,idx] = max(recon_dose(:)); 

[maxx,maxy,maxz] = ind2sub(size(recon_dose),idx); 

% recon_dose = flip(recon_dose,1); 
% recon_dose = flip(recon_dose,2); 
figure; 
subplot(2,3,1); 
imagesc([0 Nx*dx],[0 Nz*dz],squeeze(doseGrid(:,maxy,:))); 
title("Reference (Gy)")
subplot(2,3,2); 
imagesc([0 Ny*dy],[0 Nz*dz],squeeze(doseGrid(maxx,:,:))); 
subplot(2,3,3); 
imagesc([0 Nx*dx],[0 Ny*dy],squeeze(doseGrid(:,:,maxz))); 
colorbar;
subplot(2,3,4); 
imagesc([0 Nx*dx],[0 Nz*dz],squeeze(recon_dose(:,maxy,:))); 
title("Reconstruction (Gy)")
subplot(2,3,5); 
imagesc([0 Ny*dy],[0 Nz*dz],squeeze(recon_dose(maxx,:,:))); 
subplot(2,3,6); 
imagesc([0 Nx*dx],[0 Ny*dy],squeeze(recon_dose(:,:,maxz))); 
colorbar; 
