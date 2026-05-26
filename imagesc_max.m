function imagesc_max(vol, varargin)
% PURPOSE
%   Display three orthogonal imagesc slices through the maximum value of a
%   3-D array (transverse, coronal, sagittal views at the peak voxel).
%
% INPUTS
%   vol      - 3-D numeric array  [ny x nx x nz]
%   varargin - optional name passed to suptitle / figure title (char)
%
% EXAMPLE
%   imagesc_max(doseArray, 'Dose [Gy]')

if nargin < 2
    titleStr = '';
else
    titleStr = varargin{1};
end

[~, idx] = max(vol(:));
[iy, ix, iz] = ind2sub(size(vol), idx);

sliceXY = squeeze(vol(:, :, iz));   % transverse  at max z
sliceXZ = squeeze(vol(iy, :, :));   % coronal     at max y  — [nx x nz]
sliceYZ = squeeze(vol(:, ix, :));   % sagittal    at max x  — [ny x nz]

figure;

subplot(1, 3, 1);
imagesc(sliceXY);
axis image; colormap hot; colorbar;
xlabel('x'); ylabel('y');
title(sprintf('Transverse (z = %d)', iz));

subplot(1, 3, 2);
imagesc(sliceXZ.');
axis image; colormap hot; colorbar;
xlabel('x'); ylabel('z');
title(sprintf('Coronal (y = %d)', iy));

subplot(1, 3, 3);
imagesc(sliceYZ.');
axis image; colormap hot; colorbar;
xlabel('y'); ylabel('z');
title(sprintf('Sagittal (x = %d)', ix));

if ~isempty(titleStr)
    sgtitle(titleStr);
end
end
