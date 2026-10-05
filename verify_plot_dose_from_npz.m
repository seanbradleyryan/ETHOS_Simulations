%VERIFY_PLOT_DOSE_FROM_NPZ Quick visual check of RayStation NPZ field doses.
%
%   Picks 5 random plan/beam/segment combinations that exist for BOTH CT_1
%   and CT_3 in RayStationFiles/<patient_id>/<session>/ and plots each one:
%   one figure per combination, rows = CT_1 / CT_3, columns = transverse /
%   coronal / sagittal slices through that CT's max-dose voxel. Both rows
%   share one color scale so the two CTs can be compared directly.
%
%   NPZ files come from calc_beam_plan_doses.py:
%       dose_{id}_{session}_{plan_type}_{CT_1|CT_3}_{origbeam}_{seg:02d}.npz
%   containing dose.npy as float32, C-order, shape (nz, ny, nx).
%
%   See also: step14_npz_to_mat

%% ---- Settings (edit these) ----
patient_id  = '1194203';
session     = 'Session_1';
working_dir = 'C:/Users/80030361/ETHOS_Simulations';
n_pick      = 5;

%% ---- Find combinations present on both CTs ----
rs_dir = fullfile(working_dir, 'RayStationFiles', patient_id, session);
files  = dir(fullfile(rs_dir, 'dose_*.npz'));

keys_ct1 = {};
keys_ct3 = {};
for i = 1:numel(files)
    tok = regexp(files(i).name, '_(adapted|reference)_(CT_\d)_(B\d+)_(\d+)\.npz$', 'tokens', 'once');
    if isempty(tok)
        continue;
    end
    key = sprintf('%s_%s_%s', tok{1}, tok{3}, tok{4});   % e.g. adapted_B13_02
    if strcmp(tok{2}, 'CT_1')
        keys_ct1{end+1} = key; %#ok<SAGROW>
    elseif strcmp(tok{2}, 'CT_3')
        keys_ct3{end+1} = key; %#ok<SAGROW>
    end
end

common_keys = intersect(keys_ct1, keys_ct3);
fprintf('Found %d NPZ files, %d combinations on both CT_1 and CT_3.\n', ...
    numel(files), numel(common_keys));
if isempty(common_keys)
    error('verify_plot_dose_from_npz:NoPairs', 'No CT_1/CT_3 pairs found in %s', rs_dir);
end

picked = common_keys(randperm(numel(common_keys), min(n_pick, numel(common_keys))));

%% ---- Plot each picked combination ----
ct_labels = {'CT_1', 'CT_3'};
for p = 1:numel(picked)
    % key is "<plan_type>_<beam>_<seg>"; the CT label sits between plan_type and beam
    parts = strsplit(picked{p}, '_');
    figure('Name', picked{p}, 'Color', 'w');

    dose = cell(1, 2);
    for c = 1:2
        npz_name = sprintf('dose_%s_%s_%s_%s_%s_%s.npz', patient_id, session, ...
            parts{1}, ct_labels{c}, parts{2}, parts{3});
        npz_path = fullfile(rs_dir, npz_name);

        % NPZ is a zip of .npy files; pull out dose.npy and read it
        tmp_dir = tempname;
        unzip(npz_path, tmp_dir);
        fid = fopen(fullfile(tmp_dir, 'dose.npy'), 'rb');
        fread(fid, 8, 'uint8');                         % magic (6) + version (2)
        header_len = fread(fid, 1, 'uint16');           % NPY v1.0 header length
        header = char(fread(fid, header_len, 'uint8')');
        shape_str = regexp(header, '''shape''\s*:\s*\(([^)]*)\)', 'tokens', 'once');
        shape  = str2double(regexp(shape_str{1}, '\d+', 'match'));   % [nz ny nx]
        data   = fread(fid, prod(shape), 'single=>double');
        fclose(fid);
        rmdir(tmp_dir, 's');

        % C-order (nz, ny, nx) -> MATLAB (row=Y, col=X, slice=Z)
        dose{c} = permute(reshape(data, [shape(3) shape(2) shape(1)]), [2 1 3]);
        fprintf('[%d/%d] %s  size=[%d %d %d]  max=%.4g Gy\n', p, numel(picked), ...
            npz_name, size(dose{c}, 1), size(dose{c}, 2), size(dose{c}, 3), max(dose{c}(:)));
    end

    clim_max = max([dose{1}(:); dose{2}(:)]);
    if clim_max <= 0
        clim_max = 1;   % empty dose: avoid an invalid color range
    end

    for c = 1:2
        [~, idx] = max(dose{c}(:));
        [iy, ix, iz] = ind2sub(size(dose{c}), idx);

        subplot(2, 3, 3*(c-1) + 1);
        imagesc(dose{c}(:, :, iz)); axis image; caxis([0 clim_max]);
        title(sprintf('%s transverse z=%d', ct_labels{c}, iz), 'Interpreter', 'none');

        subplot(2, 3, 3*(c-1) + 2);
        imagesc(flipud(squeeze(dose{c}(iy, :, :))')); axis image; caxis([0 clim_max]);
        title(sprintf('%s coronal y=%d', ct_labels{c}, iy), 'Interpreter', 'none');

        subplot(2, 3, 3*(c-1) + 3);
        imagesc(flipud(squeeze(dose{c}(:, ix, :))')); axis image; caxis([0 clim_max]);
        title(sprintf('%s sagittal x=%d', ct_labels{c}, ix), 'Interpreter', 'none');
        colorbar;
    end
    sgtitle(sprintf('%s %s  %s', patient_id, session, picked{p}), 'Interpreter', 'none');
end
