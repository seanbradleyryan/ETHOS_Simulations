function [ssim_map, mean_over_mask] = compute_local_ssim(reference, target, mask)
%COMPUTE_LOCAL_SSIM Local (per-voxel) structural-similarity map + masked mean.
%
%   [ssim_map, mean_over_mask] = compute_local_ssim(reference, target, mask)
%
%   PURPOSE:
%   Compute the full-volume LOCAL structural similarity (SSIM) map between a
%   target dose and a reference dose, and the mean of that map over an
%   evaluation mask. This is the SSIM analogue of a gamma pass rate: the
%   per-voxel SSIM map is averaged over the same 10%-of-reference region the
%   gamma pass rate is computed on, so the two scores are directly comparable.
%
%   This matches the local-SSIM metric used by study_pass_rates_allsegments
%   (eval_method = 'ssim'); it is NOT the global / per-slice SSIM in
%   step3_analysis (compute_dose_ssim). Both are legitimate but different.
%
%   INPUTS:
%       reference - 3D (or 2D) reference dose array. Sets the SSIM dynamic
%                   range (its max), so pass the volume the mask is built from.
%       target    - Array the same size as reference (the recon or other dose).
%       mask      - Optional logical array the same size as reference. When
%                   given and non-empty, mean_over_mask is the mean SSIM over
%                   its true voxels; otherwise it is the mean over all voxels.
%
%   OUTPUTS:
%       ssim_map       - Full-size local SSIM map (double, 0..1, 1 = identical).
%                        Returns [] when neither volume has any positive signal
%                        or the built-in ssim call fails, so callers can fall
%                        back gracefully.
%       mean_over_mask - Mean of ssim_map over the mask (fraction 0..1), or NaN
%                        when ssim_map is [] or the mask selects nothing.
%
%   ALGORITHM:
%       1. Dynamic range = max(reference); fall back to max(target).
%       2. [~, ssim_map] = ssim(target, reference, 'DynamicRange', dr).
%       3. mean_over_mask = mean(ssim_map(mask)) (or over all voxels).
%
%   EXAMPLE:
%       m               = rs_CT1 >= 0.10 * max(rs_CT1(:));
%       [smap, meanval] = compute_local_ssim(rs_CT1, recon_CT1, m);
%       ssim_pct        = 100 * meanval;   % pass-rate-like percentage
%
%   DEPENDENCIES:
%       - Image Processing Toolbox (built-in ssim with a map output).
%
%   See also: ssim, CalcGamma, step25_segment_metrics, step3_analysis

    if nargin < 3
        mask = [];
    end

    reference = double(reference);
    target    = double(target);

    % Dynamic range from the reference (fall back to the target), as in the
    % study. No positive signal anywhere -> SSIM is undefined; bail out.
    dr = max(reference(:));
    if ~(dr > 0)
        dr = max(target(:));
    end
    if ~(dr > 0)
        ssim_map       = [];
        mean_over_mask = NaN;
        return;
    end

    if exist('ssim', 'file') ~= 2
        error('compute_local_ssim:NoSSIM', ...
            ['Built-in ssim not found (Image Processing Toolbox required). ' ...
             'Cannot compute the local SSIM map.']);
    end

    try
        [~, ssim_map] = ssim(target, reference, 'DynamicRange', dr);
    catch
        ssim_map       = [];
        mean_over_mask = NaN;
        return;
    end

    if ~isempty(mask) && any(mask(:))
        mean_over_mask = mean(ssim_map(mask), 'omitnan');
    else
        mean_over_mask = mean(ssim_map(:), 'omitnan');
    end
end
