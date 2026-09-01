function gamma_results = compute_gamma(original, recon, spacing_mm, varargin)
%COMPUTE_GAMMA Gamma-index analysis of a reconstruction vs a truth dose.
%
%   gamma_results = compute_gamma(original, recon, spacing_mm)
%   gamma_results = compute_gamma(original, recon, spacing_mm, Name, Value, ...)
%
%   Runs CalcGamma for each (pct, dta) criterion over the truth's low-dose
%   region and returns a struct the standalone / comparison drivers save and
%   plot. Shared copy of the gamma block previously duplicated in every driver.
%
%   Name/value options (defaults):
%       'Criteria'  { {10,10,'10%/10mm'}; {5,5,'5%/5mm'}; {3,3,'3%/3mm'} }
%                   Cell array, one row per criterion: {pct, dta_mm, label}.
%       'Cutoff'    0.10   - low-dose cutoff as a fraction of max(original).
%       'Local'     0      - CalcGamma 'local' flag (0 = global).
%       'Restrict'  1      - CalcGamma 'restrict' flag.
%
%   Returns gamma_results with fields .maps (cell), .pass_rates, .criteria,
%   .cutoff_Gy, .eval_mask. Returns [] if CalcGamma is not on the path.
%
%   See also: CalcGamma, run_standalone_field, run_standalone_simulation

    p = struct('Criteria', {{3,3,'3%/3mm'}}, ...
               'Cutoff', 0.10, 'Local', 0, 'Restrict', 1);
    for i = 1:2:numel(varargin)
        p.(varargin{i}) = varargin{i+1};
    end

    if exist('CalcGamma', 'file') ~= 2
        warning('compute_gamma:NoCalcGamma', 'CalcGamma not found; skipping gamma analysis.');
        gamma_results = [];
        return;
    end

    ref_struct = struct('start', [0,0,0], 'width', spacing_mm, 'data', double(original));
    tgt_struct = struct('start', [0,0,0], 'width', spacing_mm, 'data', double(recon));

    low_dose_cutoff = p.Cutoff * max(double(original(:)));
    eval_mask       = double(original) >= low_dose_cutoff;

    crit  = p.Criteria;
    nCrit = size(crit, 1);
    maps  = cell(nCrit, 1);
    pass_rates = zeros(nCrit, 1);

    for gc = 1:nCrit
        pct = crit{gc, 1};
        dta = crit{gc, 2};
        lbl = crit{gc, 3};
        fprintf('  [gamma %s] computing...', lbl);
        try
            gmap = CalcGamma(ref_struct, tgt_struct, pct, dta, ...
                'local', p.Local, 'limit', dta * 2, 'restrict', p.Restrict);
            maps{gc}       = gmap;
            pass_rates(gc) = 100 * mean(gmap(eval_mask) <= 1);
            fprintf(' pass %.2f%%\n', pass_rates(gc));
        catch ME
            warning('compute_gamma:CritFail', 'Gamma [%s] failed: %s', lbl, ME.message);
            maps{gc}       = [];
            pass_rates(gc) = NaN;
            fprintf(' FAILED\n');
        end
    end

    gamma_results = struct();
    gamma_results.maps       = maps;
    gamma_results.pass_rates = pass_rates;
    gamma_results.criteria   = crit;
    gamma_results.cutoff_Gy  = low_dose_cutoff;
    gamma_results.eval_mask  = eval_mask;
end
