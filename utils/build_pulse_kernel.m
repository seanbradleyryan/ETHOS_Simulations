function kernel = build_pulse_kernel(pulse_shape, pulse_width_s, dt)
%BUILD_PULSE_KERNEL  Normalized temporal radiation-pulse kernel for the sensor model.
%
%   kernel = build_pulse_kernel(pulse_shape, pulse_width_s, dt)
%
%   Returns the discrete, unit-sum pulse profile that each sensor time series is
%   convolved with (and later Wiener-deconvolved) in the acquisition model. The
%   shape represents the finite radiation-pulse duration of the linac:
%
%     'gaussian'    - Gaussian profile, sigma = pulse_width_s, truncated at 4
%                     sigma. This is the original/default behavior.
%     'rectangular' - Flat-top (boxcar) profile of TOTAL width pulse_width_s.
%                     More consistent with a clinical linac's ~few-us beam-on
%                     burst, which is closer to a top-hat than a Gaussian.
%
%   Note the sign of the width argument differs by shape: for 'gaussian' it is
%   the Gaussian sigma; for 'rectangular' it is the full pulse width. A 4 us
%   rectangle is therefore temporally narrower than a 4 us-sigma Gaussian
%   (FWHM ~= 9.4 us), so the raw (pre-deconvolution) signal is less smeared.
%
%   The kernel is returned as a column vector normalized to sum to 1 so it
%   preserves the DC level of the signal. It is not zero-lag centered; the
%   forward convolution and the matched deconvolution use the same kernel, so
%   the resulting time shift cancels exactly.
%
%   INPUTS:
%       pulse_shape   - 'gaussian' | 'rectangular' (case-insensitive)
%       pulse_width_s - Gaussian sigma (s) or rectangular full width (s)
%       dt            - sensor sampling interval (s)
%
%   OUTPUT:
%       kernel        - column vector, unit sum
%
%   See also: run_single_field_simulation, run_standalone_simulation_upgraded

    switch lower(char(pulse_shape))
        case 'gaussian'
            sigma_samples = pulse_width_s / dt;
            kernel_half   = ceil(4 * sigma_samples);
            t_kernel      = (-kernel_half : kernel_half)';
            kernel        = exp(-t_kernel.^2 / (2 * sigma_samples^2));

        case 'rectangular'
            width_samples = max(1, round(pulse_width_s / dt));
            kernel        = ones(width_samples, 1);

        otherwise
            error('build_pulse_kernel:BadShape', ...
                'pulse_shape must be ''gaussian'' or ''rectangular'' (got ''%s'').', ...
                char(pulse_shape));
    end

    kernel = kernel / sum(kernel);   % unit-sum normalization
end
