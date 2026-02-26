function N_opt = find_optimal_kwave_size(N,pml_size)
%FIND_OPTIMAL_KWAVE_SIZE Find FFT-optimal grid dimension for k-Wave
%
%   N_opt = find_optimal_kwave_size(N)
%
%   PURPOSE:
%       k-Wave uses FFTs whose performance depends on the prime
%       factorization of grid dimensions.  Dimensions that are products
%       of small primes (2, 3, 5) run significantly faster than those
%       with large prime factors.  This function finds the smallest
%       integer in [N, ceil(1.1*N)] whose greatest prime factor is
%       minimized.
%
%   INPUTS:
%       N - Original grid dimension (positive integer)
%
%   OUTPUTS:
%       N_opt - Optimal padded dimension (N <= N_opt <= ceil(1.1*N))
%               Returns N unchanged if it already has the smallest
%               greatest prime factor in the search range.
%
%   ALGORITHM:
%       1. Search all candidates in [N .. ceil(1.1*N)]
%       2. For each candidate, compute max(factor(candidate))
%       3. Return candidate with smallest max prime factor
%          (ties broken by smallest dimension)
%
%   EXAMPLE:
%       find_optimal_kwave_size(250)  % returns 256 (2^8, max factor = 2)
%       find_optimal_kwave_size(128)  % returns 128 (already 2^7)
%
%   DEPENDENCIES:
%       MATLAB factor()
%
%   See also: run_single_field_simulation, run_standalone_simulation

    N = N + 2*pml_size; 
    if N <= 0 || N ~= round(N)
        error('find_optimal_kwave_size:InvalidInput', ...
            'N must be a positive integer. Got %.4f.', N);
    end

    N_max = ceil(1.1 * N);

    best_max_factor = max(factor(N));
    best_candidate = N;

    for candidate = (N + 1):N_max
        max_pf = max(factor(candidate));
        if max_pf < best_max_factor
            best_max_factor = max_pf;
            best_candidate = candidate;
        end
    end

    N_opt = best_candidate;
    N_opt = N_opt - 2*pml_size; 

end
