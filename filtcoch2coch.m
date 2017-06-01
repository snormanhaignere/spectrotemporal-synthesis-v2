function coch = filtcoch2coch(filtcoch, P, complex_filters)

if nargin < 3
    complex_filters = false;
end

% compute the subbands
accum_FT_transfer_function = zeros(size(filtcoch,1), size(filtcoch,2));
accum_FT_subbands = zeros(size(filtcoch,1), size(filtcoch,2));
for i = 1:length(P.spec_mod_rates)
    for j = 1:length(P.temp_mod_rates)
        
        % filter transfer functions
        Hts = filt_spectemp_mod(P.spec_mod_rates(i), P.temp_mod_rates(j), ...
            size(filtcoch,2), size(filtcoch,1), P, 0, 0, 0, 0, complex_filters);
        
        % accumulate FT of subbands
        accum_FT_subbands = ...
            accum_FT_subbands + fft2(filtcoch(:,:,i,j)) .* conj(Hts);
        
        % accumulate FT of transfer functions
        accum_FT_transfer_function = accum_FT_transfer_function + Hts .* conj(Hts);
        
    end
end

% divide by accumulated transfer functions
coch = real(ifft2(accum_FT_subbands ./ accum_FT_transfer_function));
