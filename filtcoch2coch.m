function coch = filtcoch2coch(filtcoch, P)

% total number of filters
assert(length(P.spec_mod_to_match) == length(P.temp_mod_to_match));
n_filters = length(P.spec_mod_to_match);

% compute the subbands
accum_FT_transfer_function = zeros(size(filtcoch,1), size(filtcoch,2));
accum_FT_subbands = zeros(size(filtcoch,1), size(filtcoch,2));
for i = 1:n_filters
    
    % filter transfer functions
    Hts = filt_spectemp_mod(P.spec_mod_to_match(i), P.temp_mod_to_match(i), ...
        size(filtcoch,2), size(filtcoch,1), P);
    
    % accumulate FT of subbands
    accum_FT_subbands = ...
        accum_FT_subbands + fft2(filtcoch(:,:,i)) .* conj(Hts);
    
    % accumulate FT of transfer functions
    accum_FT_transfer_function = accum_FT_transfer_function + Hts .* conj(Hts);
    
end

% divide by accumulated transfer functions
coch = real(ifft2(accum_FT_subbands ./ accum_FT_transfer_function));
