function filtcoch = coch2filtcoch_allsubbands(coch, P, complex_filters)

if nargin < 3
    complex_filters = false;
end

% FT of the cochleogram
FT_coch = fft2(coch);

% total number of filters
assert(length(P.spec_mod_to_match) == length(P.temp_mod_to_match));
n_filters = length(P.spec_mod_to_match);

% compute the subbands
filtcoch = nan(size(coch,1), size(coch,2), n_filters);
for i = 1:n_filters
    Hts = filt_spectemp_mod(P.spec_mod_to_match(i), P.temp_mod_to_match(i), ...
        size(coch,2), size(coch,1), P, 0, 0, 0, 0, complex_filters);
    filtcoch(:,:,i) = ifft2(FT_coch .* Hts);
end

% ensure real (only needed because of numerical issues)
if ~complex_filters
    filtcoch = real(filtcoch);
end