function filtcoch = coch2filtcoch(coch, spec_mod_rate, temp_mod_rate, P, complex_filters)

if nargin < 5
    complex_filters = false;
end

% time x frequency
[T,F] = size(coch);

% 2D fourier transforms
FT_coch = fft2(coch);

% filter transfer function
Hts = filt_spectemp_mod(spec_mod_rate, ...
    temp_mod_rate, F, T, P, 0, 0, 0, 0, complex_filters);

% apply filter and revert to cochleogram domain
filtcoch = ifft2(FT_coch .* Hts);

% filtered cochleogram
if ~complex_filters
    filtcoch = real(filtcoch);
end