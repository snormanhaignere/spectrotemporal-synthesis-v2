function filtcoch = coch2filtcoch(coch, spec_mod_rate, temp_mod_rate, P)

% time x frequency
[T,F] = size(coch);

% 2D fourier transforms
FT_coch = fft2(coch);

% filter transfer function
Hts = filt_spectemp_mod(spec_mod_rate, temp_mod_rate, F, T, P);

% apply filter and revert to cochleogram domain
filtcoch = real(ifft2(FT_coch .* Hts));