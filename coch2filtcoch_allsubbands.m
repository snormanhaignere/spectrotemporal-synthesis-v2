function filtcoch = coch2filtcoch_allsubbands(coch, P, complex_filters, fourier_domain, causal)

if nargin < 3
    complex_filters = false;
end

if nargin < 4
    fourier_domain = false;
end

% FT of the cochleogram
FT_coch = fft2(coch);

% compute the subbands
filtcoch = nan(size(coch,1), size(coch,2), ...
    length(P.spec_mod_rates), length(P.temp_mod_rates));
for i = 1:length(P.spec_mod_rates)
    for j = 1:length(P.temp_mod_rates)
        Hts = filt_spectemp_mod(P.spec_mod_rates(i), P.temp_mod_rates(j), ...
            size(coch,2), size(coch,1), P, P.spec_mod_lowpass(i), ...
            P.temp_mod_lowpass(j), 0, 0, complex_filters, [], causal);
        if fourier_domain
            filtcoch(:,:,i,j) = FT_coch .* Hts;
        else
            filtcoch(:,:,i,j) = ifft2(FT_coch .* Hts);
        end
    end
end

% ensure real (only needed because of numerical issues)
if ~complex_filters && ~fourier_domain
    filtcoch = real(filtcoch);
end