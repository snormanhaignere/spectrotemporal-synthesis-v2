function filtcoch = coch2filtcoch_allsubbands(coch, P, complex_filters, fourier_domain)

% Computes filtered cochleograms for all of the model filters.
% 
% 2017-06-30: Modified the function to return both negative and positive rates

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
    length(P.spec_mod_rates), length(P.temp_mod_rates), 2);
for i = 1:length(P.spec_mod_rates)
    for j = 1:length(P.temp_mod_rates)
        sgn = [-1 1];
        for k = 1:2
            Hts = filt_spectemp_mod(...
                P.spec_mod_rates(i), sgn(k) * P.temp_mod_rates(j), ...
                size(coch,2), size(coch,1), P, P.spec_mod_lowpass(i), ...
                P.temp_mod_lowpass(j), 0, 0, complex_filters, [], P.causal);
            if fourier_domain
                filtcoch(:,:,i,j,k) = FT_coch .* Hts;
            else
                filtcoch(:,:,i,j,k) = ifft2(FT_coch .* Hts);
            end
        end
    end
end

% ensure real (only needed because of numerical issues)
if ~complex_filters && ~fourier_domain
    filtcoch = real(filtcoch);
end