function C = third_layer_cov(FT_filtcoch_complex, P, ti)

% Calculates the covariance of third layer temporal filters for a collection of
% second-layer filtered cochleograms
% 
% 2017-05-24: Created, Sam NH

% dimensionality of second-layer representation
% time x frequency x all-filters
[T, F, ~] = size(FT_filtcoch_complex);

% temporal indices over which to compute covariances
% useful for handling padding
if nargin < 3
    ti = 1:size(FT_filtcoch_complex,1);
end

C = cell(1, length(P.temp_mod_third_layer));
for i = 1:length(P.temp_mod_third_layer)
    
    % complex wavelet
    complex_filters = true;
    FT_wavelet = filt_spectemp_mod(...
        NaN, P.temp_mod_third_layer(i), F, T, P, 0, 0, 0, 0, complex_filters);
    
    % select filters to apply wavelet to
    xi = abs(P.temp_mod_to_match) > P.temp_mod_third_layer(i) ...
        & abs(P.temp_mod_to_match) > 0 & abs(P.spec_mod_to_match) > 0 ...
        & ~isnan(P.temp_mod_to_match) & ~isnan(P.spec_mod_to_match);
    n_filters = sum(xi);
    
    % apply wavelet to these filters
    FT_third_layer = bsxfun(@times, FT_filtcoch_complex(:,:,xi), FT_wavelet);
    clear xi;
    
    % convert back to time domain, still complex
    third_layer = ifft2(FT_third_layer);
    
    % reshape
    third_layer = reshape(third_layer, T, F*n_filters);
    
    % covariance
    C{i} = third_layer(ti,:)' * third_layer(ti,:);
    
end

