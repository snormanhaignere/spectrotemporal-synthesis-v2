function [third_layer_pca_timecourses, third_layer_pca_weights] = ...
    third_layer_pca(filtcoch, P)

% dimensionality of second-layer representation
[T, F, ~, ~] = size(filtcoch);

% 2D fourier transform of filtered envelopes
FT_filtcoch = nan(size(filtcoch));
for i = 1:length(P.spec_mod_rates)
    for j = 1:length(P.temp_mod_rates)
        FT_filtcoch(:,:,i,j) = ifft2(filtcoch(:,:,i,j));
    end
end

% wavelet in fourier domain
third_layer_pca_timecourses = cell(1, length(P.temp_mod_third_layer));
third_layer_pca_weights = cell(1, length(P.temp_mod_third_layer));
for i = 1:length(P.temp_mod_third_layer)
    
    % complex wavelet
    complex_filters = true;
    FT_wavelet = filt_spectemp_mod(...
        NaN, P.temp_mod_third_layer(i), F, T, P, 0, 0, 0, 0, complex_filters);

    % select filters to apply wavelet to
    ti = abs(P.temp_mod_rates) > P.temp_mod_third_layer(i) ...
        & abs(P.temp_mod_rates) > 0 & ~isnan(P.temp_mod_rates);
    si = abs(P.spec_mod_rates) > 0 & ~isnan(P.spec_mod_rates);
        
    % apply wavelet to these filters
    FT_third_layer = bsxfun(@times, FT_filtcoch(:,:,si,ti), FT_wavelet);

    % convert back to time domain, still complex
    third_layer = ifft2(FT_third_layer);
    
    % unwrap frequency/filter dimensions, and apply PCA via svd
    dims = size(third_layer);
    [U, S, V] = svd(reshape(third_layer, dims(1), prod(dims(2:end))), 'econ');
    
    % truncate
    U = U(:, 1:P.n_third_layer_PCs);
    S = S(1:P.n_third_layer_PCs, 1:P.n_third_layer_PCs);
    V = V(:, 1:P.n_third_layer_PCs);
    
    % scale timecourses by variances
    third_layer_pca_timecourses{i} = U*S;
    
    % reshape weights
    third_layer_pca_weights{i} = reshape(V, [dims(2:end), P.n_third_layer_PCs]);
    
end