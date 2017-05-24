function [third_layer_pca_timecourses, third_layer_pca_weights] = ...
    third_layer_pca(FT_filtcoch_complex, P, n_comp, third_layer_pca_weights)

% dimensionality of second-layer representation
% time x frequency x all-filters
[T, F, ~] = size(FT_filtcoch_complex);

% initialize pca weights and timecourses
third_layer_pca_timecourses = cell(1, length(P.temp_mod_third_layer));
if nargin < 4 || isempty(third_layer_pca_weights)
    compute_PCs = true;
    third_layer_pca_weights = cell(1, length(P.temp_mod_third_layer));
else
    compute_PCs = false;
end

% wavelet in fourier domain
% for i = 1:length(P.temp_mod_third_layer)
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
    
    if compute_PCs
        
        % unwrap frequency/filter dimensions, and apply PCA via svd
        [U, S, V] = svd(reshape(third_layer, T, F*n_filters), 'econ');
        
        % truncate
        U = U(:, 1:n_comp);
        S = S(1:n_comp, 1:n_comp);
        V = V(:, 1:n_comp);
        
        % scale timecourses by variances
        third_layer_pca_timecourses{i} = U*S;
        
        % reshape weights
        third_layer_pca_weights{i} = V;
        
    else % apply weights from a previous analysis
        
        % estimate weights
        third_layer_pca_timecourses{i} = ...
            reshape(third_layer, T, F*n_filters) * third_layer_pca_weights{i};
        
    end
    
end