function [pca_timecourses, pca_weights, pca_eigvals] = ...
    third_layer_pca_multistim(filtcoch_MAT_files, nPC, P, output_directory, varargin)

% Computes third layer PCs from the envelopes of filtered cochleograms for
% multiple stimuli. The filtered cochleograms are stored in MAT files, one per
% stimulus, and the file names are given as input to this function. 
% 
% 2017-06-01: Created, Sam NH

I.overwrite = false;
I = parse_optInputs_keyvalue(varargin, I);

% number of stimuli and third layer filters
n_stimuli = length(filtcoch_MAT_files);
n_third_layer_filters = length(P.temp_mod_third_layer);

pca_third_layer_MAT_file = [output_directory '/third-layer-pca.mat'];
if ~exist(pca_third_layer_MAT_file, 'file') || I.overwrite
        
    % initialize
    pca_weights = cell(1, n_third_layer_filters);
    pca_eigvals = cell(1, n_third_layer_filters);
    pca_timecourses = cell(1, n_third_layer_filters);
    third_layer_MAT_files = cell(n_stimuli, n_third_layer_filters);
    
    for j = 1:n_third_layer_filters
        
        fprintf('Filter %d\n', j); drawnow;
        
        for i = 1:n_stimuli
            
            fprintf('Stim %d\n', i); drawnow;
                        
            % load filtered cochleogram
            load(filtcoch_MAT_files{i}, 'F', 'ti');
            filtcoch = F; clear F;
            temporal_indices = ti; clear ti;
            
            % compute envelopes if complex
            filtcoch = abs(filtcoch);
            
            % select relevant modulation rates
            ti = abs(P.temp_mod_rates) > P.temp_mod_third_layer(j) ...
                & abs(P.temp_mod_rates) > 0 & ~isnan(P.temp_mod_rates);
            si = abs(P.spec_mod_rates) > 0 & ~isnan(P.spec_mod_rates);
            filtcoch = filtcoch(:,:,si,ti);
            
            % 2D fourier transform of filtered envelopes
            FT_filtcoch = nan(size(filtcoch));
            for k = 1:size(filtcoch,3)
                for l = 1:size(filtcoch,4)
                    FT_filtcoch(:,:,k,l) = fft2(filtcoch(:,:,k,l));
                end
            end
            
            % dimensionality of second-layer representation
            [T, F, ~, ~] = size(filtcoch);
            
            % complex wavelet
            complex_filters = true;
            FT_wavelet = filt_spectemp_mod(...
                NaN, P.temp_mod_third_layer(j), F, T, P, 0, 0, 0, 0, complex_filters);
            
            % apply wavelet to these filters
            FT_third_layer = bsxfun(@times, FT_filtcoch, FT_wavelet);
            
            % convert back to time domain, still complex
            third_layer = ifft2(FT_third_layer);
            
            % select particular timepoints
            third_layer = third_layer(temporal_indices,:,:,:);
            
            % reshape to time x filter
            dims = size(third_layer);
            third_layer = reshape(third_layer, dims(1), prod(dims(2:end)));
            
            % save
            [~,fname,~] = fileparts(filtcoch_MAT_files{i});
            third_layer_MAT_files{i,j} = [output_directory '/' ...
                'third-layer-pca-' fname '-tempmod' num2str(P.temp_mod_third_layer(j)) 'Hz.mat'];
            save(third_layer_MAT_files{i,j}, 'third_layer', 'dims');
            
            % covariance
            if i == 1
                C = zeros(size(prod(dims(2:end))));
            end
            C = C + third_layer' * third_layer;
        end
        
        % eigen vector decomposition
        [pca_weights{j}, pca_eigvals{j}] = eig(C);
        
        % sort eigenvectors by their eigen value
        pca_eigvals{j} = diag(pca_eigvals{j});
        [~,xi] = sort(pca_eigvals{j}, 'descend');
        pca_weights{j} = pca_weights{j}(:,xi);
        pca_eigvals{j} = pca_eigvals{j}(xi);
        
        % truncate eigen vectors
        if size(pca_weights{j},2) > nPC
            pca_weights{j} = pca_weights{j}(:,1:nPC);
            pca_eigvals{j} = pca_eigvals{j}(1:nPC);
            ncomp = nPC;
        else
            ncomp = size(pca_weights{j},2);
        end
        
        % initialize
        if i == 1
            pca_timecourses{j} = nan(size(third_layer,1), n_stimuli, ncomp);
        end
        
        % project onto PCs
        for i = 1:n_stimuli
            load(third_layer_MAT_files{i,j}, 'third_layer');
            pca_timecourses{j}(:,i,:) = third_layer * pca_weights{j};
        end
        
        
    end
    
    save(pca_third_layer_MAT_file, 'pca_timecourses', 'pca_weights', 'pca_eigvals');
    
else
    
    load(pca_third_layer_MAT_file, 'pca_timecourses', 'pca_weights', 'pca_eigvals');

end