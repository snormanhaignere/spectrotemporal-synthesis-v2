function [pca_timecourses, pca_weights, pca_eigvals] = ...
    third_layer_pca_multistim(filtcoch_MAT_files, nPC, P, output_directory, varargin)

% Computes third layer PCs from the envelopes of filtered cochleograms for
% multiple stimuli. The filtered cochleograms are stored in MAT files, one per
% stimulus, and the file names are given as input to this function.
%
% 2017-06-01: Created, Sam NH
%
% 2017-06-13: Added option to demean/standard features before PCA
% 
% 2017-06-16: Added spectral filters and possibility of demeaning/standardizing

I.overwrite = false;
% I.demean_feats = false;
% I.std_feats = false;
I = parse_optInputs_keyvalue(varargin, I);

% number of stimuli and third layer filters
n_stimuli = length(filtcoch_MAT_files);
n_temp_mod_filters = length(P.temp_mod_third_layer);
n_spec_mod_filters = length(P.spec_mod_third_layer);

pca_third_layer_MAT_file = [output_directory '/third-layer-pca' ...
    '-demean' num2str(I.demean_feats) '-std' num2str(I.std_feats) '.mat'];
if ~exist(pca_third_layer_MAT_file, 'file') || I.overwrite
    
    % initialize
    pca_weights = cell(n_spec_mod_filters, n_temp_mod_filters);
    pca_eigvals = cell(n_spec_mod_filters, n_temp_mod_filters);
    pca_timecourses = cell(n_spec_mod_filters, n_temp_mod_filters);
    third_layer_MAT_files = cell(n_stimuli, n_spec_mod_filters, n_temp_mod_filters);
    
    for i = 1:n_spec_mod_filters
        for j = 1:n_temp_mod_filters
            
            fprintf('%.1f cyc/oct, %d Hz\n', ...
                P.spec_mod_third_layer(i), P.temp_mod_third_layer(j)); drawnow;
            
            N_all = 0;
            for k = 1:n_stimuli
                
                fprintf('Stim %d\n', k); drawnow;
                
                % load filtered cochleogram
                load(filtcoch_MAT_files{k}, 'F', 'ti');
                filtcoch = F; clear F;
                temporal_indices = ti; clear ti;
                
                % compute envelopes if complex
                filtcoch = abs(filtcoch);
                
                % select relevant modulation rates
                si = abs(P.spec_mod_rates) > P.spec_mod_third_layer(i);
                ti = abs(P.temp_mod_rates) > P.temp_mod_third_layer(j);
                filtcoch = filtcoch(:,:,si,ti);
                clear si ti;
                
                % 2D fourier transform of filtered envelopes
                FT_filtcoch = nan(size(filtcoch));
                for m = 1:size(filtcoch,3)
                    for l = 1:size(filtcoch,4)
                        FT_filtcoch(:,:,m,l) = fft2(filtcoch(:,:,m,l));
                    end
                end
                
                % dimensionality of second-layer representation
                [T, F, ~, ~] = size(filtcoch);
                
                % complex wavelet
                complex_filters = true;
                FT_wavelet = filt_spectemp_mod(...
                    P.spec_mod_third_layer(i), P.temp_mod_third_layer(j), ...
                    F, T, P, 0, 0, 0, 0, complex_filters);
                
                % apply wavelet to these filters
                FT_third_layer = bsxfun(@times, FT_filtcoch, FT_wavelet);
                
                % convert back to time domain, still complex
                third_layer = ifft2(FT_third_layer);
                
                % select particular timepoints (e.g. to remove padding)
                third_layer = third_layer(temporal_indices,:,:,:);
                
                % reshape to time x filter
                dims = size(third_layer);
                third_layer = reshape(third_layer, dims(1), prod(dims(2:end)));
                
                % save
                [~,fname,~] = fileparts(filtcoch_MAT_files{k});
                third_layer_MAT_files{k,i,j} = ...
                    [output_directory '/third-layer-pca-' fname ...
                    '-specmod' num2str(P.spec_mod_third_layer(i)) ...
                    '-tempmod' num2str(P.temp_mod_third_layer(j)) 'Hz.mat'];
                save(third_layer_MAT_files{k,i,j}, 'third_layer', 'dims');
                
                % accumulate means and covariances
                if k == 1
                    C = zeros(prod(dims(2:end)));
                    M = zeros(1, prod(dims(2:end)));
                end
                C = C + third_layer' * third_layer;
                M = M + sum(third_layer,1);
                N_all = N_all + size(third_layer,1);
            end
            
            % normalize by number of samples
            C = C/N_all;
            M = M/N_all;
            
            keyboard;
            
            %             % standard deviations
            %             S = sqrt(diag(C)' - conj(M).*M);
            %             assert(all(isreal(S(:))));
            %
            %             % optionally remove effects of means and standard deviations from
            %             % correlation matrix
            %             if I.demean_feats
            %                 CV = C - M' * M;
            %             else
            %                 CV = C;
            %             end
            %             if I.std_feats
            %                 CV = (1./(S' * S)) .* CV;
            %             end
            CV = C;
            clear C;
            
            % eigen vector decomposition
            [pca_weights{i,j}, pca_eigvals{i,j}] = eig(CV);
            
            % sort eigenvectors by their eigen value
            pca_eigvals{i,j} = diag(pca_eigvals{i,j});
            [~,xi] = sort(pca_eigvals{i,j}, 'descend');
            pca_weights{i,j} = pca_weights{i,j}(:,xi);
            pca_eigvals{i,j} = pca_eigvals{i,j}(xi);
            
            % truncate eigen vectors
            if size(pca_weights{i,j},2) > nPC
                pca_weights{i,j} = pca_weights{i,j}(:,1:nPC);
                pca_eigvals{i,j} = pca_eigvals{i,j}(1:nPC);
                ncomp = nPC;
            else
                ncomp = size(pca_weights{i,j},2);
            end
            
            % initialize
            if k == 1
                pca_timecourses{i,j} = nan(size(third_layer,1), n_stimuli, ncomp);
            end
            
            % project onto PCs
            for k = 1:n_stimuli
                load(third_layer_MAT_files{k,i,j}, 'third_layer');
                %                 if I.demean_feats
                %                     third_layer = bsxfun(@minus, third_layer, M);
                %                 end
                %                 if I.std_feats
                %                     third_layer = bsxfun(@times, third_layer, 1./S);
                %                 end
                pca_timecourses{i,j}(:,k,:) = third_layer * pca_weights{i,j};
            end
        end
    end
    save(pca_third_layer_MAT_file, 'pca_timecourses', 'pca_weights', 'pca_eigvals');
    
else
    
    load(pca_third_layer_MAT_file, 'pca_timecourses', 'pca_weights', 'pca_eigvals');
    
end