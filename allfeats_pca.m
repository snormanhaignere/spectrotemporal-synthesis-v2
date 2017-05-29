function allfeats_pca(stimuli, input_directory, output_directory, P, nPC)

% Computes all of the features of the spectrotemporal model in a relatively
% efficient and memory-saving manner for a set of stimuli.

%% Testing arguments

% input_directory = '/Users/svnh2/Desktop/projects/naturalsound-ecog/stims/naturalsounds165';
% output_directory = '/Users/svnh2/Desktop/projects/naturalsound-ecog/analysis/stimulus-acoustics/spectrotemporal-envelopes-v2';
% addpath(genpath(fileparts(which('allfeats.m'))));
% stimuli = mydir(input_directory);
% nPC = 10;
% P = synthesis_parameters_toy;
% P.temp_mod_rates = [2 8 32];
% P.spec_mod_rates = [1 2 4];
% P.temp_pad_sec = 1;
% P.freq_pad_oct = 1;
% P.lowrate_tempfilts_flat_spec = [];
% P.lowrate_tempfilts_impulse_spec = [];
% P.match_temp_mod = false;
% P.match_spec_mod = false;
% P.match_spectemp_mod = true;
% P.n_third_layer_PCs = 10;
% P.temp_mod_third_layer = [0.5 1 2 4 8 16];
% P = determine_filters_to_match(P);

%% Optional arguments

% optional parameters
I.overwrite = true;
I = parse_optInputs_keyvalue(varargin, I);

%% Stimulus-specific features, and feature covariances

% create output directory if it doesn't exist
if ~exist(output_directory, 'dir')
    mkdir(output_directory);
end

% loop through stimuli
for i = 1:length(stimuli)
    
    fprintf('%s\n', stimuli{i});
    
    % read in stimulus
    [y, sr] = audioread([input_directory '/' stimuli{i}]);
    
    % convert to mono if necessary
    if size(y,2) > 1
        y = mean(y,2);
    end
    
    % resample
    if P.audio_sr ~= sr
        y = resample(y, P.audio_sr, sr);
    end
        
    % compute cochleogram
    coch_MAT_file = [output_directory '/coch_' stimuli{i} '.mat'];
    if ~exist(coch_MAT_file, 'file') || I.overwrite
        [coch, P] = wav2coch_without_filts(y, P);
        
        % save
        F = coch;
        ti = (1:size(F,1));
        save(coch_MAT_file, 'F', 'ti', 'P');
        clear F ti;
    end
        
    % compute second layer filters from cochleogram
    filtcoch_MAT_file = [output_directory '/filtcoch_' stimuli{i} '.mat'];
    if ~exist(filtcoch_MAT_file, 'file') || I.overwrite
        
        % load cochleogram
        if ~exist('coch', 'var')
            load(coch_MAT_file, 'F', 'P');
            coch = F; clear F;
        end
        
        % pad cochleogram
        [padded_coch, ti, fi] = pad_coch(coch,P);
        
        % filtered cochleogram
        complex_filters = true;
        fourier_domain = false;
        filtcoch = coch2filtcoch_allsubbands(...
            padded_coch, P, complex_filters, fourier_domain);   
        clear complex_filters fourier_domain;
        
        % remove frequency padding
        filtcoch = filtcoch(:,fi,:,:);
        clear fi;
        
        % take modulus
        filtcoch = abs(filtcoch);
        
        % save
        F = filtcoch;
        save(filtcoch_MAT_file, 'F', 'ti', 'P');
        clear F;
        
    end
        
    % third layer
    third_layer_PCs_MAT_file = [output_directory '/third-layer-PCs_' stimuli{i} '.mat'];
    if ~exist(third_layer_MAT_file, 'file') || I.overwrite
        
        % load filtered cochleograms
        if ~exist('filtcoch', 'var') || ~exist('ti', 'var')
            load(filtcoch_MAT_file, 'filtcoch', 'ti');
        end
        
        % third layer PCs
        [third_layer_pca_timecourses, third_layer_pca_weights] = ...
            third_layer_pca(filtcoch, P);
        
        % concatenate different rates
        third_layer_pca_timecourses = cat(2, third_layer_pca_timecourses{:});
        
        % take modulus
        third_layer_pca_timecourses = abs(third_layer_pca_timecourses);
        
        % save
        F = third_layer_pca_timecourses;
        W = third_layer_pca_weights;
        save(third_layer_PCs_MAT_file, 'F', 'W', 'ti');
        clear F W ti;
        
    end
   
end

%% Feature covariances across stimuli

model_features = {'coch', 'filtcoch', 'third-layer-PCs'};
C = cell(1, length(model_features));
for i = 1:length(model_features)
    for j = 1:length(stimuli)
        load([output_directory '/' model_features{i} '_' stimuli{j} '.mat'], 'F', 'ti');
        dims = size(F);
        F = reshape(F, dims(1), prod(dims(2:end)));
        if j == 1
            C{i} = zeros(size(F,2));
        end
        F = F(ti,:);
        C{i} = C{i} + F' * F;
    end
end

%% Eigen decomposition of covariances

eigvecs = cell(size(model_features));
for i = 1:length(model_features)
    [eigvecs{i}, ~] = eig(C{i});
    eigvecs{i} = fliplr(eigvecs{i});
    eigvecs{i} = eigvecs{i}(:,1:nPC);
end

%% Calculate PCA activations

model_features = {'coch', 'filtcoch', 'third-layer-PCs'};
pca_activations = cell(1, length(model_features));
for i = 1:length(model_features)
    for j = 1:length(stimuli)
        load([output_directory '/' model_features{i} '_' stimuli{j} '.mat'], 'F', 'ti');
        dims = size(F);
        F = reshape(F, dims(1), prod(dims(2:end)));
        pca_activations{i}(:,j,:) = F(ti,:) * eigvecs{i};
    end
end


%%
% 
% imagesc(zscore(reshape(pca_activations{3}(:,:,:), [200*165,10])))
% % imagesc(squeeze(mean(pca_activations{3}(:,:,2:10),1)))
% 
% 
% 
% %% PCA activations
% 
% coch_activations = nan(200, 10, length(stimuli));
% filtcoch_activations = nan(200, 10, length(stimuli));
% third_layer_activations = nan(200, 10, length(stimuli));
% for i = 1:length(stimuli)
%     
%     % cochlear activations
%     load([output_directory '/coch_' stimuli{i} '.mat'], 'coch');
%     coch_activations(:,:,i) = coch * pinv(eigvecs');
%     
%     % filtered cochleogram activations
%     load([output_directory '/filtcoch_' stimuli{i} '.mat'], 'filtcoch', 'ti_filtcoch');
%     [n_t, n_freq, n_filters] = size(filtcoch);
%     filtcoch = reshape(filtcoch, [n_t, n_freq * n_filters]);
%     filtcoch_activations(:,:,i) = filtcoch(ti, :) * pinv(filtcoch_eigvecs');
%     filtcoch = reshape(filtcoch, [n_t, n_freq, n_filters]);
%     
%     % third layer activations
%     load([output_directory '/third-layer-PCs_' stimuli{i} '.mat'], ...
%         'third_layer_pca_timecourses', 'ti_filtcoch');
%     
%     % activations
%     third_layer_activations(:,:,i) = ...
%         third_layer_pca_timecourses(ti, :) * third_layer_eigvecs;
%     
% end


%% Covariance

% 
%     % covariance of the envelopes of the filtered cochleograms
%     cov_filtcoch_MAT_file = [output_directory '/cov-filtcoch_' stimuli{i} '.mat'];
%     if ~exist(cov_filtcoch_MAT_file, 'file') || I.overwrite
%         
%         % load filtered cochleograms
%         if ~exist('filtcoch', 'var') || ~exist('ti', 'var')
%             load(filtcoch_MAT_file, 'filtcoch', 'ti');
%         end
%        
%         % remove padded timepoints
%         X = filtcoch(ti, :, :, :);
%         
%         % unwrap and compute covariance
%         [n_t, n_freq, n_filters] = size(X);
%         X = reshape(X, [n_t, n_freq * n_filters]);
%         C = reshape(X' * X, n_freq, n_filters);
%         clear n_t n_freq n_filters X;
%         
%         % save results
%         save(cov_filtcoch_MAT_file, 'C');
%         
%     else
%         
%         load(cov_filtcoch_MAT_file, 'C');
%         
%     end

% % covariance of cochleogram
% cov_coch_MAT_file = [output_directory '/cov-coch_' stimuli{i} '.mat'];
% if ~exist(cov_coch_MAT_file, 'file') || I.overwrite
%     if ~exist('coch', 'var')
%         load(coch_MAT_file, 'F', 'P');
%         coch = F; clear F;
%     end
%     
%     % covariance
%     C = coch' * coch;
%     save(cov_coch_MAT_file, 'C');
%     clear C;
% end
% 
%     % covariance of the third-layer filters
%     cov_third_layer_PCs_MAT_file = [output_directory '/cov-third-layer-PCs_' stimuli{i} '.mat'];
%     if ~exist(cov_third_layer_PCs_MAT_file, 'file') || I.overwrite
%         
%         % load filtered cochleograms
%         if ~exist('third_layer_pca_timecourses', 'var')
%             load(third_layer_MAT_file, 'third_layer_pca_timecourses');
%         end
%         C_third_layer = third_layer_pca_timecourses(ti, :)' ...
%             * third_layer_pca_timecourses(ti, :);
%         save(cov_third_layer_PCs_MAT_file, 'C_third_layer');
%         
%     else
%         
%         load(cov_third_layer_PCs_MAT_file, 'C_third_layer');
%         
%     end
%     
% % accumulate cochlear covariance matrices
% if i == 1
%     C_coch_all_stimuli = zeros(size(C));
% end
% C_coch_all_stimuli = C_coch_all_stimuli + C;
% clear C_coch;
% 
% % accumulate covariance matrices
% if i == 1
%     C_filtcoch_all_stimuli = zeros(size(C_filtcoch));
% end
% C_filtcoch_all_stimuli = C_filtcoch_all_stimuli + C_filtcoch;
% clear C_filtcoch;
% 
% % accumulate covariance matrices
% if i == 1
%     C_third_layer_all_stimuli = zeros(size(C_third_layer));
% end
% C_third_layer_all_stimuli = C_third_layer_all_stimuli + C_third_layer;
% clear C_third_layer;
% 
% 
% 
% %% Activations
% 
% coch_activations = nan(200, 10, length(stimuli));
% filtcoch_activations = nan(200, 10, length(stimuli));
% third_layer_activations = nan(200, 10, length(stimuli));
% for i = 1:length(stimuli)
%     
%     % cochlear activations
%     load([output_directory '/coch_' stimuli{i} '.mat'], 'coch');
%     coch_activations(:,:,i) = coch * pinv(eigvecs');
%     
%     % filtered cochleogram activations
%     load([output_directory '/filtcoch_' stimuli{i} '.mat'], 'filtcoch', 'ti_filtcoch');
%     [n_t, n_freq, n_filters] = size(filtcoch);
%     filtcoch = reshape(filtcoch, [n_t, n_freq * n_filters]);
%     filtcoch_activations(:,:,i) = filtcoch(ti, :) * pinv(filtcoch_eigvecs');
%     filtcoch = reshape(filtcoch, [n_t, n_freq, n_filters]);
%     
%     % third layer activations
%     load([output_directory '/third-layer-PCs_' stimuli{i} '.mat'], ...
%         'third_layer_pca_timecourses', 'ti_filtcoch');
%     
%     % activations
%     third_layer_activations(:,:,i) = ...
%         third_layer_pca_timecourses(ti, :) * third_layer_eigvecs;
%     
% end

% coch_activations = permute(coch_activations, [1 3 2]);
% third_layer_activations = permute(third_layer_activations, [1 3 2]);
% filtcoch_activations = permute(filtcoch_activations, [1 3 2]);


% 
% imagesc(third_layer_activations(:,:,6))

