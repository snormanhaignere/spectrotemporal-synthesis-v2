function [pca_timecourses, model_features] = allfeats_pca(...
    stimuli, input_directory, output_directory, P, nPC, varargin)

% Computes all of the features of the spectrotemporal model in a relatively
% efficient and memory-saving manner for a set of stimuli.
% 
% 2017-06-13: Added option to demean/standard features before PCA

%% Testing arguments

% input_directory = '/Users/svnh2/Desktop/projects/naturalsound-ecog/stims/naturalsounds165';
% output_directory = '/Users/svnh2/Desktop/projects/naturalsound-ecog/analysis/stimulus-acoustics/spectrotemporal-envelopes-v2';
% addpath(genpath(fileparts(which('allfeats.m'))));
% stimuli = mydir(input_directory);
% stimuli = stimuli(1:3);
% nPC = 10;
% varargin = {};
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
I.overwrite = false;
I.keyboard = false;
I.demean_feats = false;
I.std_feats = false;
I = parse_optInputs_keyvalue(varargin, I);
if I.keyboard
    keyboard;
end

%% Cochleograms and filtered cochleograms

% create output directory if it doesn't exist
if ~exist(output_directory, 'dir')
    mkdir(output_directory);
end

% loop through stimuli
for i = 1:length(stimuli)
    
    % separate file name from extension
    [~, fname, ~] = fileparts(stimuli{i});
    
    % print file name
    fprintf('\n\n%d: %s\n', i, fname); drawnow;
    
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
    coch_MAT_file = [output_directory '/coch_' fname '.mat'];
    if ~exist(coch_MAT_file, 'file') || I.overwrite
        fprintf('Computing cochleogram\n'); drawnow;
        [coch, P] = wav2coch_without_filts(y, P);
                
        % save
        F = coch;
        ti = (1:size(F,1));
        save(coch_MAT_file, 'F', 'ti', 'P');
        clear F ti;
    end
    
    % compute second layer filters from cochleogram
    filtcoch_MAT_file = [output_directory '/filtcoch_' fname '.mat'];
    if ~exist(filtcoch_MAT_file, 'file') || I.overwrite
        fprintf('Computing filtered cochleograms\n'); drawnow;
        
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
        filtcoch = filtcoch(:,fi,:,:,:);
        clear fi;
        
        % take modulus
        filtcoch = abs(filtcoch);
        
        % save
        F = filtcoch;
        save(filtcoch_MAT_file, 'F', 'ti', 'P');
        clear F;
        
    end
end

%% Principal components for cochleograms and filtered cochleograms

pca_activations_coch_filtcoch_MAT_file = ...
    [output_directory '/pca-activations-coch-filtcoch' ...
    '-demean' num2str(I.demean_feats) '-std' num2str(I.std_feats) '.mat'];

if ~exist(pca_activations_coch_filtcoch_MAT_file, 'file') || I.overwrite
        
    model_features = {'coch', 'filtcoch'};
    
    % means, correlations, standard deviations and covariances
    C = cell(1, length(model_features));
    M = cell(1, length(model_features));
    S = cell(1, length(model_features));
    CV = cell(1, length(model_features));
    for i = 1:length(model_features)
        N_all = 0;
        for j = 1:length(stimuli)
            % load and format
            [~,fname,~] = fileparts(stimuli{j});
            fprintf('feat cov: %s, %s\n', model_features{i}, stimuli{j}); drawnow;
            load([output_directory '/' model_features{i} '_' fname '.mat'], 'F', 'ti');
            dims = size(F);
            F = reshape(F, dims(1), prod(dims(2:end)));
            F = F(ti,:);
            
            % initialize
            if j == 1
                C{i} = zeros(size(F,2));
                M{i} = zeros(1,size(F,2));
            end
            
            % accumulate correlation, sum, and total number of samples
            C{i} = C{i} + F' * F;
            M{i} = M{i} + sum(F);
            N_all = N_all + size(F,1);
            if any(~isreal(C{i}))
                break
            end
        end
        
        % normalize by number of samples
        C{i} = C{i}/N_all;
        M{i} = M{i}/N_all;
        assert(all(isreal(C{i})));
        
        % standard deviation
        S{i} = sqrt(diag(C{i})' - M{i}.^2);
        
        % optionally remove effects of means and standard deviations from
        % correlation matrix
        if I.demean_feats
            CV{i} = C{i} - M{i}' * M{i};
        else
            CV{i} = C{i};
        end
        if I.std_feats
            CV{i} = (1./(S{i}' * S{i})) .* CV{i};
        else
            CV{i} = CV{i};
        end
        assert(all(isreal(CV{i})));
    end
        
    % pca weights on the filters form the covariance/correlation matrix
    pca_weights = cell(size(model_features));
    pca_eigvals = cell(size(model_features));
    for i = 1:length(model_features)
        
        % eigenvectors of covariance
        [pca_weights{i}, pca_eigvals{i}] = eig(CV{i});
        
        % sort by variance
        pca_eigvals{i} = diag(pca_eigvals{i});
        [~, xi] = sort(pca_eigvals{i}, 'descend');
        pca_weights{i} = pca_weights{i}(:,xi);
        pca_eigvals{i} = pca_eigvals{i}(xi);
        
        % select subset of components
        if size(pca_weights{i},2) > nPC
            pca_weights{i} = pca_weights{i}(:,1:nPC);
        end
        
    end
        
    % pca timecourses from the weights
    pca_timecourses = cell(size(model_features));
    for i = 1:length(model_features)
        for j = 1:length(stimuli)
            % load and format
            [~,fname,~] = fileparts(stimuli{j});
            fprintf('pca timecourses: %s, %s\n', model_features{i}, stimuli{j}); drawnow;
            load([output_directory '/' model_features{i} '_' fname '.mat'], 'F', 'ti');
            dims = size(F);
            F = reshape(F, dims(1), prod(dims(2:end)));
            F = F(ti,:);
            
            % optionally demean and standardize
            if I.demean_feats
                F = bsxfun(@minus, F, M{i});
            end
            if I.std_feats
                F = bsxfun(@times, F, 1./S{i});
            end
            
            % apply weights
            pca_timecourses{i}(:,j,:) = F * pca_weights{i};
        end
    end
        
    % save
    save(pca_activations_coch_filtcoch_MAT_file, 'pca_timecourses', 'pca_weights', 'pca_eigvals', 'model_features', '-v7.3');
    
else
    
    load(pca_activations_coch_filtcoch_MAT_file, 'pca_timecourses', 'model_features');
    
end
    

%% Third layer principal components

pca_activations_third_layer_MAT_file = ...
    [output_directory '/pca-activations-third-layer' ...
    '-demean' num2str(I.demean_feats) '.mat'];
if ~exist(pca_activations_third_layer_MAT_file, 'file') || I.overwrite
    
    filtcoch_MAT_files = cell(1, length(stimuli));
    for i = 1:length(stimuli)
        [~, fname, ~] = fileparts(stimuli{i});
        filtcoch_MAT_files{i} = [output_directory '/filtcoch_' fname '.mat'];
    end
    
    % principal components per rate
    [third_layer_timecourses, third_layer_weights, third_layer_eigvals] = ...
        third_layer_pca_multistim(filtcoch_MAT_files, nPC, P, output_directory, ...
        'overwrite', I.overwrite);
    
    % reshape to matrix
    third_layer_timecourses = cat(3, third_layer_timecourses{:});
    [n_timepoints, n_stimuli, n_comp] = size(third_layer_timecourses);
    third_layer_timecourses = ...
        reshape(third_layer_timecourses, n_timepoints * n_stimuli, n_comp);
            
    % compute envelopes
    modulus_third_layer_timecourses = abs(third_layer_timecourses);
    
    % demean envelopes
    if I.demean_feats
        modulus_third_layer_timecourses = bsxfun(@minus, ...
            modulus_third_layer_timecourses, mean(modulus_third_layer_timecourses));
    end
        
    % principal components of envelopes
    [U,S,V] = svd(modulus_third_layer_timecourses,  'econ');
    pca_modulus_third_layer_timecourses = U(:, 1:nPC) * S(1:nPC, 1:nPC);
    pca_modulus_third_layer_weights = V(:,1:nPC);
    
    % reshape timecourse
    pca_modulus_third_layer_timecourses = reshape(...
        pca_modulus_third_layer_timecourses, n_timepoints, n_stimuli, nPC);
    
    % save
    save(pca_activations_third_layer_MAT_file, ...
        'pca_modulus_third_layer_timecourses', 'pca_modulus_third_layer_weights', ...
        'third_layer_timecourses', 'third_layer_weights', 'third_layer_eigvals', '-v7.3');
    
else
    
    load(pca_activations_third_layer_MAT_file, ...
        'pca_modulus_third_layer_timecourses');
    
end

% add third layer to output
pca_timecourses = [pca_timecourses, pca_modulus_third_layer_timecourses];
model_features = [model_features, {'thirdlayer'}];
