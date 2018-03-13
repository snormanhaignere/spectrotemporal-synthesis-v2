function [pca_timecourse_MAT_files, pca_weight_MAT_files, model_features, ...
    pca_timecourses_allstim_allmodels, pca_weights, pca_eigvals] = allfeats_pca_multistim(...
    modulation_types, stimuli, nPCs, sr, input_directory, output_directory, P_orig)

% Computes PCAs of the spectrotemporal wavelet representation in a relatively
% efficient manner.
% 
% See the example below for how to use this script. Also look at
% measure_parameters_default.m and the documentation below for additional
% information about the various parameters.
% 
% % -- Example with low-res cochleogram and small number of filters --
% 
% % tempmod: only temporal modulation filters
% % specmod: only spectral modulation filters
% % spectempmod: joint spectrotemporal filters
% modulation_types = {'tempmod', 'specmod', 'spectempmod'};
% 
% % stimulus files and directory containing the stimuli
% stimuli = {'man_speaking.wav', 'orchestra_music.wav', 'walking_on_leaves.wav'};
% input_directory = [pwd '/stimuli'];
% 
% % directory to save results
% output_directory = [pwd '/feature_PCAs'];
% 
% % number of PCs to keep
% nPCs = 10;
% 
% % temporal sampling rate of the feautres in Hz
% sr = 50;
% 
% % parameters for highly reduced model
% % use measurement_parameters_default for actual analyses
% P = measurement_parameters_toy;
% P.temp_mod_rates = [2, 8];
% P.temp_mod_lowpass = [0, 0];
% P.spec_mod_rates = [0.25, 1];
% P.spec_mod_lowpass = [0, 0];
% P.causal = false;
% 
% % compute PCAs
% [pca_timecourse_MAT_files, pca_weight_MAT_files, model_features] = allfeats_pca_multistim(...
%     modulation_types, stimuli, nPCs, sr, input_directory, output_directory, P);
% 
% % Plot the PCA timecourses for the various features and stimuli
% figure;
% set(gcf, 'Position', [0 0 800 800]);
% for j = 1:length(model_features)
%     for i = 1:length(stimuli)
%         subplot(length(model_features), 3, i + (j-1)*3);
%         load(pca_timecourse_MAT_files{i, j}, 'pca_timecourses');
%         imagesc(pca_timecourses);
%         title(sprintf('%s\n%s',  strrep(model_features{j}, '_', ' '), ...
%             strrep(stimuli{i}, '_', ' ')));
%         xlabel('PC'); ylabel('Time / Samples');
%         set(gca, 'FontSize', 10);
%     end
% end
% 
% % load pca timecourses and weights for the first stimulus, run through the joint
% % spectrotemporal model with a modulus nonlinearity
% % note all of the weights for the PCs have been kept, in contrast with the
% % timecourses where only the top nPCs have been kept
% % pca_timecourses: time x components
% % pca_weights: features x components
% load(pca_timecourse_MAT_files{1, 4}, 'pca_timecourses', 'dims');
% load(pca_weight_MAT_files{1, 4}, 'pca_weights')
% 
% % reconstruct the original filtered spectrograms
% % time x frequency x spectral modulation x temporal modulation x orientation (up / down)
% F = reshape(pca_timecourses * pca_weights(:,1:10)', dims);
% 
% % plot the reconstructed spectrograms for just the first orientation
% figure;
% set(gcf, 'Position', [0 0 800 800]);
% for i = 1:length(P.temp_mod_rates)
%     for j = 1:length(P.spec_mod_rates)
%         subplot(length(P.temp_mod_rates), length(P.spec_mod_rates), j + (i-1)*length(P.spec_mod_rates));
%         imagesc(F(:,:,j,i,1));
%         title(sprintf('%.2f Hz, %.2f cyc/oct', P.temp_mod_rates(i), P.spec_mod_rates(j)));
%     end
% end
% 
% % plot eigen value spectrum
% load(pca_weight_MAT_files{1, 4}, 'pca_eigvals');
% figure;
% plot(pca_eigvals, '-o');
% xlabel('PCA'); ylabel('Eigen value');
% 
% 
% 2018-03-05: Started working on this script, copied much from
% twolayer_feats_multiple_stim.m

%% Optional arguments

addpath(genpath(fileparts(which(mfilename))));

% add the required parameters to the parameter structure
P_orig.nPCs = nPCs;
P_orig.sr = sr;
clear sr nPCs;

% nonlinearity to apply
% type of nonlinearity
% options are "modulus" and "real"
if ~isfield(P_orig, 'nonlin')
    P_orig.nonlin = 'modulus';
end

% stimuli with which to compute the PCA weights, useful for cross validation
% when you just want to compute the PCA weights using the training set
% by default all of the stimuli are used to compute the PCA weights
if ~isfield(P_orig, 'stim_to_compute_PC_weights')
    P_orig.stim_to_compute_PC_weights = stimuli;
end

% whether to return the PCA timecourses/weights as variables
% default is to save them in MAT files
if ~isfield(P_orig, 'return_PCA_timecourses')
    P_orig.return_PCA_timecourses = false;
end
    
% whether or not to demean or standarize features prior to computing PCs
% default is to both standardize and demean
if ~isfield(P_orig, 'demean_feats')
    P_orig.demean_feats = true;
end
if ~isfield(P_orig, 'std_feats')
    P_orig.std_feats = true;
end

% whether or not to overwrite previously computed results
% ** note that if you change parameters, you should overwrite the files or
% change the output directory **
% the only exception is the number of PCs and the modulation_types which you can
% safely vary without overwriting or changing the output directory
if ~isfield(P_orig, 'overwrite')
    P_orig.overwrite = false;
end

% whether or not to debug
if ~isfield(P_orig, 'keyboard')
    P_orig.keyboard = false;
end

% enter debug mode
if P_orig.keyboard
    keyboard;
end

%% Strings identifying the parameters of this analysis

coch_idstring = [...
    'coch' ...
    '_audsr' num2str(P_orig.audio_sr) ...
    '_nfilts' num2str(P_orig.n_filts) ...
    '_lofreq' num2str(P_orig.lo_freq_hz) ...
    '_over' num2str(P_orig.overcomplete) ...
    '_exp' num2str(P_orig.compression_factor) ...
    '_envsr' num2str(P_orig.env_sr) ...
    '_freqsr' num2str(1/P_orig.logf_spacing) ...
    '_tpad' num2str(P_orig.temp_pad_sec) ...
    '_outsr' num2str(P_orig.sr) ...
    ];

modulation_idstring = [...
    'mod' ...
    '_tprate' sprintf('-%.2f', P_orig.temp_mod_rates) ...
    '_tplow' sprintf('%d', P_orig.temp_mod_lowpass) ...
    '_sprate' sprintf('-%.2f', P_orig.spec_mod_rates) ...
    '_splow' sprintf('%d', P_orig.spec_mod_lowpass) ...
    '_fpad' num2str(P_orig.freq_pad_oct) ...
    '_caus' num2str(P_orig.causal) ...
    '_' coch_idstring
    ];

pca_idstring = [...
    'nPCs' num2str(P_orig.nPCs) ...
    '_demean' num2str(P_orig.demean_feats) ...
    '_std' num2str(P_orig.std_feats) ...
    ];

%% Cochleograms and filtered cochleograms

coch_output_directory = [output_directory '/' coch_idstring];
modulation_output_directory = [output_directory '/' modulation_idstring];
if ~exist(coch_output_directory, 'dir'); mkdir(coch_output_directory); end
if ~exist(modulation_output_directory, 'dir'); mkdir(modulation_output_directory); end

% loop through stimuli
for i = 1:length(stimuli)
    
    % separate file name from extension
    [~, fname, ~] = fileparts(stimuli{i});
    
    % MAT file to save cochleogram
    coch_MAT_file = [coch_output_directory '/coch_' fname '.mat'];

    % MAT file for modulation representations
    filtcoch_MAT_files = cell(1, length(modulation_types));
    for j = 1:length(modulation_types)
        filtcoch_MAT_files{j} = [modulation_output_directory '/' ...
            modulation_types{j} '_' P_orig.nonlin '_' fname '.mat'];
    end
    
    % check if the files already exist
    all_features_computed = true;
    files_to_check = [{coch_MAT_file}, filtcoch_MAT_files];
    for j = 1:length(files_to_check)
        if ~exist(files_to_check{j}, 'file');
            all_features_computed = false;
        end
    end
    
    % check if the features have already all been computed
    if all_features_computed && ~P_orig.overwrite
        continue;
    end
    
    % format cochleogram and save
    if ~exist(coch_MAT_file, 'file') || P_orig.overwrite
        
        % print file name
        fprintf('\n\n%d: %s\n', i, fname); drawnow;
        
        % read in stimulus
        [y, wav_sr] = audioread([input_directory '/' stimuli{i}]);
        
        % convert to mono if necessary
        if size(y,2) > 1
            y = mean(y,2);
        end
        
        % resample
        if P_orig.audio_sr ~= wav_sr
            y = resample(y, P_orig.audio_sr, wav_sr);
        end
        
        % pad in time
        wav_dur_sec = length(y) / P_orig.audio_sr;
        y = [y; zeros(P_orig.audio_sr * P_orig.temp_pad_sec,1)]; %#ok<AGROW>
        
        % cochleogram
        fprintf('Computing cochleogram\n'); drawnow;
        [coch, P_coch] = wav2coch_without_filts(y, P_orig);
        
        % remove temporal padding
        F = coch((1:wav_dur_sec*P_coch.env_sr), :);
        
        % resample
        if P_orig.sr ~= P_coch.env_sr
            F = resample_ndarray(F, P_orig.sr, P_coch.env_sr, 1);
        end
        
        % save
        P = P_coch;
        save(coch_MAT_file, 'F', 'P');
        clear F P;
        
    else
        
        % load
        save(coch_MAT_file, 'F', 'P');
        coch = F;
        P_coch = P;
        clear F P;

    end
    
    % compute second layer filters from cochleogram
    for j = 1:length(modulation_types)
        if ~exist(filtcoch_MAT_files{j}, 'file') || P_orig.overwrite
                           
            fprintf('Computing %s representation\n', modulation_types{j}); drawnow;
                        
            switch modulation_types{j}
                case 'specmod'
                    P_mod = P_coch;
                    P_mod.temp_mod_rates = NaN;
                    P_mod.temp_mod_lowpass = NaN;
                    separable = true;
                case 'tempmod'
                    P_mod = P_coch;
                    P_mod.spec_mod_rates = NaN;
                    P_mod.spec_mod_lowpass = NaN;
                    P_mod.freq_pad_oct = 0;
                    separable = true;
                case 'spectempmod'
                    P_mod = P_coch;
                    separable = false;
                case 'spectempmod_separable'
                    P_mod = P_coch;
                    separable = true;
                otherwise
                    error('No matching case');
            end
            
            % pad frequency
            n_freq_pad_smps = round(P_mod.freq_pad_oct / P_mod.logf_spacing);
            padded_coch = [coch, zeros(size(coch,1), n_freq_pad_smps)];

            % filtered cochleogram
            fourier_domain = false;
            switch P_orig.nonlin
                case 'modulus'
                    complex_filters = true;
                    filtcoch = coch2filtcoch_allsubbands(...
                        padded_coch, P_mod, complex_filters, fourier_domain, separable);
                    filtcoch = abs(filtcoch);
                case 'real'
                    complex_filters = false;
                    filtcoch = coch2filtcoch_allsubbands(...
                        padded_coch, P_mod, complex_filters, fourier_domain, separable);
                otherwise
                    error('No matching case');
            end
            clear complex_filters fourier_domain;

            % remove padding
            F = filtcoch((1:wav_dur_sec*P_mod.env_sr), 1:size(coch,2), :, :, :);
            clear filtcoch padded_coch;
            
            % resample
            if P_orig.sr ~= P_mod.env_sr
                F = resample_ndarray(F, P_orig.sr, P_mod.env_sr, 1);
            end
            
            % save
            P = P_mod;
            save(filtcoch_MAT_files{j}, 'F', 'P');
            clear F P;
            
        end
    end
    clear coch P_coch;
end

%% Principal components for cochleograms and filtered cochleograms

% feature file names and directory
n_model_features = length(modulation_types)+1;
model_features = cell(1, n_model_features);
feature_directories = cell(1, n_model_features);
pca_output_directories = cell(1, n_model_features);
model_features{1} = 'coch';
feature_directories{1} = coch_output_directory;
pca_output_directories{1} = [output_directory '/' coch_idstring '_' pca_idstring];
for j = 1:length(modulation_types)
    model_features{j+1} = [modulation_types{j} '_' P_orig.nonlin];
    feature_directories{j+1} = modulation_output_directory;
    pca_output_directories{j+1} = [output_directory '/' modulation_idstring '_' pca_idstring];
end

% initialize PCA matrices
pca_weights_allmodels = cell(1, n_model_features);
pca_eigvals_allmodels = cell(1, n_model_features);
if P_orig.return_PCA_timecourses
    pca_timecourses_allstim_allmodels = cell(length(stimuli), n_model_features);
else
    pca_timecourses_allstim_allmodels = {};
end
pca_weight_MAT_files = cell(1, n_model_features);
pca_timecourse_MAT_files = cell(length(stimuli), n_model_features);
for i = 1:n_model_features
    
    % PCA weights from covariance
    if ~exist(pca_output_directories{i}, 'dir'); mkdir(pca_output_directories{i}); end
    pca_weight_MAT_files{i} = [pca_output_directories{i} '/' model_features{i} '_PCA_weights'];
    if ~exist(pca_weight_MAT_files{i}, 'file') || P_orig.overwrite
        
        % un-normalized covariance matrix
        N_all = 0;
        for j = 1:length(P_orig.stim_to_compute_PC_weights)
            
            % load and format
            [~,fname,~] = fileparts(P_orig.stim_to_compute_PC_weights{j});
            fprintf('feat cov: %s, %s\n', model_features{i}, P_orig.stim_to_compute_PC_weights{j}); drawnow;
            load([feature_directories{i} '/' model_features{i} '_' fname '.mat'], 'F', 'P');
            dims = size(F);
            F = reshape(F, dims(1), prod(dims(2:end)));
            
            % initialize
            if j == 1
                C = zeros(size(F,2));
                M = zeros(1,size(F,2));
            end
            
            % accumulate correlation, sum, and total number of samples
            C = C + F' * F;
            M = M + sum(F);
            N_all = N_all + size(F,1);
            assert(all(isreal(C)))
            
        end
        
        % normalize by number of samples
        C = C/N_all;
        M = M/N_all;
        S = sqrt(diag(C)' - M.^2);
        assert(all(isreal(C)));
        
        % optionally remove effects of means and standard deviations from
        % correlation matrix
        if P_orig.demean_feats
            CV = C - M' * M;
        else
            CV = C;
        end
        if P_orig.std_feats
            CV = (1./(S' * S)) .* CV;
        else
            CV = CV;
        end
        assert(all(isreal(CV)));
        
        % eigenvectors of covariance
        [pca_weights, pca_eigvals] = eig(CV);
        
        % sort by variance
        pca_eigvals = diag(pca_eigvals);
        [~, xi] = sort(pca_eigvals, 'descend');
        pca_weights = pca_weights(:,xi);
        pca_eigvals = pca_eigvals(xi);
        
        save(pca_weight_MAT_files{i}, 'pca_weights', 'pca_eigvals', 'P', 'M', 'S')
        
    else
        
        load(pca_weight_MAT_files{i}, 'pca_weights', 'pca_eigvals', 'P', 'M', 'S')
        
    end
    pca_weights_allmodels{i} = pca_weights;
    pca_eigvals_allmodels{i} = pca_eigvals;
    
    % Apply weights to get timecourses
    fprintf('\n\n');
    for j = 1:length(stimuli)
        
        [~,fname,~] = fileparts(stimuli{j});
        pca_timecourse_MAT_files{j,i} = [pca_output_directories{i} '/' ...
            model_features{i} '_PCA_timecourses_' num2str(min(size(pca_weights,2), P_orig.nPCs)) '_' fname '.mat'];
        if ~exist(pca_timecourse_MAT_files{j,i}, 'file') || P_orig.overwrite
            
            % load and format
            fprintf('pca timecourses: %s, %s\n', model_features{i}, stimuli{j}); drawnow;
            load([feature_directories{i} '/' model_features{i} '_' fname '.mat'], 'F', 'P');
            dims = size(F);
            F = reshape(F, dims(1), prod(dims(2:end)));
            
            % optionally demean and standardize
            if P_orig.demean_feats
                F = bsxfun(@minus, F, M);
            end
            if P_orig.std_feats
                F = bsxfun(@times, F, 1./S);
            end
            
            % apply weights
            pca_timecourses = F * pca_weights(:,1:min(size(pca_weights,2), P_orig.nPCs));
            clear F;

            % save results
            save(pca_timecourse_MAT_files{j,i}, 'pca_timecourses', 'P', 'dims');
            
        elseif P_orig.return_PCA_timecourses
            
            load(pca_timecourse_MAT_files{j,i}, 'pca_timecourses', 'P', 'dims');
            
        else % do nothing
    
        end
        
        % save
        if P_orig.return_PCA_timecourses
            pca_timecourses_allstim_allmodels{j,i} = pca_timecourses;
        end
        clear pca_timecourses;
        
    end
    fprintf('\n\n');
end