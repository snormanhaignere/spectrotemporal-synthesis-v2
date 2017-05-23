function run_spectrotemporal_synthesis(...
    P, fname, input_directory, output_directory, varargin)

% Simplified version of the algorithm that does not rely directly on shihab's
% code. There is only a single subband per filter (not real and imaginary). And
% the subbands are not all returned at once in a large matrix, substantially
% lowering the memory requirements and run time. 
% 
% run_spectrotemporal_synthesis_v4(...
%     P, fname, input_directory, output_directory)
% 
% Top-level script for running synthesis
% 
% -- Input --
% 
% P: parameter structure, see default_synthesis_parameters.m
% 
% fname: name of the audio file for the input sound file to try and match
% 
% input_directory: directory containing the input file
% 
% output_directory: directory to save results to

%% Directories, Paths, Misc

% find the directory containing this file
name_of_this_file = mfilename;
directory_containing_this_file = fileparts(which(name_of_this_file));

% McDermott Texture toolbox
addpath(genpath([directory_containing_this_file ...
    '/Sound_Texture_Synthesis_Toolbox']));

% 2D FT tools
addpath(genpath([directory_containing_this_file '/plot_2DFT']));

% reset the random stream
% allows identical synthetics to be created
% if the random seed is the same
ResetRandStream2(P.random_seed);

%% Read and format input sound

% separate file name and extension
fname_split = regexp(fname, '\.', 'split');
if length(fname_split) == 1    
    error('Input file name lacks an extension');
else
    fname_without_extension = cat(2,fname_split{1:end-1});
    fname_extension = fname_split{end};
end
clear fname_split fname;

% read wave file
[wav_orig, wav_sr] = audioread(...
    [input_directory '/' fname_without_extension '.'  fname_extension]);
wav_orig = format_wav(wav_orig, wav_sr, P);

% write formatted waveform to file
audiowrite_checkclipping(...
    [output_directory '/' fname_without_extension '_orig.wav'], ...
    0.01*wav_orig / rms(wav_orig), P.audio_sr);

%% Creates synthetic or load synthetic from previous run of the algorithm

% add the random seed to the filename
fname_without_extension = [fname_without_extension '_' num2str(P.random_seed)];

% mat file that stores the synthetic and other useful infor
synth_mat_file = [output_directory '/' fname_without_extension '_synth.mat'];

if exist(synth_mat_file, 'file');
    
    load(synth_mat_file, ...
        'wav_synth', 'n_completed_iterations', 'C', 'M_synth', 'M_orig');
    starting_iteration = n_completed_iterations + 1; %#ok<NODEF>
    
else
    
    wav_synth = format_wav(randn(size(wav_orig)), P.audio_sr, P);    
    starting_iteration = 1;
    
end

% save the parameters
save([output_directory '/' fname_without_extension '_parameters.mat'], 'P');

%% Cochleograms

duration_sec = length(wav_orig) / P.audio_sr;

% cochleogram filters
if P.overcomplete==0
    [audio_filts, audio_low_cutoff] = ...
        make_erb_cos_filters(duration_sec*P.audio_sr, P.audio_sr, ...
        P.n_filts, P.lo_freq_hz, P.audio_sr/2);
    
elseif P.overcomplete==1
    [audio_filts, audio_low_cutoff] = ...
        make_erb_cos_filts_double2(duration_sec*P.audio_sr, P.audio_sr, ...
        P.n_filts, P.lo_freq_hz, P.audio_sr/2);
    
elseif P.overcomplete==2
    [audio_filts, audio_low_cutoff] = ...
        make_erb_cos_filts_quadruple2(duration_sec*P.audio_sr, P.audio_sr, ...
        P.n_filts, P.lo_freq_hz, P.audio_sr/2);
end

% remove filters below and above desired cutoffs
xi = audio_low_cutoff > P.lo_freq_hz - 1e-3 ...
    & audio_low_cutoff < P.audio_sr/2 + 1e-3;
audio_filts = audio_filts(:,xi);
audio_low_cutoff = audio_low_cutoff(xi);

fprintf('Measuring cochleograms...\n');
drawnow;

% cochleogram of original sound
[coch_orig, P.f, P.t, R_orig] = ...
    wav2coch(wav_orig, audio_filts, audio_low_cutoff, ...
    P.audio_sr, P.env_sr, P.compression_factor, P.logf_spacing);

% cochleogram of synthetic sound
[coch_synth, ~, ~, R_synth] = ...
    wav2coch(wav_synth, audio_filts, audio_low_cutoff, ...
    P.audio_sr, P.env_sr, P.compression_factor, P.logf_spacing);

% temporal indices to match, excluding the first P.buffer_sec seconds
n_t = size(coch_orig, 1);
% ti_to_match = {1:n_t};
ti_to_match = 1 : n_t;

%% Spectrotemporal features

% whether or not there are any modulation stats in the model to match
match_mod_stats = ...
    any( P.match_temp_mod || P.match_spec_mod || P.match_spectemp_mod );

if match_mod_stats
    
    % determines the filters to match 
    % given the specified parameters;
    P = determine_filters_to_match(P);
    
end

% figure;
% plot(P.temp_mod_to_match', P.spec_mod_to_match','o');
% drawnow;
% WaitSecs(1);

%% Run synthesis loop

% compute moments of the filteres for the original waveform
fprintf('Computing moments of filter responses for original waveform...\n');
M_orig = all_filter_moments_from_coch(coch_orig, P, ti_to_match);

for i = starting_iteration:P.n_iter+1
    
    fprintf('\n\nSynthesis iteration %d\n\n', i);
    
    % computes moments of the envelope of cochlear and modulation filters
    fprintf('Computing moments of filter responses for synthetic...\n');
    M_synth = all_filter_moments_from_coch(coch_synth, P, ti_to_match);
    
    % comparison of the moments for this iteration
    if ~exist('C', 'var');
        C = struct;
    end
    C = moment_comparisons(M_orig, M_synth, i, C);

    % break if the desired number of iterations have been run
    % this is useful because moments are measured n_iter+1 times
    if i > P.n_iter
        break;
    end
        
    % match modulation filter stats
    if match_mod_stats
        
        % match the cochlear envelopes before first matching the spectrotemporal
        % envelopes
        if i == 1 && P.match_coch
            % match cochlear statistics before first measuring filtered cochleograms
            if ~isempty(P.match_coch_every_Nsec)
                coch_synth = ...
                    match_coch_hists_every_Nsec(...
                    coch_orig, coch_synth, ti_to_match, P);
            else
                if P.match_just_mean_and_var
                    coch_synth = match_coch_mean_and_var(...
                        coch_orig, coch_synth, ti_to_match);
                else
                    coch_synth = match_coch_hists(...
                        coch_orig, coch_synth, ti_to_match);
                end
            end
        end
        
        % match histograms of spectrotemporal filters
        fprintf('Matching modulation filter envelopes...\n');
        ti_to_match_assuming_padding = ...
            ti_to_match + round(P.env_sr * P.temp_pad_sec);
        coch_synth = ...
            match_filtcoch(...
            pad_coch(coch_orig, P), pad_coch(coch_synth, P), ...
            P, ti_to_match_assuming_padding, P.match_just_mean_and_var);

        % remove frequency padding
        coch_synth = remove_pad(coch_synth, P);
        
    end
    
    % match cochlear histograms
    if P.match_coch
        fprintf('Matching cochlear envelopes...\n');
        if ~isempty(P.match_coch_every_Nsec)
            coch_synth = ...
                match_coch_hists_every_Nsec(...
                coch_orig, coch_synth, ti_to_match, P);
        else
            if P.match_just_mean_and_var
                coch_synth = match_coch_mean_and_var(...
                    coch_orig, coch_synth, ti_to_match);
            else
                coch_synth = match_coch_hists(...
                    coch_orig, coch_synth, ti_to_match);
            end
        end
    end
        
    % reconstruct waveform
    fprintf('Reconstructing waveform...\n');
    wav_synth = coch2wav(coch_synth, R_synth, ...
        audio_filts, audio_low_cutoff, P.audio_sr, P.env_sr);
    
    % remeasure cochleogram from reconstructed waveform
    fprintf('Recomputing cochleogram...\n');
    [coch_synth, ~, ~, R_synth] = ...
        wav2coch(wav_synth, audio_filts, audio_low_cutoff, ...
        P.audio_sr, P.env_sr, P.compression_factor, P.logf_spacing);
    
    % save information
    n_completed_iterations = i; %#ok<NASGU>
    save(synth_mat_file, 'wav_synth', 'n_completed_iterations',...
        'C', 'M_synth', 'M_orig');
    audiowrite_checkclipping(...
        strrep(synth_mat_file, '.mat', '.wav'), ...
        0.01*wav_synth / rms(wav_synth), P.audio_sr);
    
    close all;
    
    % plot cochleograms
    plot_cochleogram_orig_and_synth(coch_orig, coch_synth, ...
        P, output_directory, fname_without_extension)
    
    % summary of moment comparisons
    if i > 1
        plot_moment_comparison_summary(C, ...
            output_directory, fname_without_extension);
    end
        
end

%% Calculate histogram divergences

D = histogram_similarity(coch_orig, coch_synth, P); %#ok<NASGU>
save(synth_mat_file, 'D', '-append');

%% Plots

% plot cochleograms
plot_cochleogram_orig_and_synth(coch_orig, coch_synth, ...
    P, output_directory, fname_without_extension);
close all;

% summary of moment comparisons
plot_moment_comparison_summary(C, output_directory, fname_without_extension);
close all;

% detailed plots of the moments
plot_moments(M_orig, M_synth, P, output_directory, fname_without_extension);
close all;

% detailed plots of the moments
plot_moment_scatter(M_orig, M_synth, P, output_directory, fname_without_extension);
close all;

% filtered cochleograms
plot_filtered_cochleograms(coch_orig, coch_synth, ...
    P, output_directory, fname_without_extension);
close all;

% histograms of example cochlear filters
plot_coch_hists(coch_orig, coch_synth, ti_to_match, ...
    P, output_directory, fname_without_extension);
close all;

% histograms of example spectrotemporal filters
plot_spectrotemporal_hists(coch_orig, coch_synth, ti_to_match, ...
    P, output_directory, fname_without_extension);
close all;

% plot 2DFT of the original and synthetic cochleogram
plot_2DFT_orig_and_synth(...
    coch_orig, coch_synth, P, output_directory, fname_without_extension);
close all;

% F = 217;
% T = 500;
% Hts = nan(T,F,length(P.temp_mod_to_match));
% for i = 1:length(P.temp_mod_to_match)
%     fprintf('%d\n',i);
%     Hts(:,:,i) = filt_spectemp_mod(P.spec_mod_to_match(i), P.temp_mod_to_match(i), F, T, P);
% end

% %%
% D = nan(length(P.temp_mod_to_match), length(P.temp_mod_to_match));
% for i = 1:length(P.temp_mod_to_match)
%     fprintf('%d\n',i);
%     for j = i+1:length(P.temp_mod_to_match)
%         D(i,j) = SumAll(abs(Hts(:,:,i)-Hts(:,:,j)).^2);
%     end
% end