function S = measure_texture_stats(sample_sound, P, varargin)

% function S = MEASURE_TEXTURE_STATS(SAMPLE_SOUND, P, MEASUREMENT_WIN, ...
%   AUDIO_FILTS, AUDIO_CUTOFFS_HZ, MOD_FILTS, ENV_AC_FILTS, MOD_C1_FILTS, ...
%   MOD_C2_FILTS, SUB_HIST_BINS, SUB_ENV_HIST_BINS)
%
% Function to measure texture statistics from a sound signal.
%
% MEASUREMENT_WIN specifies the weights on each sample of 
% the signal in the average
%
% subsequent arguments (filters etc.) are optional and are created by
% function if not provided or if provided as empty ([]) arguments
%
% other parameters are specified in P (and can be set using
% SYNTHESIS_PARAMETERS.M)
%

% This code is part of an instantiation of a sound texture synthesis
% algorithm developed with Eero Simoncelli and described in this paper:
%
% McDermott, J.H., Simoncelli, E.P. (2011) Sound texture perception via
% statistics of the auditory periphery: Evidence from sound synthesis.
% Neuron, 71, 926-940. 
%
% Dec 2012 -- Josh McDermott <jhm@mit.edu>
% May 2014 -- modified to fix incompatibility with overcomplete filter
% banks -- Josh McDermott


if nargin>4 && ~isempty(varargin{2}) && ~isempty(varargin{3})
    audio_filts = varargin{2};
    audio_cutoffs_Hz = varargin{3};
else
    %generate audio filters given specified parameters
    if P.lin_or_log_filters==1 || P.lin_or_log_filters==2
        if P.use_more_audio_filters==0
            [audio_filts, audio_cutoffs_Hz] = make_erb_cos_filters(length(sample_sound), P.audio_sr, P.N_audio_channels, P.low_audio_f, P.hi_audio_f);
        elseif P.use_more_audio_filters==1
            [audio_filts, audio_cutoffs_Hz] = make_erb_cos_filts_double2(length(sample_sound), P.audio_sr, P.N_audio_channels, P.low_audio_f, P.hi_audio_f);
        elseif P.use_more_audio_filters==2
            [audio_filts, audio_cutoffs_Hz] = make_erb_cos_filts_quadruple2(length(sample_sound), P.audio_sr, P.N_audio_channels, P.low_audio_f, P.hi_audio_f);
        end
    elseif P.lin_or_log_filters==3 || P.lin_or_log_filters==4
        if P.use_more_audio_filters==0
            [audio_filts, audio_cutoffs_Hz] = make_lin_cos_filters(length(sample_sound), P.audio_sr, P.N_audio_channels, P.low_audio_f, P.hi_audio_f);
        elseif P.use_more_audio_filters==1
            [audio_filts, audio_cutoffs_Hz] = make_lin_cos_filts_double(length(sample_sound), P.audio_sr, P.N_audio_channels, P.low_audio_f, P.hi_audio_f);
        end
    end
end

%generate subbands and subband envelopes from which statistics are measured
subbands = generate_subbands(sample_sound, audio_filts);
%analytic_subbands = hilbert(subbands);
subband_envs = abs(hilbert(subbands));
if P.compression_option==1 %power compression
    subband_envs = subband_envs.^P.comp_exponent;
elseif P.compression_option==2 %log compression
    subband_envs = log10(subband_envs+P.log_constant);
end
ds_factor=P.audio_sr/P.env_sr;
subband_envs = resample(subband_envs,1,ds_factor);
subband_envs(subband_envs<0)=0;

if P.use_zp==1
    mod_filt_length = length(subband_envs)*2;
elseif P.use_zp==0
    mod_filt_length = length(subband_envs);
end

if nargin>2 && ~isempty(varargin{1})
    measurement_win = varargin{1};
else
    measurement_win = ones(length(subband_envs),1);
end

%generate modulation filters
if nargin>5 && ~isempty(varargin{4})
    mod_filts = varargin{4};
else
    if P.lin_or_log_filters==1 || P.lin_or_log_filters==3
        [mod_filts,Hz_mod_cfreqs,mod_freqs] = make_constQ_cos_filters(mod_filt_length, P.env_sr, P.N_mod_channels, P.low_mod_f, P.hi_mod_f, P.mod_filt_Q_value);
    elseif P.lin_or_log_filters==2 || P.lin_or_log_filters==4
        [mod_filts,Hz_mod_cfreqs,mod_freqs] = make_lin_cos_constQ_cntrl_filters(mod_filt_length, P.env_sr, P.N_mod_channels, P.low_mod_f, P.hi_mod_f, P.mod_filt_Q_value);
    end
end
if nargin>6 && ~isempty(varargin{5})
    env_ac_filts = varargin{5};
else
    [env_ac_filts,ac_filt_cutoffs] = make_env_acc_filters2(mod_filt_length, P.env_sr, P.env_ac_intervals_smp);
end
if nargin>6 && ~isempty(varargin{6}) && ~isempty(varargin{7})
    mod_C1_filts = varargin{6};
    mod_C2_filts = varargin{7};
else
    [mod_C2_filts,C2_cfreqs,C2_freqs] = make_octave_cos_filters2(mod_filt_length, P.env_sr, P.low_mod_f_C12, P.hi_mod_f);
    if P.lin_or_log_filters==1 || P.lin_or_log_filters==3
        mod_C1_filts = mod_C2_filts(:,2:end);
        C1_cfreqs = C2_cfreqs(2:end);
    elseif P.lin_or_log_filters==2 || P.lin_or_log_filters==4
        [mod_C1_filts,C1_cfreqs,C1_freqs] = make_lin_cos_oct_cntrl_filters(mod_filt_length, P.env_sr, size(mod_C2_filts,2)-1, C2_cfreqs(2), P.hi_mod_f);
    end
end

%compute marginals - of subbands, subband envelopes, and modulation bands
S.subband_mean = mean(subbands); %this should be zero for a bandpass filter, so not very meaningful; included for completeness
S.subband_var = var(subbands); %this is a row vector of the var of each subband
for j=1:size(audio_filts,2) %go through subbands
    S.subband_skew(j) = skewness(subbands(:,j));
    S.subband_kurt(j) = kurtosis(subbands(:,j));
    S.env_mean(j) = stat_central_moment_win(subband_envs(:,j),1,measurement_win);
    S.env_var(j) = stat_central_moment_win(subband_envs(:,j),2,measurement_win,S.env_mean(j));
    S.env_skew(j) = stat_central_moment_win(subband_envs(:,j),3,measurement_win,S.env_mean(j));
    S.env_kurt(j) = stat_central_moment_win(subband_envs(:,j),4,measurement_win,S.env_mean(j));
    if nargin>7 && ~isempty(varargin{8})
        sub_hist_bins = varargin{8}(j,:);
        [temp,bins]=hist(subbands(:,j),sub_hist_bins);
    else
        [temp,bins]=hist(subbands(:,j),P.n_hist_bins);
    end
    S.subband_hist(j,1:P.n_hist_bins)=temp/sum(temp);
    S.subband_bins(j,1:P.n_hist_bins)=bins;
    if nargin>8 && ~isempty(varargin{9})
        sub_env_hist_bins = varargin{9}(j,:);
        [temp,bins]=hist(subband_envs(:,j),sub_env_hist_bins);
    else
        [temp,bins]=hist(subband_envs(:,j),P.n_hist_bins);
    end
    S.env_hist(j,1:P.n_hist_bins)=temp/sum(temp);
    S.env_bins(j,1:P.n_hist_bins)=bins;

    S.env_ac(j,1:length(P.env_ac_intervals_smp)) = stat_env_ac_scaled_win(subband_envs(:,j), env_ac_filts, P.env_ac_intervals_smp, P.use_zp, measurement_win);
    S.mod_power(j,1:P.N_mod_channels) = stat_mod_power_win(subband_envs(:,j), mod_filts, P.use_zp, measurement_win);
    S.mod_C2(j,1:size(mod_C2_filts,2)-1,1:2) = stat_mod_C2_win(subband_envs(:,j), mod_C2_filts, P.use_zp, measurement_win);
end

%compute subband envelope, modulation band correlations
S.env_C = stat_corr_filt_win_full(subband_envs, ones(size(mod_C1_filts(:,1))), P.use_zp, measurement_win);
S.mod_C1 = stat_corr_filt_win_full(subband_envs, mod_C1_filts, P.use_zp, measurement_win);

%subband autocorrelation
sub_ac_N_smp = round(P.num_sub_ac_period./audio_cutoffs_Hz*P.audio_sr);
sub_ac_N_smp(sub_ac_N_smp > P.num_sub_ac_period/20*P.audio_sr)=P.num_sub_ac_period/20*P.audio_sr; %do not store more than the value for a 20 Hz channel
temp = autocorr_mult_zp(subbands, P.sub_ac_win_choice, P.sub_ac_undo_win);
L2 = length(subbands);
C2 = L2/2+1;
for k=1:size(audio_filts,2)
    S.subband_ac{k} = temp(C2-sub_ac_N_smp(k):C2+sub_ac_N_smp(k),k);
    S.subband_ac_power(k) = sum(S.subband_ac{k}.^2); %used in SNR calculation
end

%amplitude histogram of sample
[amp_hist,amp_bins] = hist(sample_sound, P.n_hist_bins);
S.subband_hist(size(audio_filts,2)+1,1:P.n_hist_bins)=amp_hist;
S.subband_bins(size(audio_filts,2)+1,1:P.n_hist_bins)=amp_bins;

