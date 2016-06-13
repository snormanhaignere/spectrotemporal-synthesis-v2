function [S] = edit_measured_stats(S, P)
% [S] = EDIT_MEASURED_STATS(S, P)
%
% function to modify measured statistics S (typically the target statistics
% in the synthesis procedure) in order to deal with potential instability
% in bands with low envelope variance.
% also, sets some of them to values for noise, or to be of opposite sign,
% as specified in the parameters P
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

%first, measure the statistics in a pink noise sample
N_P = P; %use same auditory model parameters, just with longer signal
N_P.orig_sound_filename = 'pink_noise_20s_20kHz.wav'; %must be in path for this to work
N_P.orig_sound_folder = '';
N_P.max_orig_dur_s = 20; %use full 20 seconds to get maximally reliable measurements
pink_noise = format_orig_sound(N_P);
measurement_win = set_measurement_window(length(pink_noise),N_P.measurement_windowing,N_P);
noise_S = measure_texture_stats(pink_noise,N_P,measurement_win);

num_channels = length(S.subband_var);

% prevent skew and kurtosis from blowing up when variance is small, by
% setting them to value for a noise signal
var_threshold =.02;
S.env_skew(S.env_var<var_threshold) = noise_S.env_skew(S.env_var<var_threshold);
S.env_kurt(S.env_var<var_threshold) = noise_S.env_kurt(S.env_var<var_threshold);

%Because modulation power statistic is normalized by the total variance in
%the envelope, when the variance is small the mod. power is set to the
%value for noise (as otherwise it can be unstable)
for j=1:num_channels
    if S.env_var(j)<var_threshold
        S.mod_power(j,:) = noise_S.mod_power(j,:);
    end
end

%Negate skew and/or C2 as specified
for j=1:num_channels
    if j>9
        chan_to_neg=[1:6]; %this assumes 6 M2 channels, deals with fact that 
    elseif j>2              % real part of M2 has to be negative in regions
        chan_to_neg=[1:5]; % where envelopes have no power
    else
        chan_to_neg=[1:4];
    end
    
    if P.neg_mod_C2==1
        S.mod_C2(j,:, 2) = -S.mod_C2(j,:, 2);
    elseif P.neg_mod_C2==2
        S.mod_C2(j,chan_to_neg, 1) = -S.mod_C2(j,chan_to_neg, 1);
    elseif P.neg_mod_C2==3
        S.mod_C2(j,:, 2) = -S.mod_C2(j,:, 2);
        S.mod_C2(j,chan_to_neg, 1) = -S.mod_C2(j,chan_to_neg, 1);
    end
end

if P.neg_env_skew
    S.env_skew = -S.env_skew;
end

%Set stats to values for noise as specified
if P.use_noise_stats.sub_var
    S.subband_var = noise_S.subband_var;
end
if P.use_noise_stats.sub_kurt
    S.subband_kurt = noise_S.subband_kurt;
end
if P.use_noise_stats.env_mean
    S.env_mean = noise_S.env_mean;
end
if P.use_noise_stats.env_var
    S.env_var = noise_S.env_var;
end
if P.use_noise_stats.env_skew
    S.env_skew = noise_S.env_skew;
end
if P.use_noise_stats.env_kurt
    S.env_kurt = noise_S.env_kurt;
end
if P.use_noise_stats.env_C
    S.env_C(:,:,end) = noise_S.env_C(:,:,end);
end
if P.use_noise_stats.env_ac
    S.env_ac = noise_S.env_ac;
end
if P.use_noise_stats.mod_pow
    S.mod_power = noise_S.mod_power;
end
if P.use_noise_stats.mod_C1
    S.mod_C1 = noise_S.mod_C1;
end
if P.use_noise_stats.mod_C2
    S.mod_C2 = noise_S.mod_C2;
end

