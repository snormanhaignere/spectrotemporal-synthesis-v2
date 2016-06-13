function [subbands,subband_envs] = generate_subbands_and_envs(sample_sound, audio_sr, env_sr, N_audio_channels,...
    low_audio_f, hi_audio_f, lin_or_log_filters, use_more_audio_filters,compression_option,comp_exponent,log_constant)
%
% [SUBBANDS,SUBBAND_ENVS] = GENERATE_SUBBANDS_AND_ENVS(SAMPLE_SOUND, AUDIO_SR, ENV_SR, N_AUDIO_CHANNELS,...
%    LOW_AUDIO_F, HI_AUDIO_F, LIN_OR_LOG_FILTERS, USE_MORE_AUDIO_FILTERS,COMPRESSION_OPTION,COMP_EXPONENT,LOG_CONSTANT)
%
% generates subbands and their envelopes for an audio signal
% generates filters within function
%
% if LIN_OR_LOG_FILTERS is <= 2, filters are equally spaced on an ERB scale
% if LIN_OR_LOG_FILTERS is > 2, filters are equally spaced on a linear frequency scale
%

% This code is part of an instantiation of a sound texture synthesis
% algorithm developed with Eero Simoncelli and described in this paper:
%
% McDermott, J.H., Simoncelli, E.P. (2011) Sound texture perception via
% statistics of the auditory periphery: Evidence from sound synthesis.
% Neuron, 71, 926-940. 
%
% Dec 2012 -- Josh McDermott <jhm@mit.edu>


%compute sample sound subbands for future use
if lin_or_log_filters==1 || lin_or_log_filters==2
    if use_more_audio_filters==0
        [audio_filts, Hz_cutoffs] = make_erb_cos_filters(length(sample_sound), audio_sr, N_audio_channels, low_audio_f, hi_audio_f);
    elseif use_more_audio_filters==1
        [audio_filts, Hz_cutoffs] = make_erb_cos_filts_double2(length(sample_sound), audio_sr, N_audio_channels, low_audio_f, hi_audio_f);
    elseif use_more_audio_filters==2
        [audio_filts, Hz_cutoffs] = make_erb_cos_filts_quadruple2(length(sample_sound), audio_sr, N_audio_channels, low_audio_f, hi_audio_f);
    end
elseif lin_or_log_filters==3 || lin_or_log_filters==4
    if use_more_audio_filters==0
        [audio_filts, Hz_cutoffs] = make_lin_cos_filters(length(sample_sound), audio_sr, N_audio_channels, low_audio_f, hi_audio_f);
    elseif use_more_audio_filters==1
        [audio_filts, Hz_cutoffs] = make_lin_cos_filts_double(length(sample_sound), audio_sr, N_audio_channels, low_audio_f, hi_audio_f);
    end
end
subbands = generate_subbands(sample_sound, audio_filts);
analytic_subbands = hilbert(subbands);
subband_envs = abs(hilbert(subbands));
subband_envs = resample(subband_envs,env_sr,audio_sr);
subband_envs(subband_envs<0)=0;
if compression_option==1 %power compression
    subband_envs = subband_envs.^comp_exponent;
elseif compression_option==2 %log compression
    subband_envs = log10(subband_envs+log_constant);
end
