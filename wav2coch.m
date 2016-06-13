function [coch_envs_compressed_dwnsmp_sr_logf, logf, t, R] = ...
    wav2coch(wav, audio_filts, audio_low_cutoff, ...
    audio_sr, dwnsmp_sr, compression_factor, logf_spacing)

% [coch_envs_dwnsmp_sr_logf, logf, t, R] = ...
%     wav2coch(wav, audio_filts, audio_low_cutoff, ...
%     audio_sr, dwnsmp_sr, compression_factor)
% 
% Calculates cochleogram envelopes from an audio waveform. The envelopes are
% downsampled (to dwnsmp_sr) and the frequency axis is interpolated to a
% logarithmic scale. A structure R is returned that along with the downsampled
% envelopes can be used to reconstruct the waveform using the function coch2wav.
% 
% -- Inputs --
% 
% wav: waveform from which to compute cochleogram
% 
% audio_filts: cosine filter parameters (see make_erb_cos_filters)
% 
% audio_low_cutoff: lower cutoff of the cosine filters (see make_erb_cos_filters)
% 
% audio_sr: sampling rate of the audio waveform
% 
% dwnsmp_sr: downsampled sampling rate of the envelopes
% 
% compression_factor: power to which subbands are raised (e.g. 0.3)
% 
% -- Outputs --
% 
% coch_envs_dwnsmp_sr_logf: cochleogram envelopes, downsampled in time and
% interpolated to a log frequency scale
% 
% logf: the logarithmic frequency scale used
% 
% t: vector of time indices for the subbands 
% 
% R: structure containing the subband phases, logarithmic frequencies,
% compression factor, and envelope errors caused by downsampling/interpolation;
% used to reconstruct the waveform (see coch2wav)
% 
% % -- Example --
% 
% % texture toolbox
% addpath(genpath([pwd '/Sound_Texture_Synthesis_Toolbox']));
% 
% % read in waveform
% [wav,sr] = audioread([pwd '/speech.wav']);
% P.max_duration_sec = 1;
% wav = 0.1 * format_wav(wav, sr, P);
% 
% % filters
% P = default_synthesis_parameters;
% [audio_filts, audio_low_cutoff] = ...
%     make_erb_cos_filters(length(wav), P.audio_sr, ...
%     P.n_filts, P.lo_freq_hz, P.audio_sr/2);
% 
% % cochleogram
% [coch, P.f, P.t, R] = ...
%     wav2coch(wav, audio_filts, audio_low_cutoff, ...
%     P.audio_sr, P.env_sr, P.compression_factor, P.logf_spacing);
% 
% % plot the cochleogram
% plot_cochleogram(coch, P.f, P.t);

% -- Analysis --

% frequency axis used for cochleograms
logf = 2.^(...
    log2(audio_low_cutoff(1)) : ...
    logf_spacing : ...
    log2(audio_low_cutoff(end)) + logf_spacing);

% save frequency axis for reconstruction
R.logf = logf;

% subbands from waveform
coch_subbands_audio_sr_erbf = generate_subbands(wav, audio_filts);

% subband envelopes
analytic_subbands_audio_sr_erbf = hilbert(coch_subbands_audio_sr_erbf);
coch_envs_audio_sr_erbf = abs(analytic_subbands_audio_sr_erbf);

% subband phases
R.coch_phases_audio_sr_erbf = analytic_subbands_audio_sr_erbf./coch_envs_audio_sr_erbf;
clear coch_subbands_audio_sr_erbf analytic_subbands_audio_sr_erbf;

% compress envelopes
R.compression_factor = compression_factor;
coch_envs_compressed_audio_sr_erbf = coch_envs_audio_sr_erbf .^ R.compression_factor;

% resample time and frequency axis
coch_envs_compressed_dwnsmp_sr_logf = resample_time_and_freq(...
    coch_envs_compressed_audio_sr_erbf, ...
    audio_sr, dwnsmp_sr, audio_low_cutoff, logf);

% time vector
n_t = size(coch_envs_compressed_dwnsmp_sr_logf,1);
t = (0:n_t-1)/dwnsmp_sr;

% truncate negatives
coch_envs_compressed_dwnsmp_sr_logf(coch_envs_compressed_dwnsmp_sr_logf<0)=0;

% -- Reverse analysis to compute errors for reconstruction --

% resample time and frequency axis
coch_envs_audio_sr_erbf_reconstructed = ...
    resample_time_and_freq(coch_envs_compressed_dwnsmp_sr_logf, ...
    dwnsmp_sr, audio_sr, logf, audio_low_cutoff);

% undo compression
coch_envs_audio_sr_erbf_reconstructed = ...
    coch_envs_audio_sr_erbf_reconstructed .^ (1/R.compression_factor);

% truncate
coch_envs_audio_sr_erbf_reconstructed(coch_envs_audio_sr_erbf_reconstructed<0)=0;

% compute error
R.recon_error_audio_sr_erbf = ...
    coch_envs_audio_sr_erbf_reconstructed - coch_envs_audio_sr_erbf;