function wav = coch2wav(coch_envs_dwnsmp_sr_logf, R, ...
    audio_filts, audio_low_cutoff, audio_sr, dwnsmp_sr)

% wav = coch2wav(coch_envs_dwnsmp_sr_logf, R, ...
%     audio_filts, audio_low_cutoff, audio_sr, dwnsmp_sr)
% 
% Reconstructs a waveform from collection of cochlear envelopes.
% 
% -- Inputs --
% 
% coch_envs_dwnsmp_sr_logf: [time x frequency] matrix of downsampled cochlear
% envelopes, transformed to a logarithmic frequency scale (see wav2coch)
% 
% R: structure containing the subband phases, logarithmic frequencies,
% compression factor, and envelope errors caused by downsampling/interpolation
% (see wav2coch)
% 
% audio_filts: cosine filter parameters (see make_erb_cos_filters)
% 
% audio_low_cutoff: lower cutoff of the cosine filters (see make_erb_cos_filters)
% 
% audio_sr: sampling rate of the audio waveform
% 
% dwnsmp_sr: downsampled sampling rate of the envelopes
% 
% -- Example --
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
% % reconstruction from cochleogram
% wav_reconstruct = coch2wav(coch, R, ...
%     audio_filts, audio_low_cutoff, P.audio_sr, P.env_sr);
% 
% % plot waveform
% subplot(2,1,1);
% plot(wav);
% subplot(2,1,2);
% plot(wav_reconstruct);

% resample back to audio sampling rate and erb scale
coch_envs_audio_sr_erbf_reconstructed = resample_time_and_freq(...
    coch_envs_dwnsmp_sr_logf, dwnsmp_sr, audio_sr, R.logf, audio_low_cutoff);

% truncate negatives
coch_envs_audio_sr_erbf_reconstructed(coch_envs_audio_sr_erbf_reconstructed<0)=0;

% undo compression
coch_envs_audio_sr_erbf_reconstructed = ...
    coch_envs_audio_sr_erbf_reconstructed .^ (1/R.compression_factor);

% remove reconstruction error caused by the above steps (see wav2coch)
coch_envs_audio_sr_erbf_reconstructed = ...
    coch_envs_audio_sr_erbf_reconstructed - R.recon_error_audio_sr_erbf;

% combine envelopes with phases
coch_subbands_audio_sr_erbf = ...
    real(coch_envs_audio_sr_erbf_reconstructed .* R.coch_phases_audio_sr_erbf);

% collapse subbands
wav = collapse_subbands(coch_subbands_audio_sr_erbf, audio_filts);