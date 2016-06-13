function M = all_filter_moments_from_wav(wav, P, ti_center)

% Calculates moments of cochlear and modulation filters from an input
% waveform. See all_filter_moments_from_coch called by this function.

% cochleogram from cosine envelopes
% -> time x frequency matrix
coch = ...
    wav2coch(wav, P.audio_filts, P.audio_cfs_Hz, P.audio_sr, P.env_sr);

% moments from the cochleogram
M = all_filter_moments_from_coch(coch, P, ti_center);

