function [wav, P, R] = coch2wav_without_filts(coch, P, R)

% same as wav2coch.m but doesn't require the filters to be specified, just a
% parameter struct
% 
% 2016-06-23: Created by Sam NH
% 
% 2017-05-26: Now relies on cochfilts

% filters
[audio_filts, audio_low_cutoff] = ...
    cochfilts(P.audio_sr * size(coch,1) / P.env_sr, P);

% convert to cochleogram
wav = coch2wav(coch, R, ...
    audio_filts, audio_low_cutoff, P.audio_sr, P.env_sr);
