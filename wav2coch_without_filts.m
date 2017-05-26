function [coch, P, R] = wav2coch_without_filts(wav, P)

% same as wav2coch.m but doesn't require the filters to be specified, just a
% parameter struct
% 
% 2016-06-23: Created by Sam NH
% 
% 2017-05-11: Very small fix, prior could request non-integer number of samples
% from the filter
% 
% 2017-05-26: Now relies on cochfilts.m

% filters
[audio_filts, audio_low_cutoff] = cochfilts(length(wav), P);

% cochleogram from filters
[coch, P.f, P.t, R] = ...
    wav2coch(wav, audio_filts, audio_low_cutoff, ...
    P.audio_sr, P.env_sr, P.compression_factor, P.logf_spacing);