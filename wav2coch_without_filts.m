function [coch, P, R] = wav2coch_without_filts(wav, P)

% same as wav2coch.m but doesn't require the filters to be specified, just a
% parameter struct
% 
% 2016-06-23: Created by Sam NH
% 
% 2017-05-11: Very small fix, prior could request non-integer number of samples
% from the filter

% cochleogram filters
if P.overcomplete==0
    [audio_filts, audio_low_cutoff] = ...
        make_erb_cos_filters(length(wav), P.audio_sr, ...
        P.n_filts, P.lo_freq_hz, P.audio_sr/2);
    
elseif P.overcomplete==1
    [audio_filts, audio_low_cutoff] = ...
        make_erb_cos_filts_double2(length(wav), P.audio_sr, ...
        P.n_filts, P.lo_freq_hz, P.audio_sr/2);
    
elseif P.overcomplete==2
    [audio_filts, audio_low_cutoff] = ...
        make_erb_cos_filts_quadruple2(length(wav), P.audio_sr, ...
        P.n_filts, P.lo_freq_hz, P.audio_sr/2);
end

% remove filters below and above desired cutoffs
xi = audio_low_cutoff > P.lo_freq_hz - 1e-3 ...
    & audio_low_cutoff < P.audio_sr/2 + 1e-3;
audio_filts = audio_filts(:,xi);
audio_low_cutoff = audio_low_cutoff(xi);

% cochleogram 
[coch, P.f, P.t, R] = ...
    wav2coch(wav, audio_filts, audio_low_cutoff, ...
    P.audio_sr, P.env_sr, P.compression_factor, P.logf_spacing);