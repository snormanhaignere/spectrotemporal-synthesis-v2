function [wav, P, R] = coch2wav_without_filts(coch, P, R)

% same as wav2coch.m but doesn't require the filters to be specified, just a
% parameter struct
% 
% 2016-06-23: Created by Sam NH
% 
% 2018-03-18: Made it possible to use ferret ERBs by adding animal field

if ~isfield(P, 'animal')
    P.animal = 'human';
end

duration_sec = size(coch,1) / P.env_sr;

% cochleogram filters
if P.overcomplete==0
    [audio_filts, audio_low_cutoff] = ...
        make_erb_cos_filters(duration_sec*P.audio_sr, P.audio_sr, ...
        P.n_filts, P.lo_freq_hz, P.audio_sr/2, P.animal);
    
elseif P.overcomplete==1
    [audio_filts, audio_low_cutoff] = ...
        make_erb_cos_filts_double2(duration_sec*P.audio_sr, P.audio_sr, ...
        P.n_filts, P.lo_freq_hz, P.audio_sr/2, P.animal);
    
elseif P.overcomplete==2
    [audio_filts, audio_low_cutoff] = ...
        make_erb_cos_filts_quadruple2(duration_sec*P.audio_sr, P.audio_sr, ...
        P.n_filts, P.lo_freq_hz, P.audio_sr/2, P.animal);
end

% remove filters below and above desired cutoffs
xi = audio_low_cutoff > P.lo_freq_hz - 1e-3 ...
    & audio_low_cutoff < P.audio_sr/2 + 1e-3;
audio_filts = audio_filts(:,xi);
audio_low_cutoff = audio_low_cutoff(xi);

% convert to cochleogram
wav = coch2wav(coch, R, ...
    audio_filts, audio_low_cutoff, P.audio_sr, P.env_sr);
