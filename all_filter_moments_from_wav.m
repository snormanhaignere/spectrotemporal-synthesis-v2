function M = all_filter_moments_from_wav(wav, P, ti)

% Calculates moments of cochlear and modulation filters from an input
% waveform. See all_filter_moments_from_coch called by this function.

duration_sec = length(wav) / P.audio_sr;

% cochleogram filters
if P.overcomplete==0
    [audio_filts, audio_low_cutoff] = ...
        make_erb_cos_filters(duration_sec*P.audio_sr, P.audio_sr, ...
        P.n_filts, P.lo_freq_hz, P.audio_sr/2);
    
elseif P.overcomplete==1
    [audio_filts, audio_low_cutoff] = ...
        make_erb_cos_filts_double2(duration_sec*P.audio_sr, P.audio_sr, ...
        P.n_filts, P.lo_freq_hz, P.audio_sr/2);
    
elseif P.overcomplete==2
    [audio_filts, audio_low_cutoff] = ...
        make_erb_cos_filts_quadruple2(duration_sec*P.audio_sr, P.audio_sr, ...
        P.n_filts, P.lo_freq_hz, P.audio_sr/2);
end

% cochleogram 
[coch, P.f, P.t] = ...
    wav2coch(wav, audio_filts, audio_low_cutoff, ...
    P.audio_sr, P.env_sr, P.compression_factor, P.logf_spacing);

% temporal indices to compute moments over
% by default all temporal indices are used
if nargin < 3
    ti = 1:size(coch,1);
end

% moments from the cochleogram
M = all_filter_moments_from_coch(coch, P, ti);

