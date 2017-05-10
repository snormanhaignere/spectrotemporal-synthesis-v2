function M = all_filter_moments_from_wav(wav, sr, P, ti, envelopes)

% Calculates moments of cochlear and modulation filters from an input
% waveform. See all_filter_moments_from_coch called by this function.
% 
% 2016-11-18: Automaticaly resamples waveform to desired rate

% average across second dimension if present
wav = mean(wav,2);

if sr ~= P.audio_sr
    wav = resample(wav, P.audio_sr, sr);
    clear sr;
end

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

% remove filters below and above desired cutoffs
xi = audio_low_cutoff > P.lo_freq_hz - 1e-3 ...
    & audio_low_cutoff < P.audio_sr/2 + 1e-3;
audio_filts = audio_filts(:,xi);
audio_low_cutoff = audio_low_cutoff(xi);

% cochleogram 
[coch, P.f, P.t] = ...
    wav2coch(wav, audio_filts, audio_low_cutoff, ...
    P.audio_sr, P.env_sr, P.compression_factor, P.logf_spacing);

% temporal indices to compute moments over
% by default all temporal indices are used
if nargin < 4 || isempty(ti)
    ti = 1:size(coch,1);
end

if nargin < 5
    envelopes = false;
end

% moments from the cochleogram
M = all_filter_moments_from_coch(coch, P, ti, envelopes);

