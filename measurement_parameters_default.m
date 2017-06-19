function P = measurement_parameters_default

% maximum duration of the input and synthesis sound in seconds
P.max_duration_sec = 12;

% center frequencies of the temporal modulation filters in Hz
% 0 indicates a filter with only power at the DC
% wether or not the rates are lowpass or bandpass (the default)
P.temp_mod_rates = [0.5,1,2,4,8,0.5,1,2,4,8,16,32,64,128];
P.temp_mod_lowpass = [ones(1,5), zeros(1,9)]; 

% center frequencies of the spectral modulation filters in cyc/octave
% 0 indicates a filter with only power at the DC
% wether or not the scales are lowpass or bandpass (the default)
P.spec_mod_rates = [0.25,0.5,1,0.25,0.5,1,2,4,8];
P.spec_mod_lowpass = [ones(1,3), zeros(1,6)];

% amount of frequency padding, twice the period of the lowest spectral scale
P.freq_pad_oct = 2/min(P.spec_mod_rates(P.spec_mod_rates>0));

% temporal padding
% 3x the longest period in the synthesis
all_temp_rates = [P.temp_mod_rates(P.temp_mod_rates>0), ...
    P.lowrate_tempfilts_flat_spec, P.lowrate_tempfilts_impulse_spec];
P.temp_pad_sec = 3/min(all_temp_rates);

% audio sampling rate
P.audio_sr = 20000;

% sampling rate of the envelope in seconds
P.env_sr = 400;

% lowest filter in the audio filter bank
% highest is the nyquist - P.audio_sr/2
P.lo_freq_hz = 50;

% number of cosine filters to use
% increasing the number of filters
% decreases the bandwidth of the filters
P.n_filts = round((freq2erb(P.audio_sr)-freq2erb(P.lo_freq_hz))/1.3581);

% whether or not the number of filters is
% complete (=0), 1x overcomplete (=1), or 2x overcomplete (=2)
% overcomplete representations typically result in slightly more compelling 
% synthetics, but require more time and memory 
P.overcomplete = 2;

% frequency spacing of samples after interpolation to a logarithmic scale
% in octaves
P.logf_spacing = 1/24;

% factor to which cochleogram envelopes are raised
P.compression_factor = 0.3;