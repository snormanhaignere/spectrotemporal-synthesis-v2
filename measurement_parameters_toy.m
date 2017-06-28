function P = measurement_parameters_toy

% Toy parameters used to test the code
% 
% 2017-06-19: Created, Sam NH

% temporal modulation rates of second layer
% and whether the filters are lowpass or bandpass (the default)
P.temp_mod_rates = [2 2 8];
P.temp_mod_lowpass = [1 0 0];

% spectral modulation rates of second layer
% and whether the filters are lowpass or bandpass (the default)
P.spec_mod_rates = [0.5 0.5 2];
P.spec_mod_lowpass = [1 0 0];

% number of PCs to compute in the third layer
P.n_third_layer_PCs = 50;

% maximum duration of the input and synthesis sound in seconds
P.max_duration_sec = 12;

% amount of frequency padding, twice the period of the lowest spectral scale
P.freq_pad_oct = 1;

% temporal padding
P.temp_pad_sec = 1;

% audio sampling rate
P.audio_sr = 20000;

% sampling rate of the envelope in seconds
P.env_sr = 50;

% lowest filter in the audio filter bank
% highest is the nyquist - P.audio_sr/2
P.lo_freq_hz = 100;

% number of cosine filters to use
% increasing the number of filters
% decreases the bandwidth of the filters
P.n_filts = round((freq2erb(P.audio_sr)-freq2erb(P.lo_freq_hz))/1.3581);

% whether or not the number of filters is
% complete (=0), 1x overcomplete (=1), or 2x overcomplete (=2)
% overcomplete representations typically result in slightly more compelling 
% synthetics, but require more time and memory 
P.overcomplete = 0;

% frequency spacing of samples after interpolation to a logarithmic scale
% in octaves
P.logf_spacing = 1/6;

% factor to which cochleogram envelopes are raised
P.compression_factor = 0.3;

% whether or note the filters are causal in time
P.causal = true;