function P = synthesis_parameters_default

% whether or not to match cochlear marginals, temporal modulation filters, or
% spectrotemporal modulation filters
P.match_coch = 1; 
P.match_temp_mod = 1;
P.match_spec_mod = 0;
P.match_spectemp_mod = 1;

% number of iterations to run algorithm for
P.n_iter = 100;

% center frequencies of the temporal modulation filters in Hz
% 0 indicates a filter with only power at the DC
P.temp_mod_rates = [0,0.125,0.25,0.5,1,2,4,8,16,32,64,128];

% center frequencies of the spectral modulation filters in cyc/octave
% 0 indicates a filter with only power at the DC
P.spec_mod_rates = [0,0.25,0.5,1,2,4,8];

% temporal rates and spectral scales for which filtered cochleograms are plotted
P.temp_mod_to_plot = [1 4 16];
P.spec_mod_to_plot = [0.25 1 4];

% amount of frequency padding, twice the period of the lowest spectral scale
P.freq_pad_oct = 2/min(P.spec_mod_rates(P.spec_mod_rates>0));

% maximum duration of the input and synthesis sound in seconds
P.max_duration_sec = 12;

% duration of buffer-zone used to avoid minimize temporal wrap-around effects
P.buffer_sec = 2;

% audio sampling rate
P.audio_sr = 20000;

% sampling rate of the envelope in seconds
P.env_sr = 400;

% lowest filter in the audio filter bank
% highest is the nyquist - P.audio_sr/2
P.lo_freq_hz = 20;

% number of cosine filters to use
% increasing the number of filters
% decreases the bandwidth of the filters
P.n_filts = 30;

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

% number to seed the random number generator
% fix this seed if you want to the synthetics for a given
% input to always be the same
P.random_seed = randi(2^32);