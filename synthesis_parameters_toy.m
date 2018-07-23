function P = synthesis_parameters_toy

% whether to just match the mean and variance of each filter response, or to
% instead match the full histogram
P.match_just_mean_and_var = false;

% whether or not to match cochlear marginals, temporal modulation filters, or
% spectrotemporal modulation filters
P.match_coch = 1; 
P.match_temp_mod = 1;
P.match_spec_mod = 1;
P.match_spectemp_mod = 1;

% number of iterations to run algorithm for
P.n_iter = 10;

% maximum duration of the input and synthesis sound in seconds
P.max_duration_sec = 2;

% center frequencies of the temporal modulation filters in Hz
% 0 indicates a filter with only power at the DC
P.temp_mod_rates = [0,2,8,32];

% center frequencies of the spectral modulation filters in cyc/octave
% 0 indicates a filter with only power at the DC
P.spec_mod_rates = [0,1,2];

% additional lowrate temporal filters these filters are used to encourage the
% stimuli to have an even distribution of energy across the duration of the
% stimulus the filters are modulated in time and can be broadband or narrowband
% in frequency
P.lowrate_tempfilts_flat_spec = [];
P.lowrate_tempfilts_impulse_spec = [0.25];

% temporal rates and spectral scales for which filtered cochleograms are plotted
P.temp_mod_to_plot = [2 32];
P.spec_mod_to_plot = [1 4];

% frequency padding
% 2x slowest period
P.freq_pad_oct = 2/min(P.spec_mod_rates(P.spec_mod_rates>0));

% temporal padding
% 3x the longest period
all_temp_rates = [P.temp_mod_rates(P.temp_mod_rates>0), ...
    P.lowrate_tempfilts_flat_spec, P.lowrate_tempfilts_impulse_spec];
P.temp_pad_sec = 3/min(all_temp_rates);

% maximum duration of the input and synthesis sound in seconds
P.max_duration_sec = 1;

% duration of buffer-zone used to avoid minimize temporal wrap-around effects
% P.buffer_sec = 0;

% audio sampling rate
P.audio_sr = 8000;

% sampling rate of the envelope in seconds
P.env_sr = 100;

% lowest filter in the audio filter bank
% highest is the nyquist - P.audio_sr/2
P.lo_freq_hz = 20;

% number of cosine filters to use
% increasing the number of filters
% decreases the bandwidth of the filters
P.n_filts = 20;

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

% can optionally match the cochleogram every N seconds
% if this variable is empty, then the cochleogram is matched
% over the entire duration of the clip
P.match_coch_every_Nsec = [];

% number to seed the random number generator
% fix this seed if you want to the synthetics for a given
% input to always be the same
P.random_seed = randi(2^32);

% bandwidths and types/shapes of the spectrotemporal filters
P.spec_BW = 1;
P.temp_BW = 1;
P.spec_wavelet = 'mexicanhat'; 
P.temp_wavelet = 'gammatone';

% can randomize the filters or their phase, not recommended
P.spec_random_phase = false;
P.temp_random_phase = false;
P.spec_random_filt = false;
P.temp_random_filt = false;