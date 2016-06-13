function P = synthesis_parameters_highres2x

P = synthesis_parameters_default;

% center frequencies of the temporal modulation filters in Hz
% 0 indicates a filter with only power at the DC
P.temp_mod_rates = [0, 2.^( log2(0.5):0.5:log2(128) )];

% center frequencies of the spectral modulation filters in cyc/octave
% 0 indicates a filter with only power at the DC
P.spec_mod_rates = [0, 2.^( log2(0.25):0.5:log2(8) )];