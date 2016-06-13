function padded_coch = pad_coch_freq(coch, P)

% Adds frequency padding to a cochleogram

% amount of padding to add in samples
n_freq_pad_smps = round(P.freq_pad_oct / P.logf_spacing);

% matrix set to the mean of the cochleogram
M = mean(coch(:)) * ones(size(coch,1), n_freq_pad_smps);

% append padding
padded_coch = [M, coch];
