function padded_coch = pad_coch_freq(coch, P, pad_value)

% Adds frequency padding to a cochleogram

if nargin < 3
    pad_value = mean(coch(:));
end

% amount of padding to add in samples
n_freq_pad_smps = round(P.freq_pad_oct / P.logf_spacing);

% matrix set to the mean of the cochleogram
M = pad_value * ones(size(coch,1), n_freq_pad_smps);

% append padding
padded_coch = [M, coch];
