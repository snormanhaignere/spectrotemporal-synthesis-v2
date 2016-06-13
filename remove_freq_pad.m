function depadded_coch = remove_freq_pad(coch, P)

% Remove frequency padding from cochleogram

% amount of padding to add in samples
n_freq_pad_smps = round(P.freq_pad_oct / P.logf_spacing);

% remove padding
depadded_coch = coch(:, n_freq_pad_smps+1:end);
