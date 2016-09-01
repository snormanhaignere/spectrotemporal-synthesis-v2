function depadded_coch = remove_pad(coch, P)

% Remove frequency padding from cochleogram
% 
% -- Example --
% P = synthesis_parameters_default;
% coch = randn(P.env_sr*P.max_duration_sec, 9*round(1/P.logf_spacing));
% padded_and_depadded_coch = remove_pad(pad_coch(coch, P), P);
% figure;
% plot(coch(:), padded_and_depadded_coch(:));

% amount of padding to add in samples
n_freq_pad_smps = round(P.freq_pad_oct / P.logf_spacing);

% amount of temporal padding to add in samples
n_temp_pad_smps = round(P.env_sr * P.temp_pad_sec);

% remove padding
depadded_coch = coch(n_temp_pad_smps+1:end, n_freq_pad_smps+1:end);
