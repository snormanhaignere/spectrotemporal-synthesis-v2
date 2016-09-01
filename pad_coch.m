function temp_and_freq_padded_coch = pad_coch(coch, P)

% Adds temporal and frequency padding to the cochleogram
% 
% -- Example -- 
% P = synthesis_parameters_default;
% coch = randn(P.env_sr*P.max_duration_sec, 9*round(1/P.logf_spacing));
% padded_coch = pad_coch(coch, P);
% figure;
% imagesc(coch);
% figure;
% imagesc(padded_coch)

% amount of frequency padding to add in samples
n_freq_pad_smps = round(P.freq_pad_oct / P.logf_spacing);

% matrix set to the mean of the cochleogram
F = mean(coch(:)) * ones(size(coch,1), n_freq_pad_smps);

% append frequency padding padding
freq_padded_coch = [F, coch];
clear F n_freq_pad_smps;

% amount of temporal padding to add in samples
n_temp_pad_smps = round(P.env_sr * P.temp_pad_sec);

% matrix set to the mean of the cochleogram
T = mean(coch(:)) * ones(n_temp_pad_smps, size(freq_padded_coch,2));

% append frequency padding padding
temp_and_freq_padded_coch = [T; freq_padded_coch];
clear T n_temp_pad_smps;



