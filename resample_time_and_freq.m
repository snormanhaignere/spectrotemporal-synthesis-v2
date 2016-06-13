function Ytf = resample_time_and_freq(X, t1_sr, t2_sr, f1, f2)

% Resamples time axis and interpolates frequency axis of cochleogram X.
% 
% t1_sr: input temporal sampling rate
% t2_sr: output temporal sampling rate
% 
% f1: input frequencies
% f2: output frequencies

% resample time axis
Yt = resample(X, t2_sr, t1_sr);

% interpolate frequency axis
n_t = size(Yt,1);
Ytf = nan(n_t, length(f2));
for i = 1:n_t
    Ytf(i,:) = interp1(log2(f1), Yt(i,:)', log2(f2), 'pchip', 'extrap');
end