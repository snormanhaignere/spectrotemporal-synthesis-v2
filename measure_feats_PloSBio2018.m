% Shows how to measure the features used to predict voxel responses in
% Norman-Haignere et al. (2018). PLoS Biology.

% stimulus name
stim = 'woman_speaking_original.wav';
% stim = 'woman_speaking_modelmatched.wav';

% read waveform
[wav, sr] = audioread(stim);

% take one channel, if not mono
wav = wav(:,1);

% resample if needed
if sr ~= P.audio_sr
    wav = resample(wav, P.audio_sr, sr);
end

% load parameter structure
load('parameters_PLoSBio2018.mat', 'P');

%% To speed things up and reduce memory:
% - take just the first second
% - reduce padding
% comment this section out to run the full analysis
% wav = wav(1:P.audio_sr);
% P.temp_pad_sec = 0;
% P.freq_pad_oct = 0;

%% Cochleogram (stage 1)

% Need the scripts from texture toolbox
addpath('Sound_Texture_Synthesis_Toolbox');

% cochleogram filters
duration_sec = length(wav)/P.audio_sr;
if P.overcomplete==0
    [audio_filts, audio_low_cutoff] = ...
        make_erb_cos_filters(duration_sec*P.audio_sr, P.audio_sr, ...
        P.n_filts, P.lo_freq_hz, P.audio_sr/2);
    
elseif P.overcomplete==1
    [audio_filts, audio_low_cutoff] = ...
        make_erb_cos_filts_double2(duration_sec*P.audio_sr, P.audio_sr, ...
        P.n_filts, P.lo_freq_hz, P.audio_sr/2);
    
elseif P.overcomplete==2
    [audio_filts, audio_low_cutoff] = ...
        make_erb_cos_filts_quadruple2(duration_sec*P.audio_sr, P.audio_sr, ...
        P.n_filts, P.lo_freq_hz, P.audio_sr/2);
end

% remove filters below and above desired cutoffs
xi = audio_low_cutoff > P.lo_freq_hz - 1e-3 ...
    & audio_low_cutoff < P.audio_sr/2 + 1e-3;
audio_filts = audio_filts(:,xi);
audio_low_cutoff = audio_low_cutoff(xi);

% cochleogram of original sound
[coch, P.f, P.t, R_orig] = ...
    wav2coch(wav, audio_filts, audio_low_cutoff, ...
    P.audio_sr, P.env_sr, P.compression_factor, P.logf_spacing);

% plot cochleogram
figure;
plot_cochleogram(coch, P.f, P.t);

%% Stats of cochleogram and filtered cochleograms (stage 2)

% computes the first four moments of the filter responses:
% (1) mean (2) variance (3) skew (4) kurtosis
M = all_filter_moments_from_coch(coch, P, 1:size(coch,1));

% pick out mean of cochlear, standard deviation of all other feats
F.coch_env = M.coch_env(:,1);
F.temp_mod = sqrt(M.temp_mod(:,:,2));
F.spec_mod = sqrt(M.spec_mod(:,:,2));
F.spectemp_mod = sqrt(M.spectemp_mod(:,:,:,2));

% split out negative and positive temporal rates
% corresponding to upward and downward modulated ripples
% for prediction negative and positive rates were averaged
dims = size(F.spectemp_mod);
F.spectemp_mod = reshape(F.spectemp_mod, [dims(1), dims(2)/2, 2, dims(3)]);

%% Plot results

% plot average envelopes (spectrum-like measure)
figure;
semilogx(P.f, F.coch_env);
xlim([50, 10e3]);
xlabel('Frequency (Hz)');
set(gca, 'FontSize', 20);

% plot temporal modulation
figure;
imagesc(F.temp_mod');
temp_mod_rates_without_DC = P.temp_mod_rates(P.temp_mod_rates>0);
freqs_to_plot = [100 400 1600 6400];
fticks = interp1(P.f, length(P.f):-1:1, freqs_to_plot);
set(gca, 'YTick', fliplr(fticks), 'YTickLabel', fliplr(freqs_to_plot)/1000);
set(gca, 'XTick', [2,4,6,8], 'XTickLabel', round(temp_mod_rates_without_DC([2,4,6,8])))
set(gca, 'FontSize', 20);
ylabel('Audio frequency (kHz)');
xlabel('Rate (Hz)')
title('Temporal modulation');

% plot spectral modulation
figure;
imagesc(F.spec_mod');
freqs_to_plot = [100 400 1600 6400];
fticks = interp1(P.f, length(P.f):-1:1, freqs_to_plot);
set(gca, 'YTick', fliplr(fticks), 'YTickLabel', fliplr(freqs_to_plot)/1000);
set(gca, 'XTick', [2,4,6], 'XTickLabel', P.spec_mod_rates([2,4,6]))
set(gca, 'FontSize', 20);
ylabel('Audio frequency (kHz)');
xlabel('Scale (cyc/oct)');
title('Spectral modulation');

% plot spectrotemporal modulation
% for a given audio frequency
figure;
audiofreq = 200;
[~,xi] = min(abs(P.f-audiofreq));
X = cat(2, fliplr(F.spectemp_mod(:,:,2,xi)), F.spectemp_mod(:,:,1,xi));
imagesc(flipud(X));
spec_mod_rates_flip = fliplr(P.spec_mod_rates);
temp_mod_rates_neg_pos = [-fliplr(temp_mod_rates_without_DC), temp_mod_rates_without_DC];
set(gca, 'YTick', [1, 3, 5], 'YTickLabel', spec_mod_rates_flip([1 3 5]));
set(gca, 'XTick', [3, 7, 12, 16], 'XTickLabel', temp_mod_rates_neg_pos([3, 7, 12, 16]))
set(gca, 'FontSize', 20);
ylabel('Spectral scale (cyc/oct)');
xlabel('Temporal rate (Hz)');
title('Spectrotemporal modulation (200 Hz)');