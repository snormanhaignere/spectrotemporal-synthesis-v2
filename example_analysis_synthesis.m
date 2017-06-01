% Simple code illustrating how to do analysis and synthesis using the filters
% from the spectrotemporal model.
% 
% 2017-05-10: Created, Sam NH

%% Setup

% Parameters
P = synthesis_parameters_default;

% read in waveform
[wav,sr] = audioread('speech.wav');
wav = wav(1:sr);

%% Cochleogram

% cochleogram
[coch, P, R] = wav2coch_without_filts(wav, P);

% plot cochleogram
figure;
plot_cochleogram(coch, P.f, P.t)

%% Single subband

% single subband
spec_mod_rate = 2;
temp_mod_rate = 2;

% one-liner
coch_subband = coch2filtcoch(coch, spec_mod_rate, temp_mod_rate, P);

% plot subband
figure;
plot_cochleogram(coch_subband, P.f, P.t);

%% Multiple subbands

% include negative and positive rates
P.temp_mod_rates = unique([-P.temp_mod_rates, P.temp_mod_rates]);

% compute subbands
% time x frequency x spectral modulation x temporal modulation
filtcoch_allsubbands = coch2filtcoch_allsubbands(coch, P);

%% Invert back to cochleogram

coch_recon = filtcoch2coch(filtcoch_allsubbands, P);

% plot subband
figure;
subplot(2,1,1);
plot_cochleogram(coch, P.f, P.t);
subplot(2,1,2);
plot_cochleogram(coch_recon, P.f, P.t);

%% Reconstruct waveform

wav_recon = coch2wav_without_filts(coch_recon, P, R);

figure;
subplot(2,1,1);
plot(wav);
subplot(2,1,2);
plot(wav_recon);

