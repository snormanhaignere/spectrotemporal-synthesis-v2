% Simple code illustrating how to do analysis and synthesis using the filters
% from the spectrotemporal model.
% 
% 2017-05-10: Created, Sam NH

%% Setup

% Parameters
P = measurement_parameters_default;

% read in waveform
[wav,sr] = audioread('speech.wav');
wav = resample(wav, P.audio_sr, sr);
wav = wav(1:P.audio_sr);

%% Cochleogram

addpath('Sound_Texture_Synthesis_Toolbox');

% cochleogram
[coch, P, R] = wav2coch_without_filts(wav, P);

% plot cochleogram
figure;
plot_cochleogram(coch, P.f, P.t);

%% Single subband

clc;
% close all;
addpath('2DFT');

% parameters
spec_mod_rate = 2; % cycles per octave
temp_mod_rate = 8; % Hertz
spec_mod_lowpass = false;
temp_mod_lowpass = false;
complex_filters = false;
causal = true;
fourier_domain = false;
separable = false;
spec_BW = 1;
temp_BW = 1;
spec_wavelet = 'morlet';

% result of convolution with the impulse response
coch_subband = coch2filtcoch(coch, spec_mod_rate, temp_mod_rate, P, ...
    spec_mod_lowpass, temp_mod_lowpass, complex_filters, separable, causal, ...
    fourier_domain, spec_BW, temp_BW, spec_wavelet);

% plot subband
figure;
plot_cochleogram(real(coch_subband), P.f, P.t);

% impulse response
Hts_FT = filt_spectemp_mod(...
    spec_mod_rate, temp_mod_rate, ...
    size(coch,2), size(coch,1), P, spec_mod_lowpass, ...
    temp_mod_lowpass, 0, 0, complex_filters, separable, causal, ...
    spec_BW, temp_BW, spec_wavelet);

% impulse response in signal domain
Hts_signal = ifft2(Hts_FT)';
Hts_signal = flipud(circshift(Hts_signal, [ceil(size(coch,2)/2-1), 0]));

% plot impulse response
figure;
imagesc((real(Hts_signal)), ...
    max(abs(real(Hts_signal(:))))*[-1 1]);

%% Multiple subbands

% include negative and positive rates
% P.temp_mod_rates = unique([-P.temp_mod_rates, P.temp_mod_rates]);

% compute subbands
% time x frequency x spectral modulation x temporal modulation x orientation
filtcoch_allsubbands = coch2filtcoch_allsubbands(coch, P, true);


% time/sound x feature
% F
% [U, S, V] = svd(F, 'econ');
% U*S*V' = F
% U*S = F * V
% U*S is the feature that you want
% eigenvalues = diag(S).^2
% eig_cum = cumsum(eigenvalues) / sum(eigenvalues)

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

