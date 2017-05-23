% Simple script showing how to turn real filters into complex-valued filters,
% and how to use the complex-valued filters to compute envelopes
% 
% 2017-05-23

% parameters
P = synthesis_parameters_default;
T = 100;
F = 24*2;
spec_mod_rate = 2;
temp_mod_rate = 8;

% open a figure
figure;
set(gcf, 'Position', [200 200 300 800]);

% FT of real filter
FT_real_filter = filt_spectemp_mod(...
    spec_mod_rate, temp_mod_rate, F, T, P, 0, 0, 0, 0, false);

% FT of complex filter
FT_complex_filter = filt_spectemp_mod(...
    spec_mod_rate, temp_mod_rate, F, T, P, 0, 0, 0, 0, true);

% plot
subplot(4,1,1);
imagesc(real(ifft2(FT_real_filter)));
title('Real filter');
subplot(4,1,2);
imagesc(real(ifft2(FT_complex_filter)));
title('Real Part of Complex filter');
subplot(4,1,3);
imagesc(imag(ifft2(FT_complex_filter)));
title('Imaginary Part of Complex filter');
subplot(4,1,4);
imagesc(abs(ifft2(FT_complex_filter)));
title('Envelope');

% speech waveform
[wav, sr] = audioread('speech.wav');
wav = format_wav(wav, sr, P);

% cochleogram
[coch, P] = wav2coch_without_filts(wav, P);

% filtered cochleogram, subband and envelopes
spec_mod_rate = 2;
temp_mod_rate = 8;
filtcoch_complex = coch2filtcoch(coch, spec_mod_rate, temp_mod_rate, P, true);

% plot subband and envelopes
figure;
subplot(2,1,1);
plot_cochleogram(real(filtcoch_complex), P.f, P.t);
title('Subband');
subplot(2,1,2);
plot_cochleogram(abs(filtcoch_complex), P.f, P.t);
title('Envelope');