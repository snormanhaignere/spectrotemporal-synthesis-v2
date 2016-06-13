function plot_2DFT(FTX, sr)

% Plots the 2D fourier transform of a signal. 
% 
% -- Input --
% 
% FTX: Complex-valued 2D fourier transform of a signal (e.g. fft2(randn(10)))
% 
% sr (optional): [sr_col, sr_row] vector specifying the sampling rate of the column and
% rows, default = [1 1];
% 
% -- Example --
% 
% plot_2DFT(fft2(randn(6,6)))

if nargin < 2
    sr = [1,1];
end

% dimensions of X
[M,N] = size(FTX);

% shift FT such that DC is in the middle and the nyquist is last
FTX_shifted = fftshift_nyqlast(FTX);

% shifted frequency vectors for the columns and rows
f_M = fftshift_nyqlast(fft_freqs_from_siglen(M, sr(1)));
f_N = fftshift_nyqlast(fft_freqs_from_siglen(N, sr(2)));

% plot magnitude spectrum
figure;
set(gcf, 'Position', [300 300 800 300]);

% magnitude
subplot(1,2,1);
imagesc(f_N, f_M, 20*log10(abs(FTX_shifted)));
title('Magnitude');

% phase
subplot(1,2,2);
imagesc(f_N, f_M, angle(FTX_shifted));
title('Phase');

% format axes
for i = 1:2
    subplot(1,2,i);
    set(gca,'YDir','normal');
    xlabel('Freq'); ylabel('Freq');
    colorbar;
    x_ticks_to_plot = flip(N:-round(N/5):1);
    y_ticks_to_plot = flip(M:-round(M/5):1);
    set(gca, 'XTick', f_N(x_ticks_to_plot), 'YTick', f_M(y_ticks_to_plot));
    colormap('gray')
end