function plot_2DFT_orig_and_synth(...
    coch_orig, coch_synth, P, output_directory, fname_without_extension)

% add directory with useful 2D FT scripts
if ~exist('fft_freqs_from_siglen.m', 'file')
    directory_containing_this_file = fileparts(which(mfilename));
    addpath(genpath([directory_containing_this_file '/2DFT']));
end

FT_coch_orig = fft2(pad_coch_freq(coch_orig,P));
FT_coch_synth = fft2(pad_coch_freq(coch_synth,P));

% dimensions of X
[T,S] = size(FT_coch_orig);

% shift FT such that DC is in the middle and the nyquist is last
FT_coch_orig_shifted = fftshift_nyqlast(FT_coch_orig);
FT_coch_synth_shifted = fftshift_nyqlast(FT_coch_synth);

FT_orig_mag = 20*log10(abs(FT_coch_orig_shifted))';
FT_synth_mag = 20*log10(abs(FT_coch_synth_shifted))';
FT_orig_phase = unwrap(angle(FT_coch_orig_shifted)');
FT_synth_phase = unwrap(angle(FT_coch_synth_shifted)');
max_amp = max(max(FT_orig_mag(:)), max(FT_synth_mag(:)));
max_phase = max(max(FT_orig_phase(:)), max(FT_synth_phase(:)));
min_phase = min(min(FT_orig_phase(:)), min(FT_synth_phase(:)));

% shifted frequency vectors for the columns and rows
f_T = fftshift_nyqlast(fft_freqs_from_siglen(T, P.env_sr));
f_S = fftshift_nyqlast(fft_freqs_from_siglen(S, 1/P.logf_spacing));

% plot magnitude spectrum
figure;
set(gcf, 'Position', [300 300 800 600]);

% magnitude
subplot(2,2,1);
imagesc(f_T, f_S, FT_orig_mag);
title('Orig Magnitude');
caxis(max_amp + [-100 0]);

% magnitude
subplot(2,2,2);
imagesc(f_T, f_S, FT_synth_mag);
title('Synth Magnitude');
caxis(max_amp + [-100 0]);

% phase
subplot(2,2,3);
imagesc(f_T, f_S, FT_orig_phase);
title('Orig Phase');
caxis([min_phase, max_phase]);

% phase
subplot(2,2,4);
imagesc(f_T, f_S, FT_synth_phase);
title('Synth Phase');
caxis([min_phase, max_phase]);

% format axes
for i = 1:4
    subplot(2,2,i);
    set(gca,'YDir','normal');
    xlabel('Temp Mod Rate (Hz)'); ylabel('Spec Mod Scale (cyc/oct)');
    colorbar;
    x_ticks_to_plot = flip(T:-round(T/5):1);
    y_ticks_to_plot = flip(S:-round(S/5):1);
    set(gca, 'XTick', f_T(x_ticks_to_plot), 'YTick', f_S(y_ticks_to_plot));
    colormap('gray')
end

fig_fname = [output_directory '/' fname_without_extension '_2DFT'];
set(gcf, 'PaperSize', [10 8]);
set(gcf, 'PaperPosition', [0.25 0.25 9.5 7.5]);
print([fig_fname '.pdf'],'-dpdf');
print([fig_fname '.png'],'-dpng', '-r100');
