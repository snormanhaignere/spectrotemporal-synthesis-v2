function plot_cochleogram(C, f, t)    

% plot_cochleogram(C, f, t) 
% 
% Plots matrix C as a cochleogram
% C is assumed to a time x frequency matrix
% t is a corresponding vector of temporal indices
% f is a corresponding vector of frequency values
% 
% % -- Example --
% 
% % texture toolbox
% addpath(genpath([pwd '/Sound_Texture_Synthesis_Toolbox']));
% 
% % read in waveform
% [wav,sr] = audioread([pwd '/speech.wav']);
% P.max_duration_sec = 1;
% wav = 0.1 * format_wav(wav, sr, P);
% 
% % filters
% P = default_synthesis_parameters;
% [audio_filts, audio_low_cutoff] = ...
%     make_erb_cos_filters(length(wav), P.audio_sr, ...
%     P.n_filts, P.lo_freq_hz, P.audio_sr/2);
% 
% % cochleogram
% [coch, P.f, P.t, R] = ...
%     wav2coch(wav, audio_filts, audio_low_cutoff, ...
%     P.audio_sr, P.env_sr, P.compression_factor, P.logf_spacing);
% 
% % plot the cochleogram
% plot_cochleogram(coch, P.f, P.t);

% dimensionality of the cochleogram
dims = size(C);
n_t = dims(1);
n_f = dims(2);

% plot
figure(gcf);
imagesc(flipud(C'));

% y-axis ticks
yi = round( linspace( 1, n_f, 5 ) );
f_flip = flip(f);
set(gca, 'YTick', yi, 'YTickLabel', num2cellstr( f_flip(yi), '%.0f' ) );
clear yi f_flip;

% x-axis ticks
xi = round( linspace(1, n_t, 5) );
set(gca, 'XTick', xi, 'XTickLabel',  num2cellstr( t(xi) , '%.1f') );
clear xi

% x-label and y-label
xlabel('Time (s)'); ylabel('Frequency (Hz)');
box off;