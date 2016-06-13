function plot_freq_by_tempmod_image(...
    freq_by_tempmod_matrix, audio_freqs, temp_mod_rates,  varargin)

% plots an image of a statistic parameterized by audio frequency and temporal
% modulation
% 
% audio_freqs = logspace(log10(20),log10(10000),100);
% temp_mod_rates = logspace(log10(0.5),log10(128),9);
% im = randn( length(audio_freqs), length(temp_mod_rates) );
% plot_freq_by_tempmod_image(im, audio_freqs, temp_mod_rates)

% bring up a figure if it doesn't exist
gcf;

% plot
imagesc(flipud(freq_by_tempmod_matrix));

% xtick values
if optInputs(varargin, 'XTick')
    xtick = varargin{optInputs(varargin, 'XTick')+1};
    set(gca, 'XTick', xtick);
end

% format x axis
xtick = round(get(gca, 'XTick'));
xtick = intersect(xtick, 1:length(temp_mod_rates));
set(gca, 'XTick', xtick, 'XTickLabel', temp_mod_rates(xtick))
xlabel('Rate (Hz)');

% format y axis
freqs_flip = flip(audio_freqs);
ytick = round(get(gca, 'YTick'));
ytick = intersect(ytick, 1:length(freqs_flip));
set(gca, 'YTick', ytick, 'YTickLabel', num2cellstr( freqs_flip(ytick), '%.0f' ));
ylabel('Frequency (Hz)');

set(gca, 'FontSize', 8);



