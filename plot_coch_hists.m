function plot_coch_hists(coch_orig, coch_synth, ti_center, ...
    P, output_directory, fname_without_extension)

% plot_coch_hists(coch_orig, coch_synth, ...
%     P, output_directory, fname_without_extension)
% 
% Plots histograms for example cochlear filters

freqs_to_plot = [100, 200, 400, 800, 1600, 3200];

figure;
for i = 1:length(freqs_to_plot)
    
    f = freqs_to_plot(i);
    [~,f_index] = min(abs(P.f - f));
    
    n_bins = length(ti_center)/10;
    [counts,bins] = hist(...
        [coch_orig(ti_center,f_index), coch_synth(ti_center,f_index)], n_bins);
    pmf = counts/length(ti_center);
    
    subplot(2,3,i);
    plot(bins, pmf, 'LineWidth', 2);
    xlabel('Envelope Magnitude'); ylabel('Probability');
    set(gca, 'FontSize', 10);
    xlim([min(bins), max(bins)]);
    
    legend('Orig', 'Synth');
    
    title(sprintf('Coch Env Hist\n%.0f Hz', P.f(f_index)));
    
end

fig_fname = [output_directory '/' fname_without_extension ...
    '_coch_hists'];
set(gcf, 'PaperSize', [11 8]);
set(gcf, 'PaperPosition', [0.25 0.25 10.5 7.5]);
print([fig_fname '.pdf'],'-dpdf');
print([fig_fname '.png'],'-dpng', '-r100');

