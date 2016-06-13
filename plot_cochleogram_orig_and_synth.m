function plot_cochleogram_orig_and_synth(coch_orig, coch_synth, ...
    P, output_directory, fname_without_extension)

% cochleogram color bounds
X = [coch_synth(:); coch_orig(:)];
bounds = [min(X(:)), max(X(:))];
clear X;

% plot
figure;
set(gcf, 'Position', [0 0 1200 300]);
subplot(1,2,1);
plot_cochleogram(coch_orig,P.f,P.t);
caxis(bounds);
title('Original');
subplot(1,2,2);
plot_cochleogram(coch_synth,P.f,P.t);
caxis(bounds);
title('Synthetic');
drawnow;

fname_full_path = [output_directory '/' fname_without_extension '_cochleograms'];
set(gcf, 'PaperSize', [10 3]);
set(gcf, 'PaperPosition', [0.25 0.25 9.5 2.5]);
print([fname_full_path '.pdf'], '-dpdf');
print([fname_full_path '.png'], '-dpng', '-r100');

