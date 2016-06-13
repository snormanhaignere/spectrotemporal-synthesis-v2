function plot_moment_scatter(...
    M_orig, M_synth, P, output_directory, fname_without_extension)

n_moments = size(M_orig.coch_env,2);
moment_names = {'Mean', 'Var', 'Skew', 'Kurt'};

% cochlear marginals
figure;
set(gcf, 'Position', [0 0 800 800]);
for i = 1:n_moments
    subplot(4,4,i);
    plot_synth_vs_orig_scatter(...
        M_orig.coch_env(:,i), M_synth.coch_env(:,i), ...
        ['Coch ' moment_names{i}]);
end

% fig_fname = [output_directory '/' fname_without_extension '_moments_coch_scatter'];
% save_moment_figure(gcf, fig_fname);

% temporal modulation
for i = 1:n_moments
    subplot(4,4,i + 4);
    plot_synth_vs_orig_scatter(...
        M_orig.temp_mod(:,:,i), M_synth.temp_mod(:,:,i), ...
        ['Tempmod ' moment_names{i}]);
end
% fig_fname = [output_directory '/' fname_without_extension '_moments_tempmod_scatter'];
% save_moment_figure(gcf, fig_fname);

% spectral modulation
for i = 1:n_moments
    subplot(4,4,i + 8);
    plot_synth_vs_orig_scatter(...
        M_orig.spec_mod(:,:,i), M_synth.spec_mod(:,:,i), ...
        ['Specmod ' moment_names{i}]);
end
% fig_fname = [output_directory '/' fname_without_extension '_moments_specmod_scatter'];
% save_moment_figure(gcf, fig_fname);

% spectrotemporal modulation
for i = 1:n_moments
    subplot(4,4,i + 12);
    plot_synth_vs_orig_scatter(...
        M_orig.spectemp_mod(:,:,:,i), M_synth.spectemp_mod(:,:,:,i), ...
        ['Spectempmod ' moment_names{i}]);
end
fig_fname = [output_directory '/' fname_without_extension '_moments_scatter'];
save_moment_figure(gcf, fig_fname);



function save_moment_figure(figh, fig_fname)

set(figh, 'PaperSize', [11 11]);
set(figh, 'PaperPosition', [0.25 0.25 10.5 10.5]);
print([fig_fname '.pdf'],'-dpdf');
print([fig_fname '.png'],'-dpng', '-r100');

function plot_synth_vs_orig_scatter(orig, synth, titlestring)

orig = orig(:);
synth = synth(:);

% plot
plot(orig, synth, 'o');

% bounds
X = [orig, synth];
bounds = [min(X(:)), max(X(:))];
clear X;
xlim(bounds); ylim(bounds);

% x = y line
hold on;
plot(bounds, bounds, 'r--');

% axis labels
xlabel('Orig Mag');
ylabel('Synth Mag');

% title with correlation
r = corr(orig, synth);
signed_r2 = sign(r) .* r.^2;
title(sprintf('%s\nr^2 = %.2f', titlestring, signed_r2));
clear bounds r signed_r2;