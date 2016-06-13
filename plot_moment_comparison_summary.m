function plot_moment_comparison_summary(C, output_directory, fname_without_extension)

% Plots stats comparing moments of the filter outputs for original and synthetic
% stimuli. These figures plot the median SNR or the correlation across all filters
% for different stimulus moments. Values are plotted as a function of the
% synthesis iteration, to assess the improvement in matching.

% moments
moment_names = {'Mean', 'Var', 'Skew', 'Kurt'};

% number of iterations run by the algorithm
n_iterations = size(C.coch_env_SNR_median,1) - 1;

% type of statistic and corresponding units
comparison_statistic = {'SNR_median', 'corr'};
stat_units = {'dB', 'signed r^2'};

% remove spaces and capitalize
comparison_statistic_formatted = upper(strrep(comparison_statistic, '_', ' '));

% where to place the legend on the figure
legend_location = 'EastOutside';

figure;
set(gcf, 'Position', [0 0 1000 800]);
n_comparison_statistics = length(comparison_statistic);
n_subplots_per_comparison_stat = 4;
for i = 1:n_comparison_statistics
    
    %% -- cochlear stats --
    
    subplot_index = i;
    subplot(n_subplots_per_comparison_stat, n_comparison_statistics, subplot_index);
    plot(0:n_iterations, C.(['coch_env_' comparison_statistic{i}]), ...
        'LineWidth', 2 );
    legend(moment_names, 'Location', legend_location);
    title(['Coch Env: ' comparison_statistic_formatted{i}]);
    
    %% -- temporal modulation stats --
    
    subplot_index = i + n_comparison_statistics;
    subplot(n_subplots_per_comparison_stat, n_comparison_statistics, subplot_index);
    hold on;
    
    % envelopes stats
    plot(0:n_iterations, C.(['temp_mod_' comparison_statistic{i}]), ...
        'LineWidth', 2 );
    
    % legend and tile
    legend(moment_names, 'Location', legend_location);
    title(['Temp Mod: ' comparison_statistic_formatted{i}]);
    
    %% -- spectral modulation stats --
    
    subplot_index = i + 2*n_comparison_statistics;
    subplot(n_subplots_per_comparison_stat, n_comparison_statistics, subplot_index);
    hold on;
    
    % envelopes stats
    plot(0:n_iterations, C.(['spec_mod_' comparison_statistic{i}]), ...
        'LineWidth', 2 );
    
    % legend and tile
    legend(moment_names, 'Location', legend_location);
    title(['Spec Mod: ' comparison_statistic_formatted{i}]);
    
    %% -- spectrotemporal modulation stats --
    
    subplot_index = i + 3*n_comparison_statistics;
    subplot(n_subplots_per_comparison_stat, n_comparison_statistics, subplot_index);
    hold on;
    
    % envlope stats
    plot(0:n_iterations, C.(['spectemp_mod_' comparison_statistic{i}]), ...
        'LineWidth', 2);
    
    
    % legend and title
    legend(moment_names, 'Location', legend_location);
    title(['SpecTemp Mod: ' comparison_statistic_formatted{i}]);
    
    %% generic formatting for all subplots
    
    for j = 1:n_subplots_per_comparison_stat
        subplot_index = i + (j-1)*2;
        subplot(n_subplots_per_comparison_stat, n_comparison_statistics, subplot_index);
        xlabel('Iterations');
        ylabel(stat_units{i});
        if n_iterations == 0
            xlim([-0.5 0.5]);
        else
            xlim([0, n_iterations]);
        end
        box off;
    end
    

end

% save figure
fig_fname = [output_directory '/' fname_without_extension '_summary_moment_comparisons'];
set(gcf, 'PaperSize', [11 8.5]);
set(gcf, 'PaperPosition', [0.25 0.25 10.5 8]);
print([fig_fname '.pdf'],'-dpdf');
print([fig_fname '.png'],'-dpng', '-r100');



