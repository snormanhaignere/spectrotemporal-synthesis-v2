function plot_moments(M_orig, M_synth, P, output_directory, fname_without_extension)

% Plots moments of cochlear and modulation filters for both the original and
% synthetic sounds. Moments are assumed to be computed by the function
% filter_moments_from_coch.

n_moments = size(M_orig.coch_env,2);
moment_names = {'Mean', 'Var', 'Skew', 'Kurt'};

% combine fields for original and synthetic
allfields = fieldnames(M_orig);
for i = 1:length(allfields)
    dims = size(M_orig.(allfields{i}));
    if dims(end) == 1
        ndims = length(dims)-1;
    else
        ndims = length(dims);
    end
    
    M.(allfields{i}) = cat(ndims+1, M_orig.(allfields{i}), M_synth.(allfields{i}));
end
clear dims ndims;

%% Cochlear envelopes

figure;
set(gcf, 'Position', [0 0 1200 250]);
for i = 1:n_moments
    subplot(1,4,i);
    semilogx(P.f, [M_orig.coch_env(:,i), M_synth.coch_env(:,i)], ...
        'LineWidth',2);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    xlim([min(P.f), max(P.f)]);
    title(moment_names{i});
    legend({'Orig', 'Synth'});
end

drawnow;
fig_fname = [output_directory '/' fname_without_extension '_moments_coch'];
set(gcf, 'PaperSize', [11 3]);
set(gcf, 'PaperPosition', [0.25 0.25 10.5 2.5]);
print([fig_fname '.pdf'],'-dpdf');
print([fig_fname '.png'],'-dpng', '-r100');

%% Temporal modulation

% combine all of the temporal modulation stats
% useful to plot them all in one loop
tempmod_stats = {M.temp_mod};
tempmod_stat_names = {'Temp Mod'};

pos_temp_mod_rates = P.temp_mod_rates(P.temp_mod_rates>0);

figure;
set(gcf, 'Position', [0 0 600 800]);
for i = 1:length(tempmod_stats)
    for j = 1:n_moments
        orig_or_synth_flag = {'Orig', 'Synth'};
        for k = 1:length(orig_or_synth_flag)
            
            % subplot index
            n_rows = n_moments;
            n_cols = length(orig_or_synth_flag) * length(tempmod_stats);
            subplot_index = subplot_sub2ind(...
                n_rows, n_cols, j, k + (i-1) * length(orig_or_synth_flag));
                        
            % plot
            subplot(n_rows, n_cols, subplot_index);
            plot_freq_by_tempmod_image(...
                tempmod_stats{i}(:,:,j,k)', P.f, pos_temp_mod_rates);
            
            % font size
            set(gca, 'FontSize', 7);
            
            % title
            title([tempmod_stat_names{i} ' ' moment_names{j} ...
                ' ' orig_or_synth_flag{k}]);
            
            % fix color range for original and synthetic
            X = tempmod_stats{i}(:,:,j,:); % all values for this moment
            caxis( [min(X(:)), max(X(:))] );
            clear X;
            
        end
    end
end
clear tempmod_stats tempmod_stats;
drawnow;

fig_fname = [output_directory '/' fname_without_extension '_moments_tempmod'];
set(gcf, 'PaperSize', [5 8]);
set(gcf, 'PaperPosition', [0.25 0.25 4.5 7.5]);
print([fig_fname '.pdf'],'-dpdf');
print([fig_fname '.png'],'-dpng', '-r100');

%% Spectrotemporal Modulation, 1 figure per spectral scale

% combine all of the temporal modulation stats
% useful to plot them all in one loop
spectemp_stats = {M.spectemp_mod};
spectemp_stat_names = {'Spectemp Mod'};

for q = 1:length(P.spec_mod_rates)
    figure;
    set(gcf, 'Position', [0 0 600 800]);
    
    for j = 1:n_moments
        for i = 1:length(spectemp_stats)
            
            orig_or_synth_flag = {'Orig', 'Synth'};
            for k = 1:length(orig_or_synth_flag)
                
                % subplot index
                n_rows = n_moments;
                n_cols = length(orig_or_synth_flag) * length(spectemp_stats);
                subplot_index = subplot_sub2ind(...
                    n_rows, n_cols, j, k + (i-1) * length(orig_or_synth_flag));
                
                
                % xticks to plot
                hop = round(length(pos_temp_mod_rates)/4);
                xticks = hop:hop:length(pos_temp_mod_rates)-1;
                xticks = [xticks, xticks + length(pos_temp_mod_rates)]; %#ok<AGROW>
                clear hop;
                
                % positive and negative rates
                temp_mod_rates_pos_and_neg = ...
                    [pos_temp_mod_rates, -pos_temp_mod_rates];
                
                % plot
                % spectemp_stats{i}:
                % scale x rate * 2 x freq x moments x orig/synth
                subplot(n_rows, n_cols, subplot_index);
                plot_freq_by_tempmod_image(...
                    squeeze(spectemp_stats{i}(q,:,:,j,k))', ...
                    P.f, temp_mod_rates_pos_and_neg, ...
                    'XTick', xticks);
                clear xticks
                
                % font size
                set(gca, 'FontSize', 7);
                
                % title
                title(sprintf([spectemp_stat_names{i} ' ' moment_names{j} ...
                    ' ' orig_or_synth_flag{k} ...
                    '\n' num2str(P.spec_mod_rates(q)) ' cyc/oct']));
                
                % fix color range for original and synthetic
                X = spectemp_stats{i}(q,:,:,j,:); % all values for this moment
                caxis( [min(X(:)), max(X(:))] );
                
            end
        end
    end
    
    drawnow;
    
    fig_fname = [output_directory '/' fname_without_extension ...
        '_moments_spectemp_' num2str(P.spec_mod_rates(q)) 'cyc'];
    set(gcf, 'PaperSize', [5 8]);
    set(gcf, 'PaperPosition', [0.25 0.25 4.5 7.5]);
    print([fig_fname '.pdf'],'-dpdf');
    print([fig_fname '.png'],'-dpng', '-r100');
    
end

clear spectemp_stats spectemp_stat_names;


%% Spectrotemporal Modulation, 1 figure per frequency

% combine all of the temporal modulation stats
% useful to plot them all in one loop
spectemp_stats = {M.spectemp_mod};
spectemp_stat_names = {'Spectemp Mod'};

freqs = [100 200 400 800 1600 3200 6400];

for q = 1:length(freqs)
    
    [~,freq_index] = min(abs(freqs(q) - P.f));
    
    figure;
    set(gcf, 'Position', [0 0 600 800]);
    
    for j = 1:n_moments
        
        for i = 1:length(spectemp_stats)
            
            orig_or_synth_flag = {'Orig', 'Synth'};
            for k = 1:length(orig_or_synth_flag)
                
                % subplot index
                n_rows = n_moments;
                n_cols = length(orig_or_synth_flag) * length(spectemp_stats);
                subplot_index = subplot_sub2ind(...
                    n_rows, n_cols, j, k + (i-1) * length(orig_or_synth_flag));
                
                % xticks to plot
                hop = round(length(pos_temp_mod_rates)/4);
                xticks = hop:hop:length(pos_temp_mod_rates)-1;
                xticks = [xticks, xticks + length(pos_temp_mod_rates)]; %#ok<AGROW>
                clear hop;
                
                % positive and negative rates
                temp_mod_rates_pos_and_neg = ...
                    [pos_temp_mod_rates, -pos_temp_mod_rates];
                
                % remove DC
                matrix_to_plot = spectemp_stats{i}(:,:,freq_index,j,k);
%                 matrix_to_plot(P.spec_mod_rates==0, temp_mod_rates_pos_and_neg==0) = NaN;
                
                % plot
                % spectemp_stats{i}:
                % scale x rate * 2 x freq x moments x orig/synth
                subplot(n_rows, n_cols, subplot_index);
                plot_spectemp_matrix(matrix_to_plot, P.spec_mod_rates,...
                    temp_mod_rates_pos_and_neg,'XTick', xticks);
                
                % font size
                set(gca, 'FontSize', 7);
                
                % title
                title(sprintf([spectemp_stat_names{i} ' ' moment_names{j} ...
                    ' ' orig_or_synth_flag{k} ...
                    '\n' num2str(freqs(q)) ' Hz']));
                
                % fix color range for original and synthetic
                % spectemp_stats{i}:
                % scale x rate * 2 x freq x moments x orig/synth
                X = spectemp_stats{i}(:,:,freq_index,j,:);
                caxis( [min(X(:)), max(X(:))] );
                
            end
        end
    end
    
    % save
    fig_fname = [output_directory '/' fname_without_extension ...
        '_moments_spectemp_' num2str(freqs(q)) 'Hz'];
    set(gcf, 'PaperSize', [5 8]);
    set(gcf, 'PaperPosition', [0.25 0.25 4.5 7.5]);
    print([fig_fname '.pdf'],'-dpdf');
    print([fig_fname '.png'],'-dpng', '-r100');
    
end

clear spectemp_stats spectemp_stat_names;




