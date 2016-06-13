function plot_filtered_cochleograms(coch_orig, coch_synth, ...
    P, output_directory, fname_without_extension)

% 2DFT of original and synthetic cochleogram
FT_padded_coch_orig = fft2(pad_coch_freq(coch_orig, P));
FT_padded_coch_synth = fft2(pad_coch_freq(coch_synth, P));
[T,F] = size(FT_padded_coch_orig);

figure;
set(gcf, 'Position', [0 0 1200 800]);

spec_mod_to_plot = select_mod_rates(P.spec_mod_rates);
temp_mod_to_plot = select_mod_rates(P.temp_mod_rates);

for i = 1:length(spec_mod_to_plot)
    for j = 1:length(temp_mod_to_plot)
        
        % transfer function of spectrotemporal filter
        Hts = filt_spectemp_mod(...
            spec_mod_to_plot(i), temp_mod_to_plot(j), F, T, P);
                                
        % apply transfer function
        filtcoch_padded_orig = real(ifft2(FT_padded_coch_orig .* Hts));
        filtcoch_orig = remove_freq_pad(filtcoch_padded_orig, P);
        filtcoch_padded_synth = real(ifft2(FT_padded_coch_synth .* Hts));
        filtcoch_synth = remove_freq_pad(filtcoch_padded_synth, P);
        filtcochs = cat(3, filtcoch_orig, filtcoch_synth);
        
        orig_or_synth = {'orig', 'synth'};
        for k = 1:2
                        
            % select subplot
            n_row = length(spec_mod_to_plot);
            n_col = length(temp_mod_to_plot) * 2;
            subplot_index = subplot_sub2ind(n_row, n_col, i, (j-1) * 2 + k);
            subplot(n_row, n_col, subplot_index);
            
            % plot
            plot_cochleogram(filtcochs(:,:,k), P.f, P.t);
            set(gca, 'FontSize', 6);
            clear X;
            
            % set color bound
            caxis([min(filtcochs(:)), max(filtcochs(:))]);
            clear X;
            
            % title
            % if there is no spectral modulation just plot temporal rate
            if isnan(spec_mod_to_plot(i))
                title(sprintf([...
                    orig_or_synth{k} ', \n' ...
                    num2str(temp_mod_to_plot(j)) ' Hz']));
            else
                title(sprintf([...
                    orig_or_synth{k} ', \n' ...
                    num2str(spec_mod_to_plot(i)) ' cyc/oct, '...
                    num2str(temp_mod_to_plot(j)) ' Hz']));
            end
        end
    end
end

% save
fig_fname = [output_directory '/' fname_without_extension '_filt_coch'];
set(gcf, 'PaperSize', [11 8]);
set(gcf, 'PaperPosition', [0.25 0.25 10.5 7.5]);
print([fig_fname '.pdf'],'-dpdf');
print([fig_fname '.png'],'-dpng','-r200');

function mod_rates = select_mod_rates(mod_rates)

mod_rates = mod_rates(mod_rates>0);
n_mod_rates = length(mod_rates);
mod_rates = mod_rates([1 round(n_mod_rates/2), end]);

