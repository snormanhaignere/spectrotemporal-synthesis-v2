function plot_filtcoch_subbands(coch_orig, coch_synth, ti_center, ...
    P, output_directory, fname_without_extension)

% 2DFT of original and synthetic cochleogram
FT_padded_coch_orig = fft2(pad_coch_freq(coch_orig, P));
FT_padded_coch_synth = fft2(pad_coch_freq(coch_synth, P));
[T,F] = size(FT_padded_coch_orig);

temp_mods_to_plot = select_mod_rates(P.temp_mod_rates);
spec_mods_to_plot = select_mod_rates(P.spec_mod_rates);
freqs_to_plot = [200 1000 6000];

% audio frequency, temporal modulation rate, and spectral scale
example_filt_params(1,:) = [freqs_to_plot(1), temp_mods_to_plot(1), spec_mods_to_plot(1)];
example_filt_params(2,:) = [freqs_to_plot(1), temp_mods_to_plot(2), spec_mods_to_plot(2)];
example_filt_params(3,:) = [freqs_to_plot(1), temp_mods_to_plot(3), spec_mods_to_plot(3)];
example_filt_params(4,:) = [freqs_to_plot(2), temp_mods_to_plot(1), spec_mods_to_plot(1)];
example_filt_params(5,:) = [freqs_to_plot(2), temp_mods_to_plot(2), spec_mods_to_plot(2)];
example_filt_params(6,:) = [freqs_to_plot(2), temp_mods_to_plot(3), spec_mods_to_plot(3)];
example_filt_params(7,:) = [freqs_to_plot(3), temp_mods_to_plot(1), spec_mods_to_plot(1)];
example_filt_params(8,:) = [freqs_to_plot(3), temp_mods_to_plot(2), spec_mods_to_plot(2)];
example_filt_params(9,:) = [freqs_to_plot(3), temp_mods_to_plot(2), spec_mods_to_plot(3)];

figure;
set(gcf, 'Position', [200 200 800 800]);
for i = 1:size(example_filt_params,1)
    
    filt_params = example_filt_params(i,:);
    
    % transfer function of spectrotemporal filter
    Hts = filt_spectemp_mod(...
        filt_params(3), filt_params(2), F, T, P);
    
    % apply transfer function
    filtcoch_padded_orig = real(ifft2(FT_padded_coch_orig .* Hts));
    filtcoch_orig = remove_freq_pad(filtcoch_padded_orig, P);
    filtcoch_padded_synth = real(ifft2(FT_padded_coch_synth .* Hts));
    filtcoch_synth = remove_freq_pad(filtcoch_padded_synth, P);
    
    % filter response matrix
    [~,freq_index] = nearest_elem(P.f, filt_params(1));
    filt_values = [filtcoch_orig(:,freq_index), filtcoch_synth(:,freq_index)];
    
    % measure the histogram
    % conver to probability mass function
    n_bins = length(ti_center)/10;
    [counts,bins] = hist(filt_values, n_bins);
    pmf = counts/length(ti_center);
    
    % bring up the appropriate subplot
    subplot(3, 3, i);
    
    % plot the histogram
    plot(bins, pmf, 'LineWidth', 2);
    xlabel('Subband value'); ylabel('Probability');
    set(gca, 'FontSize', 10);
    xlim([min(bins), max(bins)]);
    
    % legend
    legend('Orig', 'Synth');
    
    % title of subplot
    title(sprintf(...
        '%.0f Hz, %.3f Hz, %.3f cyc/oct', ...
        P.f(freq_index), filt_params(2), filt_params(3) ));
    
end

clear filt_params;


fig_fname = [output_directory '/' fname_without_extension ...
    '_spectemp_hists'];
set(gcf, 'PaperSize', [10 10]);
set(gcf, 'PaperPosition', [0.25 0.25 9.5 9.5]);
print([fig_fname '.pdf'],'-dpdf');
print([fig_fname '.png'],'-dpng', '-r100');

function mod_rates = select_mod_rates(mod_rates)

mod_rates = mod_rates(mod_rates>0);
n_mod_rates = length(mod_rates);
mod_rates = mod_rates([1 round(n_mod_rates/2), end]);