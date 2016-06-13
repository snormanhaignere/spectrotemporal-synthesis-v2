function plot_spectrotemporal_IRFs

clc;
close all;

% find the directory containing this file
name_of_this_file = 'plot_spectrotemporal_IRFs';%mfilename;
directory_containing_this_file = fileparts(which(name_of_this_file));

% McDermott Texture toolbox
addpath(genpath([directory_containing_this_file ...
    '/Sound_Texture_Synthesis_Toolbox']));

% color maps
addpath(genpath([directory_containing_this_file '/cbrewer']));

% directory to save results
figure_directory = [directory_containing_this_file '/IRFs'];
if ~exist(figure_directory,'dir')
    mkdir(figure_directory);
end

P = synthesis_parameters_default;

% positive and negative temporal modulation rates
temp_mod_rates_pos_and_neg = ...
    [NaN, P.temp_mod_rates, -P.temp_mod_rates(P.temp_mod_rates>0)];
spec_mod_rates = [NaN, P.spec_mod_rates];


for i = 1:length(spec_mod_rates)
    for j = 1:length(temp_mod_rates_pos_and_neg)
        
        
        % number of samples
        % temporal sampling rate in Hz
        n_sec = 1;
        temporal_sr_Hz = P.env_sr;
        n_temporal_smps = temporal_sr_Hz*n_sec; % number of samples
        
        % number of samples
        % spectral sampling rate in Hz
        n_octaves = 5;
        spectral_sr_cy_per_oct = 1/P.logf_spacing;
        n_spectral_smps = spectral_sr_cy_per_oct*n_octaves; % number of samples
        
        % spectral filter transfer function
        spectemp_filt_tf = filt_spectemp_mod(...
            spec_mod_rates(i), temp_mod_rates_pos_and_neg(j),...
            n_spectral_smps, n_temporal_smps, P, 0, 0, 0, 0);
        
        % impulse response
        spectemp_filt_irf = ifft2(spectemp_filt_tf)';
        spectemp_filt_irf = circshift(spectemp_filt_irf, [n_spectral_smps/2-1, 0]);
        
        % corresponding frequencies for the spectral filter
        freqs = ...
            [(0:1:n_spectral_smps/2), (-(n_spectral_smps/2-1):1:-1)]' ...
            /spectral_sr_cy_per_oct;
        freqs = circshift(freqs, n_spectral_smps/2-1);
        
        % plot image
        figure;
        imagesc(real(spectemp_filt_irf), ...
            max(abs(real(spectemp_filt_irf(:))))*[-1 1]);
        
        % modify color map
        colormap(cbrewer('div','RdBu',256));
        
        % time axis
        set(gca,'XTick',[0 n_temporal_smps/2 n_temporal_smps],...
            'XTickLabel',[0 n_temporal_smps/2 n_temporal_smps]/temporal_sr_Hz);
        xlabel('Time (seconds)');
        
        % frequency axis
        freqs_to_plot = -floor(n_octaves/2):1:floor(n_octaves/2);
        yticks = nan(size(freqs_to_plot));
        for k = 1:length(freqs_to_plot)
            [~,yticks(k)] = min(abs(freqs - freqs_to_plot(k)));
        end
        set(gca, 'YTick', yticks, ...
            'YTickLabel', freqs_to_plot);
        ylabel('Frequency (Octaves)');
        
        % set font size
        set(gca, 'FontName', 'Helvetica','FontSize',24)
        
        % sign string
        if sign(temp_mod_rates_pos_and_neg) < 0
            sign_string = 'neg';
        elseif sign(temp_mod_rates_pos_and_neg) > 0
            sign_string = 'pos';
        else
            sign_string = '';
        end
        
        figure_string = ['spectemp-filt-'...
            sign_string num2str(abs(temp_mod_rates_pos_and_neg(j))) ...
            'Hz-' num2str(spec_mod_rates(i)) 'cyc-per-oct'];
        
        % save figure;
        figure_size = [10 8];
        set(gcf, 'PaperSize', figure_size);
        set(gcf, 'PaperPosition', [0.25, 0.25, figure_size-0.5]);
        print([figure_directory '/' figure_string '.pdf'],'-dpdf');
        print([figure_directory '/' figure_string '.png'],'-dpng', '-r100');
        
        close all;
    end
end


