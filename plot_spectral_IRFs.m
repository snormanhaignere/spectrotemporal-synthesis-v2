function plot_spectral_IRFs

% Plots spectral impulse responses for the filters from the model

clc;
close all;

% find the directory containing this file
name_of_this_file = 'plot_spectral_IRFs';%mfilename; % 'plot_temporal_irfs';%mfilename;
directory_containing_this_file = fileparts(which(name_of_this_file));

% Modified Shamma toolbox
addpath(genpath([directory_containing_this_file ...
    '/nsltools_texture_synthesis']));

% McDermott Texture toolbox
addpath(genpath([directory_containing_this_file ...
    '/Sound_Texture_Synthesis_Toolbox']));

% color maps 
addpath(genpath([directory_containing_this_file ...
    '/cbrewer']));

% directory to save results
figure_directory = [directory_containing_this_file '/IRFs'];
if ~exist(figure_directory,'dir')
    mkdir(figure_directory);
end

P = synthesis_parameters_default;
P.spec_mod_rates = P.spec_mod_rates(P.spec_mod_rates>0);

%% Temporal impulse responses

for i = 1:length(P.spec_mod_rates)
    
    figure;
    
    % period of the spectral modulation rate in octaves
    period_oct = 1/P.spec_mod_rates(i);
    
    % number of samples
    % temporal sampling rate in Hz
    specral_sr_cyc_per_oct = 100;
    n_spectral_smps = ceil(specral_sr_cyc_per_oct*period_oct*4); % number of samples
    
    % transfer function of the filter
    filt_tf = filt_spec_mod(...
        P.spec_mod_rates(i), n_spectral_smps, specral_sr_cyc_per_oct, 0, 0);
    
    % impulse response
    filt_irf = circshift(real(ifft(filt_tf)), n_spectral_smps/2-1);
    
    % corresponding frequencies for the spectral filter
    freqs = ...
        [(0:1:n_spectral_smps/2), (-(n_spectral_smps/2-1):1:-1)]' ...
        /specral_sr_cyc_per_oct;
    freqs = circshift(freqs, n_spectral_smps/2-1);
    
    % plot
    plot(freqs, filt_irf);
    
    % y-limits
    yL = ylim;
    yL = max(abs(yL))*[-1 1];
    ylim(yL);
    
    % plot period lines
    hold on;
    for j = -2:1:2
        plot(j*period_oct * [1 1], yL, 'k--')
    end
    
    % x-limits
    xlim([freqs(1),freqs(end)]);

    % font size
    set(gca, 'FontName', 'Helvetica','FontSize',20);
    
    % axis labels
    xlabel('Frequency (octaves)'); ylabel('Filter Response');
    
    % title string and figure string
    title_string = 'Period = %.3f oct';
    figure_string = ['spec-filt-' num2str(P.spec_mod_rates(i))  'cyc_per_oct'];
    if lowpass
        title_string = [title_string, '\nLowpass']; %#ok<*AGROW>
        figure_string = [figure_string '_lowpass'];
    end
    if highpass
        title_string = [title_string, '\nLowpass'];
        figure_string = [figure_string '_highpass'];
    end
    
    % title
    title(sprintf(title_string, period_oct));
    
    % save figure;    
    figure_size = [10 8];
    set(gcf, 'PaperSize', figure_size);
    set(gcf, 'PaperPosition', [0.25, 0.25, figure_size-0.5]);
    print([figure_directory '/' figure_string '.pdf'],'-dpdf');
    print([figure_directory '/' figure_string '.png'],'-dpng', '-r100');   
    
end


