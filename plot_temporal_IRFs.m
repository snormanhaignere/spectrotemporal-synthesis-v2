function plot_temporal_IRFs

% Plots temporal impulse responses for the filters from the model

clc;
close all;

% find the directory containing this file
name_of_this_file = 'plot_temporal_IRFs';%mfilename;
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
P.temp_mod_rates = P.temp_mod_rates(P.temp_mod_rates>0);
lowpass = false;
highpass = false;
BW = 2;

%% Temporal impulse responses

for i = 1:length(P.temp_mod_rates)
    
    figure;
    
    % period of the modulation rate in seconds
    period_sec = 1/P.temp_mod_rates(i);
    
    % number of samples
    % temporal sampling rate in Hz
    temporal_sr_Hz = 1200;
    n_temporal_smps = round(period_sec * 4 * temporal_sr_Hz); % number of samples
   
    % transfer function of the filter
    filt_tf = filt_temp_mod(...
        P.temp_mod_rates(i), n_temporal_smps, temporal_sr_Hz, 0, 0, true, BW);
        
    % impulse response
    filt_irf = ifft(filt_tf);
    
    % time vector
    t = (1:n_temporal_smps)/temporal_sr_Hz;
    
    % plot
    plot(t, filt_irf);
    
    % y-limits
    yL = ylim;
    yL = max(abs(yL))*[-1 1];
    ylim(yL);
    
    % plot period lines
    hold on;
    for j = 1:5
        plot(j*period_sec * [1 1], yL, 'k--')
    end
    
    % x-limits
    xlim([t(1),t(end)]);

    % font size
    set(gca, 'FontName', 'Helvetica','FontSize',20);
    
    % axis labels
    xlabel('Time (sec)'); ylabel('Filter Response');
    
    % title string and figure string
    title_string = 'Period = %.3f sec';
    figure_string = ['temp-filt-' num2str(P.temp_mod_rates(i))  'Hz'];
    if lowpass
        title_string = [title_string, '\nLowpass']; %#ok<*AGROW>
        figure_string = [figure_string '_lowpass'];
    end
    if highpass
        title_string = [title_string, '\nLowpass'];
        figure_string = [figure_string '_highpass'];
    end
    
    % title
    title(sprintf(title_string, period_sec));    
    
    % save figure;    
    figure_size = [10 8];
    set(gcf, 'PaperSize', figure_size);
    set(gcf, 'PaperPosition', [0.25, 0.25, figure_size-0.5]);
    print([figure_directory '/' figure_string '.pdf'],'-dpdf');
    print([figure_directory '/' figure_string '.png'],'-dpng', '-r100');    
    
end


