function plot_spectral_TFs

% Plots temporal impulse responses for the filters from the model

clc;
% close all;

% find the directory containing this file
name_of_this_file = 'plot_spectral_TFs'; %mfilename; % 'plot_temporal_irfs';%mfilename;
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
figure_directory = [directory_containing_this_file '/TFs'];
if ~exist(figure_directory,'dir')
    mkdir(figure_directory);
end

P = measurement_parameters_default;
wavelet = 'mexicanhat';
BW = 1;

%% Temporal impulse responses

figure;
for i = 1:length(P.spec_mod_rates)
            
    % number of samples
    % temporal sampling rate in Hz
    specral_sr_cyc_per_oct = 24;
    n_oct = log2(10e3/20)*2;
    n_spectral_smps = ceil(specral_sr_cyc_per_oct*n_oct); % number of samples
    
    % transfer function of the filter
    filt_tf = filt_spec_mod(...
        P.spec_mod_rates(i), n_spectral_smps, specral_sr_cyc_per_oct, ...
        P.spec_mod_lowpass(i), 0, BW, wavelet);
    
    % nyquist
    max_pos_freq = ceil((n_spectral_smps+1)/2);
    
    % modulation frequency vector in Hz
    mod_freqs_cyc_per_oct = specral_sr_cyc_per_oct*(0:max_pos_freq-1)/n_spectral_smps;
    
    % power of the filter
    filt_pow_dB = 20*log10(abs(filt_tf(1:max_pos_freq)));
    filt_phase = unwrap(angle(filt_tf(1:max_pos_freq)));
    
    % calculate 3dB down bandwidth
    if true%~P.spec_mod_lowpass(i) %~lowpass && ~highpass
        mod_freqs_above_cf = mod_freqs_cyc_per_oct(mod_freqs_cyc_per_oct > P.spec_mod_rates(i));
        filt_pow_above_cf = filt_pow_dB(mod_freqs_cyc_per_oct > P.spec_mod_rates(i));
        mod_freqs_below_cf = mod_freqs_cyc_per_oct(mod_freqs_cyc_per_oct < P.spec_mod_rates(i));
        filt_pow_below_cf = filt_pow_dB(mod_freqs_cyc_per_oct < P.spec_mod_rates(i));
        
        % frequencies 3dB above and below the cf
        [~,xi] = min(abs(filt_pow_above_cf - (-3)));
        down_3dB_above = mod_freqs_above_cf(xi);
        
        [~,xi] = min(abs(filt_pow_below_cf - (-3)));
        down_3dB_below = mod_freqs_below_cf(xi);
        
        % Q factor
        Q = P.spec_mod_rates(i) / (down_3dB_above - down_3dB_below);
        fprintf('%.3f, Q = %.3f\n', P.spec_mod_rates(i), Q);
    end
    
    % plot
    for k = 1:2
        
        subplot(2,1,k);
        if k == 1
            h = semilogx(mod_freqs_cyc_per_oct, filt_pow_dB);
            hold on;
        elseif k == 2;
            h = semilogx(mod_freqs_cyc_per_oct, filt_phase);
            hold on;
        end
        
        % y-limits
        if k == 1
            yL = [-20 0];
            ylim(yL);
        elseif k == 2
            yL = [-6 6];
            ylim(yL);
        end
               
        % plot octave lines
        plot( P.spec_mod_rates(i) * [1 1], yL, '--', 'Color', get(h, 'Color'));
        
        % x-limits
        xlim([mod_freqs_cyc_per_oct(2),mod_freqs_cyc_per_oct(end)]);
        
        
        if k == 1
            xL = xlim;
            plot(xL, [-3 -3], 'k--');
        end
        
        
        set(gca, 'XTick', [0.1 1 10]);
        
        % font size
        set(gca, 'FontName', 'Helvetica','FontSize',20);
        
        % axis labels
        xlabel('Modulation Frequency (Hz)');
        
        if k == 2
            ylabel('Phase (Rad)');
        else
            ylabel('Magnitude (dB)');
        end
        
    end
    
end

figure_size = [10 8];
fname_full_path = [figure_directory '/spec-filts'];
set(gcf, 'PaperSize', figure_size);
set(gcf, 'PaperPosition', [0.25, 0.25, figure_size-0.5]);
print([fname_full_path, '.pdf'],'-dpdf');
print([fname_full_path, '.png'],'-dpng', '-r100');

% 
% %%
%     % title string and figure string
%     title_string = 'Period = %.3f sec';
%     figure_string = ['temp-filt-' num2str(P.spec_mod_rates(i))  'Hz'];
%     if lowpass
%         title_string = [title_string, '\nLowpass']; %#ok<*AGROW>
%         figure_string = [figure_string '_lowpass'];
%     end
%     if highpass
%         title_string = [title_string, '\nLowpass'];
%         figure_string = [figure_string '_highpass'];
%     end
%     
%     % title
%     title(sprintf(title_string, period_sec));    
%     
%     % save figure;
%     box off;
%     export_fig(...
%         [figure_directory '/' figure_string '.pdf'],...
%         '-nocrop','-transparent', '-pdf');
%     
% 
