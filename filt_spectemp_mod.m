function Hts = filt_spectemp_mod(...
    spec_mod_rate, temp_mod_rate, F, T, P, ...
    lowpass_specmod, lowpass_tempmod, highpass_specmod, highpass_tempmod, ...
    complex_filters, separable, causal, spec_BW, temp_BW, spec_wavelet)

% Returns 2D transfer function for a spectrotemporal filter
% 
% -- Example --
% 
% temp_mod_rate = 10;
% temp_mod_lowpass = 0;
% spec_mod_rate = 1;
% spec_mod_lowpass = 0;
% 
% % number of samples
% % temporal sampling rate in Hz
% n_sec = 1;
% P.env_sr = 400;
% n_temporal_smps = P.env_sr*n_sec; % number of samples
% 
% % number of samples
% % spectral sampling rate in Hz
% n_octaves = 5;
% P.logf_spacing = 1/24;
% spectral_sr_cy_per_oct = 1/P.logf_spacing;
% n_spectral_smps = spectral_sr_cy_per_oct*n_octaves; % number of samples
% 
% % spectral filter transfer function
% complex_filters = false;
% separable = true;
% causal = false;
% spec_mod_highpass = false; 
% temp_mod_highpass = false;
% spectemp_filt_tf = filt_spectemp_mod(...
%     spec_mod_rate, temp_mod_rate,...
%     n_spectral_smps, n_temporal_smps, P, spec_mod_lowpass, temp_mod_lowpass, ...
%     spec_mod_highpass, temp_mod_highpass, complex_filters, separable, causal);
% 
% % impulse response
% spectemp_filt_irf = ifft2(spectemp_filt_tf)';
% spectemp_filt_irf = circshift(spectemp_filt_irf, [n_spectral_smps/2-1, 0]);
% 
% % corresponding frequencies for the spectral filter
% freqs = ...
%     [(0:1:n_spectral_smps/2), (-(n_spectral_smps/2-1):1:-1)]' ...
%     /spectral_sr_cy_per_oct;
% freqs = circshift(freqs, n_spectral_smps/2-1);
% 
% % plot image
% figure;
% imagesc(real(spectemp_filt_irf), ...
%     max(abs(real(spectemp_filt_irf(:))))*[-1 1]);
% 
% % modify color map
% colormap(cbrewer('div','RdBu',256));
% 
% % time axis
% set(gca,'XTick',[0 n_temporal_smps/2 n_temporal_smps],...
%     'XTickLabel',[0 n_temporal_smps/2 n_temporal_smps]/P.env_sr);
% xlabel('Time (seconds)');
% 
% % frequency axis
% freqs_to_plot = -floor(n_octaves/2):1:floor(n_octaves/2);
% yticks = nan(size(freqs_to_plot));
% for k = 1:length(freqs_to_plot)
%     [~,yticks(k)] = min(abs(freqs - freqs_to_plot(k)));
% end
% set(gca, 'YTick', yticks, ...
%     'YTickLabel', freqs_to_plot);
% ylabel('Frequency (Octaves)');
% 
% % set font size
% set(gca, 'FontName', 'Helvetica','FontSize',24)
% 
% 2017-06-19: Made it possible to make the filters separable
% 
% 2017-06-19: Added an example
% 
% 2017-06-28: Added non-causal filters

% add directory with useful 2D FT scripts
if ~exist('fft_freqs_from_siglen.m', 'file')
    directory_containing_this_file = fileparts(which(mfilename));
    addpath(genpath([directory_containing_this_file '/2DFT']));
end

if nargin < 6 || isempty(lowpass_specmod)
    lowpass_specmod = 0;
end

if nargin < 7 || isempty(lowpass_tempmod)
    lowpass_tempmod = 0;
end

if nargin < 8 || isempty(highpass_specmod)
    highpass_specmod = 0;
end

if nargin < 9 || isempty(highpass_tempmod)
    highpass_tempmod = 0;
end

if nargin < 10 || isempty(complex_filters)
    complex_filters = false;
end

if nargin < 11 || isempty(separable)
    separable = false;
end

if isnan(temp_mod_rate) || isnan(spec_mod_rate) || lowpass_specmod || lowpass_tempmod
    separable = true;
end

if nargin < 12 || isempty(causal)
    causal = true;
end

if nargin < 13 || isempty(spec_BW)
    spec_BW = 1;
end

if nargin < 14 || isempty(temp_BW)
    temp_BW = 1;
end

if nargin < 15 || isempty(spec_wavelet)
    spec_wavelet = 'mexicanhat';
end

if ~isnan(temp_mod_rate)

    % TF of temporal modulation filter
    Ht = filt_temp_mod(...
        abs(temp_mod_rate), T, P.env_sr, ...
        lowpass_tempmod, highpass_tempmod, causal, temp_BW);
        
else
    
    Ht = ones(T,1);
    separable = true;
    
end

if ~isnan(spec_mod_rate)
        
    % TF of spectral modulation filter
    Hs = filt_spec_mod(...
        spec_mod_rate, F, (1/P.logf_spacing),...
        lowpass_specmod, highpass_specmod, spec_BW, spec_wavelet);
        
else
    
    Hs = ones(F,1);
    separable = true;
    
end

% 2D transfer function
Hts = Ht * transpose(Hs);

% if the filter is a spectrotemporal filter
% orient the filter by removing quadrants
if ~separable && ~isnan(temp_mod_rate) && ~isnan(spec_mod_rate)
    
    % FFT frequencies excluding DC and nyq
    [f_spec, nyq_index] = fft_freqs_from_siglen(F,1);
    f_spec([1,nyq_index]) = NaN;
    [f_temp, nyq_index] = fft_freqs_from_siglen(T,1);
    f_temp([1,nyq_index]) = NaN;
    clear nyq_index
    
    % expand to matrices
    % T x F
    f_spec = ones(T,1) * f_spec';
    f_temp = f_temp * ones(1,F);
    
    % quadrants to zero out
    first_quad_to_zero = ...
        sign(f_temp) == -sign(temp_mod_rate) & sign(f_spec)==1;
    second_quad_to_zero = ...
        sign(f_temp) == sign(temp_mod_rate) & sign(f_spec)==-1;
    Hts(first_quad_to_zero | second_quad_to_zero) = 0;
    
end

% create complex-valued filters via analytic signal
if complex_filters
    temp_mod_bandpass = ~isnan(temp_mod_rate) && ~isnan(lowpass_tempmod) && ~lowpass_tempmod;
    spec_mod_bandpass = ~isnan(spec_mod_rate) && ~isnan(lowpass_specmod) && ~lowpass_specmod;
    
    if temp_mod_bandpass && ~spec_mod_bandpass % temporal
        Hts = analytic_from_spectrum_2D(Hts, 1);
        
    elseif spec_mod_bandpass && ~temp_mod_bandpass % spectral
        Hts = analytic_from_spectrum_2D(Hts, 2);
        
    elseif temp_mod_bandpass && spec_mod_bandpass % spectrotemporal
        Hts = analytic_from_spectrum_2D(Hts, 1); % direction doesn't matter for this
        
    else
        assert(~temp_mod_bandpass && ~spec_mod_bandpass);
    end
end
    
    
% FlushEvents('keyDown');
% X = ifft2(Hts);
% fprintf('%.3f Hz, %.3f cyc/oct, %.2f\n',...
%     temp_mod_rate, spec_mod_rate, ...
%     log10(std(real(X(:)))) - log10(std(imag(X(:)))))

% plot_2DFT(Hts', [1/P.logf_spacing, P.env_sr]);
% set(gcf, 'Position', [200 200 1000 600]);

% GetChar;
% close all;

