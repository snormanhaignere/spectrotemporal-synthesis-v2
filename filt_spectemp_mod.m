function Hts = filt_spectemp_mod(...
    spec_mod_rate, temp_mod_rate, F, T, P, ...
    lowpass_specmod, lowpass_tempmod, highpass_specmod, highpass_tempmod, ...
    complex_filters, )

% Returns 2D transfer function for a spectrotemporal filter

% add directory with useful 2D FT scripts
if ~exist('fft_freqs_from_siglen.m', 'file')
    directory_containing_this_file = fileparts(which(mfilename));
    addpath(genpath([directory_containing_this_file '/2DFT']));
end

if nargin < 6 
    lowpass_specmod = 0;
    lowpass_tempmod = 0;
end

if nargin < 8
    highpass_specmod = 0;
    highpass_tempmod = 0;
end

if nargin < 10
    complex_filters = false;
end

if ~isnan(temp_mod_rate)

    % TF of temporal modulation filter
    Ht = filt_temp_mod(...
        abs(temp_mod_rate), T, P.env_sr, ...
        lowpass_tempmod, highpass_tempmod);
        
else
    
    Ht = ones(T,1);
    
end

if ~isnan(spec_mod_rate)
        
    % TF of spectral modulation filter
    Hs = filt_spec_mod(...
        spec_mod_rate, F, (1/P.logf_spacing),...
        lowpass_specmod, highpass_specmod);
        
else
    
    Hs = ones(F,1);
    
end

% 2D transfer function
Hts = Ht * transpose(Hs);

% if the filter is a spectrotemporal filter
% orient the filter by removing quadrants
if ~isnan(temp_mod_rate) && ~isnan(spec_mod_rate)
    
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
    if ~isnan(temp_mod_rate) && isnan(spec_mod_rate) % temporal
        Hts = analytic_from_spectrum_2D(Hts, 1);
        
    elseif ~isnan(spec_mod_rate) && isnan(temp_mod_rate) % spectral
        Hts = analytic_from_spectrum_2D(Hts, 2);
        
    elseif ~isnan(spec_mod_rate) && ~isnan(temp_mod_rate) % spectrotemporal
        Hts = analytic_from_spectrum_2D(Hts, 1); % direction doesn't matter for this
        
    else
        assert(isnan(temp_mod_rate) && isnan(spec_mod_rate));
        
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

