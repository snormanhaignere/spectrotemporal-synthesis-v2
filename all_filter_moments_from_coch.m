function M = all_filter_moments_from_coch(coch, P, ti, envelopes)

% Calculates moments of cochlear and modulation filters from an input
% cochleogram. 
% 
% -- Inputs --
% 
% coch: cochleogram (e.g. from wav2coch)
% 
% P: parameter structure (see default_synthesis_parameters)
% 
% ti: indices of the time-points over which to average when measuring the
% moments, useful for excluding boundary time points; e.g. for computing
% cochlear marginals the moments are computed as moment(coch(ti,:))
% 
% -- Example --
% 
% P = toy_synthesis_parameters;
% coch = randn(P.env_sr*4, 9*round(1/P.logf_spacing));
% ti = 1:size(coch,1);
% M = all_filter_moments_from_coch(coch, P, ti)
% 
% 2016-12-09: Optional choice to compute moments of filter envelopes

if nargin < 3
    ti = 1:size(coch,1);
end

if nargin < 4
    envelopes = false;
end

% cochleogram moments
M.coch_env = moment_measures(coch(ti,:),1);

% 2D fourier transform of the padded cochleogram
FT_padded_coch = fft2(pad_coch(coch,P));

% remove zeros if present
% don't need to measure the temporal statistics of a signal with no temporal
% modulation
pos_temp_mod_rates = P.temp_mod_rates(P.temp_mod_rates>0);

% temporal modulation moments
M.temp_mod = filtcoch_moments(FT_padded_coch, NaN, pos_temp_mod_rates, P, ti, envelopes);

% spectral modulation moments
M.spec_mod = filtcoch_moments(FT_padded_coch, P.spec_mod_rates, NaN, P, ti, envelopes);

% spectrotemporal modulation moments
temp_mod_rates_pos_and_neg = [pos_temp_mod_rates, -pos_temp_mod_rates];
M.spectemp_mod = ...
    filtcoch_moments(FT_padded_coch, P.spec_mod_rates, temp_mod_rates_pos_and_neg, P, ti, envelopes);

function filter_moments = ...
    filtcoch_moments(FT_padded_coch, spec_mod_rates, temp_mod_rates, P, ti, envelopes)

% amount of frequency padding
n_freq_pad_smps = round(P.freq_pad_oct / P.logf_spacing);

% dimensions of frequency-padded cochleogram
[T,F] = size(FT_padded_coch);

n_moments = size(moment_measures(1,1),2);
n_spec_mod_rates = length(spec_mod_rates);
n_temp_mod_rates = length(temp_mod_rates);
filter_moments = nan(n_spec_mod_rates, n_temp_mod_rates, F-n_freq_pad_smps, n_moments);
for i = 1:n_spec_mod_rates
    for j = 1:n_temp_mod_rates
        
        % transfer function of spectrotemporal filter
        Hts = filt_spectemp_mod(...
            spec_mod_rates(i), temp_mod_rates(j), F, T, P,...
            [], [], [], [], [], [], [], ...
            P.spec_BW, P.temp_BW, P.spec_wavelet);
        
        % whether or not to compute envelopes of the signal
        if envelopes
            
            if ~isnan(temp_mod_rates(j)) % compute temporal envelopes if possible
                Hts = analytic_from_spectrum_2D(Hts, 1);
            elseif ~isnan(spec_mod_rates(i)) % otherwise compute spectral envelopes
                Hts = analytic_from_spectrum_2D(Hts, 2);
            else
                error('Envelopes cannot be computed for filter without any modulation');
            end
            
            % apply transfer function, absolute value of analytic signal
            filtcoch_padded = abs(ifft2(FT_padded_coch .* Hts));
            filtcoch = remove_pad(filtcoch_padded, P);
            
        else
            
            % apply transfer function
            filtcoch_padded = real(ifft2(FT_padded_coch .* Hts));
            filtcoch = remove_pad(filtcoch_padded, P);
        
        end
        
        % imagesc(filtcoch);
        % plot_2DFT(Hts', fliplr([P.env_sr, 1/P.logf_spacing]));
                       
        % moments of filtered cochleogram
        filter_moments(i,j,:,:) = moment_measures(filtcoch(ti,:),1);
        
    end
end

% remove singular dimensions
if n_spec_mod_rates == 1
    filter_moments = squeeze_dim(filter_moments, 1);
end
if n_temp_mod_rates == 1
    filter_moments = squeeze_dim(filter_moments, 2);
end





