function coch_synth_matched = match_filtcoch(...
    coch_orig, coch_synth, P, ti, match_mean_and_var)

% Matches the subbands of filtered cochleograms. Much simpler than version 1 and
% doesn't depend on the nsl toolbox.
%
% -- Example: Test Reconstruction --
%
% P = synthesis_parameters_default;
% 
% % add spectrotemporal filters
% P.temp_mod_to_match = [];
% P.spec_mod_to_match = [];
% if P.match_spectemp_mod
%     for i = 1:length(P.temp_mod_rates)
%         for j = 1:length(P.spec_mod_rates);
%             if P.temp_mod_rates(i) == 0
%                 P.temp_mod_to_match = ...
%                     [P.temp_mod_to_match, 0];
%                 P.spec_mod_to_match = ...
%                     [P.spec_mod_to_match, P.spec_mod_rates(j)];
%             else
%                 P.temp_mod_to_match = ...
%                     [P.temp_mod_to_match, ...
%                     P.temp_mod_rates(i), -P.temp_mod_rates(i)];
%                 P.spec_mod_to_match = ...
%                     [P.spec_mod_to_match, ...
%                     P.spec_mod_rates(j) * ones(1,2)];
%             end
%         end
%     end
% end
% 
% % test reconstruction (i.e. match cochleogram to itself)
% coch = randn(P.env_sr*4, 9*round(1/P.logf_spacing));
% match_mean_and_var = true;
% coch_recon = ...
%     match_filtcoch(coch, coch, P, 1:P.env_sr*4, match_mean_and_var);
% figure;
% plot(coch, coch_recon);

% add directory with useful 2D FT scripts
if ~exist('fft_freqs_from_siglen.m', 'file')
    directory_containing_this_file = fileparts(which(mfilename));
    addpath(genpath([directory_containing_this_file '/2DFT']));
end

% check matched sizes
assert(all(size(coch_orig) == size(coch_synth)));

% dimensions
% time x frequency
[T,F] = size(coch_orig);

% 2D fourier transforms
FT_coch_orig = fft2(coch_orig);
FT_coch_synth = fft2(coch_synth);

% total number of filters
assert(length(P.spec_mod_to_match) == length(P.temp_mod_to_match));
n_filters = length(P.spec_mod_to_match);

% accumulated FT of filtered cochleograms for matched synthetics
% transfer functions are also accumulated
accum_FT_filtcoch_synth_matched = zeros(T,F);
accum_Hts = zeros(T,F);
for i = 1:n_filters
    
    % filter transfer function
    Hts = filt_spectemp_mod(...
        P.spec_mod_to_match(i), P.temp_mod_to_match(i), F, T, P);
    
    % apply filter and revert to cochleogram domain
    filtcoch_orig = real(ifft2(FT_coch_orig .* Hts));
    filtcoch_synth = real(ifft2(FT_coch_synth .* Hts));
    
    % match mean/variance or full histogram for subset of temporal indices
    if match_mean_and_var
        filtcoch_synth = ...
            match_coch_mean_and_var(filtcoch_orig, filtcoch_synth, ti);
    else
        filtcoch_synth = ...
            match_coch_hists(filtcoch_orig, filtcoch_synth, ti);
    end
    
    % accumulate FT of matched cochleograms
    accum_FT_filtcoch_synth_matched = ...
        accum_FT_filtcoch_synth_matched ...
        + fft2(filtcoch_synth) .* conj(Hts);
    
    % accumulate FT of transfer functions
    accum_Hts = accum_Hts + Hts .* conj(Hts);
    
end

if any(accum_Hts(:) < 1e-10)
    error('Accumulated transfer function has entries with close-to-zero values');
end

% divide by accumulated transfer functions
FT_coch_synth_matched = accum_FT_filtcoch_synth_matched ./ accum_Hts;

% invert back to the cochleogram domain
coch_synth_matched = real(ifft2(FT_coch_synth_matched));