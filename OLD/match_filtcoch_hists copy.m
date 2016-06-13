function coch_synth_matched = match_filtcoch_hists(coch_orig, coch_synth, P, ti)

% Matches the subbands of filtered cochleograms. Much simpler than version 1 and
% doesn't depend on the nsl toolbox.
% 
% -- Example: Test Reconstruction --
% 
% P = toy_synthesis_parameters_v2;
% P.temp_mod_to_match = [P.temp_mod_rates, -P.temp_mod_rates(P.temp_mod_rates>0)];
% P.spec_mod_to_match = [P.spec_mod_rates];
% 
% % test reconstruction (i.e. match cochleogram to itself)
% coch = randn(P.env_sr*4, 9*round(1/P.logf_spacing));
% coch_recon = ...
%     match_filtcoch_hists_v2(coch, coch, P, 1:P.env_sr*4);
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

% number of spectral and temporal modulation rates
n_spec_mod_rates = length(P.spec_mod_to_match);
n_temp_mod_rates = length(P.temp_mod_to_match);

% accumulated FT of filtered cochleograms for matched synthetics
% transfer functions are also accumulated
accum_FT_filtcoch_synth_matched = zeros(T,F);
accum_Hts = zeros(T,F);
for i = 1:n_spec_mod_rates
    for j = 1:n_temp_mod_rates
        
        % no need to match negative temporal modulation filters
        % which are redudant with positive temporal modulation filters
        % NaN value for spectral modulation is a flag, indicating that the
        % filters are not modulated in frequency
        isa_negative_temp_mod_filt = ...
            P.temp_mod_to_match(j) < 0 && isnan(P.spec_mod_to_match(i));
        
        if ~isa_negative_temp_mod_filt
            
            % filter transfer function
            Hts = filt_spectemp_mod(...
                P.spec_mod_to_match(i), P.temp_mod_to_match(j), F, T, P);
            
            % apply filter and revert to cochleogram domain
            filtcoch_orig = real(ifft2(FT_coch_orig .* Hts));
            filtcoch_synth = real(ifft2(FT_coch_synth .* Hts));
            
            % histogram match, all temporal indices
            filtcoch_synth = ...
                match_coch_hists(filtcoch_orig, filtcoch_synth, 1:T);
            
            % histogram match, specified subset of temporal indices
            filtcoch_synth = ...
                match_coch_hists(filtcoch_orig, filtcoch_synth, ti);
            
            % accumulate FT of matched cochleograms
            accum_FT_filtcoch_synth_matched = ...
                accum_FT_filtcoch_synth_matched ...
                + fft2(filtcoch_synth) .* conj(Hts);
            
            % accumulate FT of transfer functions
            accum_Hts = accum_Hts + Hts .* conj(Hts);
            
        end
    end
end

if any(accum_Hts(:) < 1e-10)
    error('Accumulated transfer function has entries with close-to-zero values');
end

% divide by accumulated transfer functions
FT_coch_synth_matched = accum_FT_filtcoch_synth_matched ./ accum_Hts;

% invert back to the cochleogram domain
coch_synth_matched = real(ifft2(FT_coch_synth_matched));