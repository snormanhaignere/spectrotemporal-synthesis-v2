function C = moment_comparisons(M_orig, M_synth, iteration, C)

% C = moment_comparisons(M_orig, M_synth, iteration, C)
% 
% Compares moments using two metrics (SNR and correlation). 

if nargin < 3
    iteration = 1;
end

if nargin < 4;
    C = struct;
end

all_fields = fieldnames(M_orig);
for i = 1:length(all_fields);
    
    % SNR: abs(orig) / abs(syn - orig)
    syn = M_synth.(all_fields{i});
    orig = M_orig.(all_fields{i});
    C.([all_fields{i} '_SNR']) = 20*log10(abs(orig)) - 20*log10(abs(syn-orig));
    
    % unwrap SNR values across all filters
    % dim1 x dim2 ... x moment
    % -> dim1 * dim2 ... x moment
    SNR_unwrapped = unwrap_allbutlast( C.([all_fields{i} '_SNR']) );
    
    % median SNR across all filters
    C.([all_fields{i} '_SNR_median'])( iteration, : ) = nanmedian( SNR_unwrapped );
    
    % unwrap filter values
    % dim1 x dim2 ... x moment
    % -> dim1 * dim2 ... x moment
    syn_unwrapped = unwrap_allbutlast( syn );
    orig_unwrapped = unwrap_allbutlast( orig );
    
    % correlate values across all filters
    R = nanfastcorr(syn_unwrapped, orig_unwrapped);
    
    % signed R2
    R2 = sign(R) .* R.^2;
    C.([all_fields{i} '_corr'])( iteration, : ) = R2;
        
    
end
clear X Y;

%% Old code
% % frequency-specific moment correlation
% dims = size(C.spectemp_orig_moments);
% assert_equal_dimensions(C.spectemp_orig_moments, C.spectemp_synth_moments);
% syn = reshape(C.spectemp_orig_moments, [prod(dims(1:2)), prod(dims(3:4))]);
% orig = reshape(C.spectemp_synth_moments, [prod(dims(1:2)), prod(dims(3:4))]);
% C.spectemp_moment_corr_freqspec = reshape( fastcorr(syn, orig), dims(3:4) );
% 
% %%
% 
% 
% 
% % SNRs
% C.coch_moment_SNR = ...
%     20*log10(abs(M_synth.coch - M_orig.coch)) ...
%     - 20*log10(abs(M_orig.coch));
% 
% C.temp_moment_SNR = ...
%     20*log10(abs(M_synth.temp - M_orig.temp)) ...
%     - 20*log10(abs(M_orig.temp));
% 
% C.spectemp_moment_SNR = ...
%     20*log10(abs(M_synth.spectemp - M_orig.spectemp)) ...
%     - 20*log10(abs(M_orig.spectemp));
% 
% % average SNRs
% C.coch_moment_SNR_mean( iteration, : ) = ...
%     mean( unwrap_allbutlast( C.coch_moment_SNR ) );
% 
% C.temp_moment_SNR_mean( iteration,: ) = ...
%     mean( unwrap_allbutlast( C.temp_moment_SNR ) );
% 
% C.spectemp_moment_SNR_mean( iteration,: ) = ...
%     mean( unwrap_allbutlast( C.spectemp_moment_SNR ) );
% 
% % moment correlations across all features
% C.coch_moment_corr( iteration, : ) = ...
%     unwrap_and_corr(C.coch_orig_moments, C.coch_synth_moments);
% 
% C.temp_moment_corr( iteration, : ) = ...
%     unwrap_and_corr(C.temp_orig_moments, C.temp_synth_moments);
% 
% C.spectemp_moment_corr( iteration, : ) = ...
%     unwrap_and_corr(C.spectemp_orig_moments, C.spectemp_synth_moments);
% 
% % frequency-specific moment correlation
% dims = size(C.spectemp_orig_moments);
% assert_equal_dimensions(C.spectemp_orig_moments, C.spectemp_synth_moments);
% syn = reshape(C.spectemp_orig_moments, [prod(dims(1:2)), prod(dims(3:4))]);
% orig = reshape(C.spectemp_synth_moments, [prod(dims(1:2)), prod(dims(3:4))]);
% C.spectemp_moment_corr_freqspec = reshape( fastcorr(syn, orig), dims(3:4) );
% 
% % average across frequency-specific correlation
% C.spectemp_moment_corr_freqspec_mean( iteration, : ) = ...
%     mean(C.spectemp_moment_corr_freqspec, 1);



