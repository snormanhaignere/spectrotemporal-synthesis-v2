% CBC_VALUE = STAT_CORR_FILT_WIN(VECTOR1, VECTORS2, MOD_FILTS, USE_ZP, WIN)
%
% computes correlation coefficients between VECTOR1 and VECTORS2 after they
% have been filtered by the columns of MOD_FILTS 
%
% if USE_ZP is 1, vectors are zero-padded prior to filtering
%
% the correlation is weighted by an optional window WIN
%
% typically VECTORS is a set of subband envelopes, and MOD_FILTS are
% modulation filters
%
%

% This code is part of an instantiation of a sound texture synthesis
% algorithm developed with Eero Simoncelli and described in this paper:
%
% McDermott, J.H., Simoncelli, E.P. (2011) Sound texture perception via
% statistics of the auditory periphery: Evidence from sound synthesis.
% Neuron, 71, 926-940. 
%
% Dec 2012 -- Josh McDermott <jhm@mit.edu>

function cbc_value = stat_corr_filt_win(vector1, vectors2, mod_filts, use_zp, win)

if nargin==4
    win = ones(length(vector1),1);
end
win=win/sum(win);

for k=1:size(mod_filts,2)
    f_env1 = apply_filter(mod_filts(:,k),vector1,use_zp);
    f_envs2 = apply_filter_mult(mod_filts(:,k),vectors2,use_zp);
    
    env_std1 = sqrt(sum(win.*(f_env1.^2))); %assumes filters are bandpass, so zero-mean
    env_stds2 = sqrt(sum( (win*ones(1,size(f_envs2,2))).*(f_envs2.^2))); %assumes filters are bandpass, so zero-mean
        
    cbc_value(:,k) = ( f_envs2' * (win.*f_env1))./(env_stds2'*env_std1);
end
