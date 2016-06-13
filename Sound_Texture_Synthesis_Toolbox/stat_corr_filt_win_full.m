% CBC_VALUE = STAT_CORR_FILT_WIN_FULL(VECTORS, MOD_FILTS, USE_ZP, WIN)
%
% computes correlation coefficients between a set of vectors (columns of
% VECTORS) after they have been filtered by the columns of MOD_FILTS
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

function cbc_value = stat_corr_filt_win_full(vectors, mod_filts, use_zp, win)

if nargin==3
    win = ones(size(vectors,1),1);
end
win=win/sum(win);

for k=1:size(mod_filts,2)
    f_envs = apply_filter_mult(mod_filts(:,k),vectors,use_zp);
    
    meanf_envs = mean(f_envs);
    mf_envs = f_envs-ones(size(f_envs,1),1)*meanf_envs;
    env_stds = sqrt(mean(mf_envs.^2));
        
    cbc_value(:,:,k) = ( ( (win*ones(1,size(f_envs,2))) .*mf_envs)' *mf_envs)./(env_stds'*env_stds);
end
