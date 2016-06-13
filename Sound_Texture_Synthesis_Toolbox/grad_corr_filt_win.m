% GRAD = GRAD_CORR_FILT_WIN(VECTOR1, VECTOR2, CURRENT_CBC, MOD_FILT, USE_ZP, WIN)
%
% computes gradient of correlations between VECTOR1 and VECTORS2 after they
% have been filtered MOD_FILT, as computed by the function stat_corr_win
%
% if USE_ZP is 1, vectors are zero-padded prior to filtering
%
% the correlation is weighted by an optional window WIN
%
% gradient is computed with respect to VECTOR1
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

function grad = grad_corr_filt_win(vector1, vector2, current_cbc, mod_filt, use_zp, win)

if nargin==5
    win = ones(length(vector1),1);
end
win=win/sum(win);

f_env1 = apply_filter(mod_filt,vector1,use_zp);
f_env2 = apply_filter(mod_filt,vector2,use_zp);

env_std1 = sqrt(sum(win.*(f_env1.^2)));
env_std2 = sqrt(sum(win.*(f_env2.^2)));

grad = apply_filter(mod_filt,win.*f_env2,use_zp)/(env_std2*env_std1) - current_cbc/(env_std1^2)*apply_filter(mod_filt,win.*f_env1,use_zp);
