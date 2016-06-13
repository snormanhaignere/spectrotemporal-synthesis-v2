% GRADS = GRAD_CORR_WIN(VECTOR1, VECTORS2, CURRENT_CBC, WIN)
%
% computes gradient of correlations between VECTOR1 and VECTORS2, as
% computed by the function stat_corr_win
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

function grads = grad_corr_win(vector1, vectors2, current_cbc, win)

if nargin==4
    win = ones(length(vector1),1);
end
win=win/sum(win);
L=length(vector1);

f_env1 = vector1;
f_env2 = vectors2;

meanf_env1 = mean(f_env1);
mf_env1 = f_env1-meanf_env1;
env_std1 = sqrt(mean(mf_env1.^2));

meanf_env2 = mean(f_env2);
mf_env2 = f_env2-meanf_env2;
env_std2 = sqrt(mean(mf_env2.^2));

%compute additional stuff for gradient
smf_env2 = mf_env2;
swsmf_env2 = win.*smf_env2;
mean_swsmf_env2 = mean(swsmf_env2);
mswsmf_env2 = swsmf_env2-mean_swsmf_env2;
fmswsmf_env2 = mswsmf_env2;
    
fmmf_env1 = mf_env1; %mean is already zero, so no need to subtract again

grads = fmswsmf_env2./(ones(L,1)*env_std1*env_std2) - fmmf_env1/(L*env_std1^2)*current_cbc;


