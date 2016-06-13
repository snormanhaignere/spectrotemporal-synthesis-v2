% CBC_VALUE = STAT_CORR_WIN(VECTOR1, VECTORS2, WIN)
%
% computes correlation coefficients between one subband envelope (VECTOR1)
% and some others (contained in VECTORS2)
%
% the input format facilitates gradient checks.
%
% computes correlation as windowed average with weights given by WIN
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

function cbc_value = stat_corr_win(vector1, vectors2, win)

if nargin==4
    win = ones(length(vector1),1);
end
win=win/sum(win);

f_env1 = vector1;
f_env2 = vectors2;

meanf_env1 = mean(f_env1);
mf_env1 = f_env1-meanf_env1;
env_std1 = sqrt(mean(mf_env1.^2));

meanf_env2 = mean(f_env2);
mf_env2 = f_env2-meanf_env2;
env_std2 = sqrt(mean(mf_env2.^2));

cbc_value = sum(win.*mf_env1.*mf_env2)/(env_std1*env_std2);

