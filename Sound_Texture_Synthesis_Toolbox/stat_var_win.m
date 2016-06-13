% W_VAR = STAT_VAR_WIN(S, WIN)
%
% returns variance of S weighted by WIN, a vector the same length as S.
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

function w_var = stat_var_win(s, win)

if nargin==1
    win = ones(length(s),1);
end
win=win/sum(win);

w_var = sum(win.*(s-sum(win.*s)).^2);