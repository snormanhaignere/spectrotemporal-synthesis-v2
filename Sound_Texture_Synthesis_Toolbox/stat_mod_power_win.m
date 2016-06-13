% MP = STAT_MOD_POWER_WIN(S, MOD_FILTS, USE_ZP, WIN)
%
% Returns vector of power in modulation bands of an envelope S normalized
% by the envelope variance
%
% Modulation bands are computed using the filters in the columns of
% MOD_FILTS.
%
% If bands have been pre-computed, can be passed as S.
%
% If USE_ZP is 1, vectors are zero-padded prior to filtering.
%
% The power averaging is weighted by an optional window WIN.
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

function  mp = stat_mod_power_win(s, mod_filts, use_zp, win)

if nargin==3
    win = ones(length(s),1);
end
win=win/sum(win);

if min(size(s))==1 %have to generate mod subbands
    if use_zp
        mod_subbands = generate_subbands_zp(s, mod_filts);
    else
        mod_subbands = generate_subbands(s, mod_filts);
    end
else
    mod_subbands=s;
end
s_var = stat_var_win(s,win);

mp = sum((win*ones(1,size(mod_subbands,2))).*(mod_subbands.^2))/s_var;
