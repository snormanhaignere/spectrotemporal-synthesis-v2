% MP_GRAD = GRAD_MOD_POWER_WIN(S, MOD_FILT, USE_ZP, WIN)
%
% returns gradient of power in modulation band of an envelope S normalized
% by the envelope variance, as computed by the function stat_mod_power
% using MOD_FILT as the filter
%
% if USE_ZP is 1, vectors are zero-padded prior to filtering
%
% the power averaging is weighted by an optional window WIN
%
% gradient is with respect to S
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

function  mp_grad = grad_mod_power_win(s, mod_filt, use_zp, win)

if nargin==3
    win = ones(length(s),1);
end
win=win/sum(win);


v = stat_var_win(s,win);
mp_unnorm = sum(win.*(apply_filter(mod_filt,s,use_zp).^2));
fwfe = apply_filter(mod_filt,win.*apply_filter(mod_filt,s,use_zp),use_zp);
me = s - sum(win.*s);
wme = win.*me;
mwme = wme;

mp_grad = (2*v*fwfe - 2*mp_unnorm*mwme)/(v^2);
