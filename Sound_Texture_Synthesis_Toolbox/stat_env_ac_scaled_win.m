% AC_VALUES = STAT_ENV_AC_SCALED_WIN(ENV, FILTS, SAMPLE_SPACING, USE_ZP, WIN)
%
% computes auto-correlation coefficient values (subtracts mean and divides
% out variance) of ENV at lags specified in SAMPLE_SPACING (in samples)
%
% in this version (10/10/09) the lags to be imposed are spaced 
% logarithmically, and there is a filter (a column vector of FILTS)
% for each lag that low pass filters in proportion to the lag. 
%
% the optional variable WIN is used to compute the correlation as a
% weighted average (typically with a window that fades to zero at the
% edges)
%
% for now, the global mean and variance are used to normalize the
% correlation coefficient (didn't get around to coding up the windowed mean
% and variance in the gradient calculation)
%
% if USE_ZP is set to 1, ENV is zero-padded before filtering (filters must be
% twice as long in order for this to work)
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

function ac_values = stat_env_ac_scaled_win(env, filts, sample_spacing, use_zp, win)

if nargin==4
    win = ones(length(env),1);
end
win = win/sum(win);

for p=1:length(sample_spacing)
    num_samp = sample_spacing(p);
    
    f_env = apply_filter(filts(:,p),env,use_zp);
    meanf_env = mean(f_env);
    mf_env = f_env-meanf_env;
    
    env_var = mean(mf_env.^2);
    
    ac_values(p) = sum(win.*(shift_s(mf_env,-num_samp).*shift_s(mf_env,num_samp)))/env_var;
end

