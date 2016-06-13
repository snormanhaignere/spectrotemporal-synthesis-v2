% AC_GRADS = GRAD_ENV_AC_SCALED_WIN(ENV, FILTS, SAMPLE_SPACING, USE_ZP, WIN)
%
% computes gradient of auto-correlation coefficient of envelope ENV at lags
% specified in SAMPLE_SPACNG (in samples)
%
% gradients for each lag are returned as column vectors of AC_GRADS
%
% the autocorrelation is that computed in stat_env_ac_scaled_win - there is
% a filter (a column vector of FILTS) for each lag that low pass filters in
% proportion to the lag.
%
% the other arguments should be the same as whatever was passed to
% stat_env_ac_scaled_win for calculation of the autocorrelation statistic
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


function ac_grads = grad_env_ac_scaled_win(env, filts, sample_spacing, use_zp, win)

if nargin==4
    win = ones(length(env),1);
end
win=win/sum(win);
L=length(env);
num_sub = size(filts,2);

ac_grads=[];
for p=1:length(sample_spacing)
    num_samp = sample_spacing(p);
    
    f_env = apply_filter(filts(:,p),env,use_zp);
    meanf_env = mean(f_env);
    mf_env = f_env-meanf_env;
    
    env_var = mean(mf_env.^2);
    
    temp = shift_s(win.*shift_s(mf_env,num_samp),num_samp) + shift_s(win.*shift_s(mf_env,-num_samp),-num_samp);
    %temp2 = win.*(temp-mean(temp));
    temp2 = (temp-mean(temp));
    
    temp3 = sum(win.*(shift_s(mf_env,-num_samp).*shift_s(mf_env,num_samp)));
    
    grad = (env_var*apply_filter(filts(:,p),temp2,use_zp) - ...
        2/L*temp3*apply_filter(filts(:,p), mf_env,use_zp))/(env_var^2);
    
    ac_grads = [ac_grads grad];
end

