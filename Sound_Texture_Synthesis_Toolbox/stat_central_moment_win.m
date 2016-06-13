% M = STAT_CENTRAL_MOMENT_WIN(X, N, WIN, X_MEAN)
%
% computes central moment of a vector X for use in texture synthesis
%
% the moment to be computed is specified by the variable N
%
% the moment is a weighted average with weights specified by WIN, a
% vector of the same length as X
%
% X_MEAN can be optionally provided if it has been computed elsewhere
%
% **Note - the 2nd moment is not the variance, but rather the standard
% deviation divided by the mean. This was to make it unitless like the skew
% and kurtosis.
%
% assumes X_MEAN is the windowed mean
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

function m = stat_central_moment_win(x,n,win,x_mean)

win = win/sum(win);
if nargin==3
    x_mean = sum(win.*x);
end
if n==1
    m = x_mean;
elseif n==2 %variance
    m = sum(win.*((x-x_mean).^2));
    m = sqrt(m)/x_mean;
elseif n==3 %skew
    m2 = sum(win.*((x-x_mean).^2));
    m = sum(win.*((x-x_mean).^3))/(m2^(3/2));
    %m = sum(win.*((x-x_mean).^3))/(central_moment_win(x,2,win,x_mean)^(3/2));
elseif n==4 %kurtosis
    m2 = sum(win.*((x-x_mean).^2));
    m = sum(win.*((x-x_mean).^4))/(m2^2);
    %m = sum(win.*((x-x_mean).^4))/(central_moment_win(x,2,win,x_mean)^2);
end