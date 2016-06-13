% G = GRAD_CENTRAL_MOMENT_WIN(X, N, WIN, X_MEAN)
% 
% returns gradient of Nth central moment of X computed by
% stat_central_moment_win with weighting window WIN
% 
% gradient is with respect to X
%
% Note: second moment is not the variance, but rather the standard
% deviation divided by the mean (to produce a unitless statistic analogous
% to skew and kurtosis)
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

function g = grad_central_moment_win(x, n, win, x_mean)
win=win/sum(win);
if nargin==3
    x_mean = sum(win.*x);
end
if n==3 | n==4
    m2 = sum(win.*((x-x_mean).^2));
    m3 = sum(win.*((x-x_mean).^3));
end
if n==1
    g = win;
elseif n==2
    %g = 2*win.*(1-win).*(x-x_mean);
    m2 = sum(win.*((x-x_mean).^2));
    g = win.*(1-win).*(x-x_mean)/(x_mean*sqrt(m2)) - win*sqrt(m2)/(x_mean^2);
elseif n==3
    g = 3/(m2^1.5)*win.*((x-x_mean).^2 - m2 - m3/m2*(x-x_mean));
elseif n==4
    m4 = sum(win.*((x-x_mean).^4));
    g = 4*win.*( (x-x_mean).^3 - m3 - m4/m2*(x-x_mean))/(m2^2);
end