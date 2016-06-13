% C2 = STAT_MOD_C2_WIN(S, FILTS, USE_ZP, WIN)
%
% returns phase-adjusted "C2" correlations of analytic subbands produced by
% filtering S with the filters in the columns of FILTS
%
% if USE_ZP is 1, S is zero-padded prior to filtering
%
% the correlation is weighted by an optional window WIN
%
% the correlation is complex-valued, and the real and complex parts are
% returned in the two columns of C2
%
% typically S is a subband envelope and FILTS are modulation filters.
% Correlations are computed between adjacent modulation bands of S. The
% filters must be octave spaced, such that the center frequency of adjacent
% filters differ by a factor of 2.
%
% can be passed the pre-filtered subbands in S
%
% The formulation is derived from the relative phase statistic introduced
% by Portilla and Simoncelli for visual image analysis and synthesis. This
% version differs from theirs in using a proper correlation coefficient (by
% dividing out the standard deviation of the modulation bands), so as to
% make it independent of the modulation power
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

function C2 = stat_mod_C2_win(s, filts, use_zp, win)

if nargin==3
    win = ones(length(s),1);
end
win=win/sum(win);
L=length(s);

if min(size(s))==1 %have to generate subbands
    if use_zp
        subbands = generate_subbands_zp(s, filts);
    else
        subbands = generate_subbands(s, filts);
    end
else
    subbands=s;
end

%make subbands analytic
analytic_subbands = hilbert(subbands);

N = size(analytic_subbands,2);
o=1;
for k = 1:N-o
    c = (analytic_subbands(:,k).^2)./abs(analytic_subbands(:,k)); %phase adjusted coarse subband, normalized
    sig_cw = sqrt(sum(win.*(real(c).^2)));
    sig_fw = sqrt(sum(win.*(real(analytic_subbands(:,k+o)).^2)));
    
    C2(k,1) = sum(win.*real(c).*real(analytic_subbands(:,k+o)))/(sig_cw*sig_fw);
    C2(k,2) = sum(win.*real(c).*imag(analytic_subbands(:,k+o)))/(sig_cw*sig_fw);    
end

