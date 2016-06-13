% GRAD = GRAD_MOD_C2_WIN(S, FILTS, WHICH_PART, USE_ZP, WIN)
%
% returns gradient of either the real or imaginary part of the
% phase-adjusted "C2" correlation of analytic subbands produced by
% filtering S with the filters in the columns of FILTS, as computed by the
% function stat_mod_C2_win.
%
% WHICH_PART can be 1 (for real part) or 2 (for imaginary part)
%
% if USE_ZP is 1, S is zero-padded prior to filtering
%
% the correlation is weighted by an optional window WIN
%
% assumes it gets passed two filters in FILTS
%
% can be passed the pre-filtered subbands in S
%
% The gradient is only correct for octave-spaced filters, as it relies on the
% use of the double angle formula to simplify the angle in the complex plane
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

function grad = grad_mod_C2_win(s, filts, which_part, use_zp, win)

if nargin==4
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

sig_cw = sqrt(sum(win.*real(analytic_subbands(:,1)).^2));
sig_fw = sqrt(sum(win.*real(analytic_subbands(:,2)).^2));

c_orig = analytic_subbands(:,1);

u = 2*(real(c_orig).^2)./abs(c_orig) - abs(c_orig); %true because of double angle formula (only applies when f2/f1 = 2)

if which_part==1
    f = real(analytic_subbands(:,2));
elseif which_part==2
    f = imag(analytic_subbands(:,2));
end

top_right_first_pt = ...
    (apply_filter(filts(:,1), 4*real(c_orig)./abs(c_orig).*win.*f, use_zp) - ...
    apply_filter(filts(:,1), 2*(real(c_orig).^2)./ (abs(c_orig).^3).*real(c_orig).*win.*f, use_zp) + ...
    imag(hilbert(apply_filter(filts(:,1), 2*(real(c_orig).^2)./ (abs(c_orig).^3).*imag(c_orig).*win.*f, use_zp))) - ...
    apply_filter(filts(:,1), real(c_orig)./abs(c_orig).*win.*f, use_zp) + ... %should cancel with first part...
    imag(hilbert(apply_filter(filts(:,1), imag(c_orig)./abs(c_orig).*win.*f, use_zp))))*sig_cw*sig_fw; %sign changes because of hilbert

if which_part==1
    top_right_second_pt = apply_filter(filts(:,2), win.*u, use_zp)*sig_cw*sig_fw;
elseif which_part==2
    top_right_second_pt = -imag(hilbert(apply_filter(filts(:,2), win.*u, use_zp)))*sig_cw*sig_fw;
end

top_left = sum(u.*win.*f)*(sig_fw/sig_cw *apply_filter(filts(:,1), win.*real(c_orig), use_zp) + sig_cw/sig_fw *apply_filter(filts(:,2), win.*real(analytic_subbands(:,2)), use_zp));

grad = (top_right_first_pt + top_right_second_pt - top_left)/(sig_cw^2 * sig_fw^2);
