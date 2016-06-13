% [AC] = AUTOCORR_MULT_ZP(S, WIN_CHOICE, UNDO_WIN)
%
%
% applies a window to signal and zero-pads before taking autocorr
%
% works on arrays of column vectors
%
% Dec 2012 -- Josh McDermott

function [AC] = autocorr_mult_zp(S, win_choice, undo_win)

N=size(S,2)-2;
%begin by windowing and zero-padding
L = length(S);
wt=[1:L]/L;
if win_choice==1 %hanning
    w = 0.5-0.5*cos(2*pi*wt);
elseif win_choice==2 %rect
    w = ones(size(wt));
elseif win_choice==3 %hamming
    w = 0.54-0.46*cos(2*pi*wt);
elseif win_choice==4 %hamming
    w = 0.6-0.4*cos(2*pi*wt);
elseif win_choice==5 %welch
    w = sin(pi*wt);
end

S_w = S.*(w'*ones(1,N+2));
S_wp = [zeros(L/2,N+2); S_w; zeros(L/2, N+2)];

w_p = [zeros(1,length(w)/2) w zeros(1,length(w)/2)]';

AC = autocorr_mult(S_wp);

if undo_win
    w_ac = autocorr_mult(w_p);
    AC = AC./(w_ac*ones(1,N+2));
end

AC = AC(L/2+1:3*L/2,:);