% FUNCTION CX = AUTOCORR_MULT(X)
%
%
% works on matrices of column vectors X
%
% Dec 2012 -- Josh McDermott

function Cx = autocorr_mult(X)

Xf=fft(X);
Xf2=abs(Xf).^2;
Cx2=real(ifft(Xf2));
for j=1:size(Cx2,2)
    Cx(:,j) = fftshift(Cx2(:,j));
end