function f = fftshift_nyqlast(f)

% same as fftshift but puts the nyquist last in the array
% 
% -- Examples --
% 
% f = fft_freqs_from_siglen(6);
% fftshift(f)
% fftshift_nyqlast(f)
% 
% N = 10;
% x = randn(N,1);
% ftx = fft(x);
% f = fft_freqs_from_siglen(N);
% plot(fftshift_nyqlast(f), fftshift_nyqlast(abs(ftx).^2));


f = circshift(f, -ceil((size(f)+1)/2));
