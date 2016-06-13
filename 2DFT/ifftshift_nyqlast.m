function f = ifftshift_nyqlast(f)

% inverts fftshift_nyqlast
% 
% -- Example --
% 
% f = fft_freqs_from_siglen(6)
% fftshift_nyqlast(f)
% ifftshift_nyqlast(fftshift_nyqlast(f))

f = circshift(f, ceil((size(f)+1)/2));
