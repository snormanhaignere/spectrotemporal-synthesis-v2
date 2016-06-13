function [f, nyq_index] = fft_freqs_from_siglen(N, sr)

% Returns a vector of frequencies corresponding to the FT of a signal of length
% N with sampling rate sr.
% 
% -- Examples --
% 
% fft_freqs_from_siglen(5, 1)
% fft_freqs_from_siglen(6, 1)
% fft_freqs_from_siglen(6, 100)
%
% N = 10;
% x = randn(N,1);
% ftx = fft(x);
% f = fft_freqs_from_siglen(N,1);
% plot(fftshift_nyqlast(f), fftshift_nyqlast(abs(ftx).^2));

if nargin < 2
    sr = 1;
end

% nyquist information if present
nyq_present = mod(N,2)==0;
if nyq_present
    nyq = sr/2;
    nyq_index = N/2+1;
else
    nyq = [];
    nyq_index = [];
end

% strictly positive frequencies
% excluding the DC and nyquist
n_positive_freqs = (N - 1 - nyq_present)/2;
pos_freqs = sr * (1:n_positive_freqs)/N;

% compose them
f = transpose([0, pos_freqs, nyq, -flip(pos_freqs)]);
