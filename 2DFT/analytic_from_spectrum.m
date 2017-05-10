function [ftx, freqs_to_zero, freqs_to_double] = analytic_from_spectrum(ftx, sgn)

% Creates an analytic signal from an input spectrum, zeroing negative
% frequencies and doubling positive frequencies. 
% 
% -- Inputs --
% 
% ftx: Fourier transform of a signal
% 
% sgn: if sgn == -1, doubles negative frequencies and zeros positive frequencies
% (default value: 1)
% 
% -- Example showing equivalence with time-domain function "hilbert.m" --
% 
% x = randn(5,1);
% hilbert(x)
% ifft(analytic_from_spectrum(fft(x)))

if nargin < 2
    sgn = 1;
end

% frequencies for a signal of length(x)
[freqs, nyq_index] = fft_freqs_from_siglen(length(ftx),1);

% frequencies to zero out, excluding DC and nyquist
freqs_to_zero = find(sign(freqs) == -sgn);
freqs_to_zero = setdiff(freqs_to_zero, [1, nyq_index]);
ftx(freqs_to_zero) = 0;

% frequencies to zero out, excluding DC and nyquist
freqs_to_double = find(sign(freqs) == sgn);
freqs_to_double = setdiff(freqs_to_double, [1, nyq_index]);
ftx(freqs_to_double) = 2*ftx(freqs_to_double);