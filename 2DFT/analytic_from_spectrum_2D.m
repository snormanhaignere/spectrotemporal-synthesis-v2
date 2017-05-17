function [FT, freqs_to_zero, freqs_to_double] = ...
    analytic_from_spectrum_2D(FT, dim, sgn)

% plot_2DFT(analytic_from_spectrum_2D(fft2(randn(6,6)),1))
% plot_2DFT(analytic_from_spectrum_2D(fft2(randn(6,6)),2))

if nargin < 3
    sgn = 1;
end

[~, freqs_to_zero, freqs_to_double] = ...
    analytic_from_spectrum(1:size(FT,dim), sgn);

if dim == 1
    FT(freqs_to_zero,:) = 0;
    FT(freqs_to_double,:) = 2*FT(freqs_to_double,:);
else
    FT(:,freqs_to_zero) = 0;
    FT(:,freqs_to_double) = 2*FT(:,freqs_to_double);
end