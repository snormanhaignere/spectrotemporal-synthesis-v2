F = 5;
T = 5;
P.env_sr = 10;
P.logf_spacing = P.env_sr;
Hts = filt_spectemp_mod(0,2,F,T,P)

% ifft2()


% plot_2DFT(abs(Hts),[1 1]*P.env_sr);