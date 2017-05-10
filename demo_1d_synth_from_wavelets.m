% Illustrates how to do analysis and synthesis from 1D temporal filters
% 
% 2017-05-10

% sampling rate
sr_Hz = 100;

% rates of the wavelets in Hz
fc_Hz = [0 0.5 1 2 4 8 16 32];

% create a 4-second long signal
N = 4*sr_Hz;
sig = rand(N,1)-0.5;

% compute subbands
subb = nan(N, length(fc_Hz));
for i = 1:length(fc_Hz)
    H = filt_temp_mod(fc_Hz(i), N, sr_Hz, 0, 0);
    subb(:,i) = real(ifft(fft(sig) .* H));
end
imagesc(subb);

% collapse subbands to reconstruct
FT_sig_recon = zeros(N, 1);
H_accum = zeros(N, 1);
for i = 1:length(fc_Hz)
    H = filt_temp_mod(fc_Hz(i), N, sr_Hz, 0, 0);
    FT_sig_recon = FT_sig_recon + fft(subb(:,i)) .* conj(H);
    H_accum = H_accum + H .* conj(H);
end
recon = real(ifft(FT_sig_recon ./ H_accum));



