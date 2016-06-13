function H = filt_spec_mod(fc_cycPoct, N, sr_oct, LOWPASS, HIGHPASS)

% H = filt_spec_mod(fc_cycPoct, N, sr_oct, LOWPASS, HIGHPASS)
% 
% Spectral modulation filters
% derived from gen_corf.m in the nsl toolbox
% 
% -- Example --
% 
% fc_cycPoct = 4;
% sr_oct = 50;
% n_spectral_smps = 4 * sr_oct / fc_cycPoct; % number of samples
% lowpass = 0;
% highpass = 0;
% 
% % transfer function of the filter
% H = filt_spec_mod(fc_cycPoct, n_spectral_smps, sr_oct, lowpass, highpass);
% 
% % plot TF
% figure;
% subplot(2,1,1);
% mod_f = fft_freqs_from_siglen(n_spectral_smps, sr_oct);
% plot(fftshift_nyqlast(mod_f), fftshift_nyqlast(abs(H)));
% xlabel('Mod Freq (Hz)'); ylabel('Magnitude');
% subplot(2,1,2);
% plot(fftshift_nyqlast(mod_f), unwrap(fftshift_nyqlast(phase(H))));
% xlabel('Mod Freq (Hz)'); ylabel('Magnitude');
% 
% % plot IRF
% figure;
% f = (0:n_spectral_smps-1)/sr_oct;
% plot(f, ifft(H));
% xlabel('Freq (Hz)'); ylabel('Resp');

% delta transfer function
% constant impulse response
if fc_cycPoct == 0
    H = zeros(N,1);
    H(1) = 1;
    return;
end

% lowpass / highpass flags
if nargin < 4, 
    LOWPASS = 0;
end

if nargin < 5
    HIGHPASS = 0;
end

% index of the nyquist if present (i.e. if N is even)
% otherwise maximum positive frequency
nyq_index = ceil((N+1)/2);

% positive frequencies
pos_freqs = sr_oct*(0:nyq_index-1)'/N;

% check center frequency is below the maximum positive frequency
if fc_cycPoct > nyq_index
    error('Error in gen_corf_nopad.m: center frequency exceeds nyquist')
end

% Gabor function
% R1 = freqs/abs(fc);
% C1 = 1/2/.3/.3;
% H0 = exp(-C1*(R1-1).^2) + exp(-C1*(R1+1).^2);

% transfer function for positive frequencies
% fc is the center frequency
R1 = (pos_freqs/abs(fc_cycPoct)) .^ 2;
H_pos_freqs = R1 .* exp(1-R1); 
clear R1;

% passband
if LOWPASS
    [~, maxi] = max(H_pos_freqs);
    H_pos_freqs(1:maxi-1) = ones(maxi-1, 1);
end
if HIGHPASS
    [~, maxi] = max(H_pos_freqs);
    H_pos_freqs(maxi+1:length(H_pos_freqs)) = ...
        ones(length(H_pos_freqs)-maxi, 1);
end

% negative frequencies
if mod(N,2)==0 % nyquist present
    H_neg_freqs = conj(flip(H_pos_freqs(2:nyq_index-1)));
else % nyquist absent
    H_neg_freqs = conj(flip(H_pos_freqs(2:nyq_index)));
end

% combine positive and negative frequencies
H = [H_pos_freqs; H_neg_freqs];
