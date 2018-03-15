function H = filt_spec_mod(fc_cycPoct, N, sr_oct, LOWPASS, HIGHPASS, BW, ...
    WAVELET, RANDOM_PHASE, RANDOM_FILT, RANDOM_SEED)

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
% 
% 2018-03-13: Updated to allow for a morlet wavelet instead of a mexican hat
% wavelet. Also made it possible to alter the bandwidth for the morlet wavelet

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

if nargin < 6
    BW = 1;
end

if nargin < 7
    WAVELET = 'mexicanhat';
end

if nargin < 8
    RANDOM_PHASE = false;
end

if nargin < 9
    RANDOM_FILT = false;
end

if nargin < 10
    RANDOM_SEED = 1;
end

ResetRandStream2(RANDOM_SEED + round(fc_cycPoct*1000) + 1);

if strcmp(WAVELET, 'mexicanhat')
    if BW ~= 1
        error('spectemp:input',...
        'Can''t change bandwidth of the mexican hat\nBW must be 1');
    end
end

if RANDOM_PHASE || RANDOM_FILT
    N_orig = N;
    N = ceil(sr_oct * (2/fc_cycPoct)/BW);
end

if RANDOM_FILT
    if N >= N_orig
        h = randn(N_orig,1);
    else
        N_pad = N_orig-N;
        h = randn(N,1);
        h = [h(1:ceil(N/2)); zeros(ceil(N_pad),1); h(ceil(N/2)+1:end)]; % not exact
    end
    H = fft(h);
    return
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

% Note that in the signal domain the mexican hat filter equals:
% (1 - 2*(fc_cycPoct*pi*pos_freqs).^2) .* exp(-(fc_cycPoct*pi*pos_freqs).^2)

switch WAVELET

    case 'mexicanhat'
        
        % transfer function for positive frequencies
        % fc is the center frequency
        if LOWPASS
            per = 1/abs(fc_cycPoct);
            sig = per/2; % period = 2*sigma
            a = 1/(sig.^2);
            H_pos_freqs = exp(-(pi^2) * pos_freqs.^2 / a);
            % see: http://mathworld.wolfram.com/FourierTransformGaussian.html
        else
            f2 = (pos_freqs).^2 / abs(fc_cycPoct).^2;
            H_pos_freqs = f2 .* exp(-f2);
        end
        clear R1;
        
    case 'morlet'
        
        % signal domain
        center_bin = (floor(N/2)+1);
        f = ((1:N) - center_bin)/sr_oct;
        if LOWPASS
            h = exp(-(fc_cycPoct*f*BW).^2);
        else
            h = cos(2*pi*f*fc_cycPoct) .* exp(-(fc_cycPoct*f*BW).^2);
        end
        
        % shift to create causal filter
        h = circshift(h, [0, -(center_bin-1)]);
        
        % frequency domain magnitude and phase of gammatone
        H0 = real(fft(h));
        H_pos_freqs = H0(1:nyq_index)';        
        
    otherwise
        
        error('No matching wavelet for %s', WAVELET);
        
end

% normalize maximum frequency
H_pos_freqs = H_pos_freqs / max(H_pos_freqs);

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

if RANDOM_PHASE
    H_pos_freqs = abs(H_pos_freqs) .* exp(sqrt(-1) * unifrnd(0, 2*pi, size(H_pos_freqs)));
end

% negative frequencies
if mod(N,2)==0 % nyquist present
    H_neg_freqs = conj(flip(H_pos_freqs(2:nyq_index-1)));
else % nyquist absent
    H_neg_freqs = conj(flip(H_pos_freqs(2:nyq_index)));
end

% combine positive and negative frequencies
H = [H_pos_freqs; H_neg_freqs];


if RANDOM_PHASE
    h = real(ifft(H));
    if N >= N_orig
        h = h(1:N_orig);
    else
        N_pad = N_orig-N;
        h = [h(1:ceil(N/2)); zeros(ceil(N_pad),1); h(ceil(N/2)+1:end)]; % not exact
    end
    H = fft(h);
end



