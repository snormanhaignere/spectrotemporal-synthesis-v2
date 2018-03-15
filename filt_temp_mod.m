function H = filt_temp_mod(fc_Hz, N, sr_Hz, LOWPASS, HIGHPASS, CAUSAL, BW, ...
    WAVELET, RANDOM_PHASE, RANDOM_FILT, RANDOM_SEED)

% Temporal modulation filters
% derived from gen_cort.m in the nsl toolbox
%
% -- Example --
%
% fc_Hz = 2;
% sr_Hz = 100;
% N = 6 * sr_Hz / fc_Hz; % number of samples
% lowpass = 0;
% highpass = 0;
% causal = 0;
%
% % transfer function of the filter
% H = filt_temp_mod(fc_Hz, N, sr_Hz, lowpass, highpass, causal);
% h = real(ifft(H));
%
% % plot IRF
% figure;
% t = (0:N-1)/sr_Hz;
% plot(t, h);
% xlabel('Freq (Hz)'); ylabel('Resp');
%
% % plot TF
% figure;
% subplot(2,1,1);
% mod_f = fft_freqs_from_siglen(N, sr_Hz);
% plot(fftshift_nyqlast(mod_f), fftshift_nyqlast(abs(H)));
% xlabel('Mod Freq (Hz)'); ylabel('Magnitude');
% subplot(2,1,2);
% plot(fftshift_nyqlast(mod_f), unwrap(fftshift_nyqlast(phase(H))));
% xlabel('Mod Freq (Hz)'); ylabel('Magnitude');
%
% 2017-06-28: Added non-causal filters
%
% 2018-03-13: Made it possible to change the bandwidths of the filter by
% altering the duration/extent of the envelope

% delta transfer function
% constant impulse response
if fc_Hz == 0
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
    CAUSAL = true;
end

if nargin < 7
    BW = 1;
end

if nargin < 8
    WAVELET = 1;
end

if nargin < 9
    RANDOM_PHASE = false;
end

if nargin < 10
    RANDOM_FILT = false;
end

if nargin < 11
    RANDOM_SEED = 1;
end

ResetRandStream2(RANDOM_SEED + round(fc_Hz*1000));

if RANDOM_PHASE || RANDOM_FILT
    N_orig = N;
    N = ceil(sr_Hz * (2/fc_Hz)/BW);
end

if RANDOM_FILT
    if N >= N_orig
        h = randn(N_orig,1);
    else
        if ~CAUSAL
            N_pad = N_orig-N;
            h = randn(N,1);
            h = [h(1:ceil(N/2)); zeros(ceil(N_pad),1); h(ceil(N/2)+1:end)]; % not exact
        else
            h = [randn(N,1); zeros(N_orig-N,1)];
        end
    end
    H = fft(h);
    return
end

% irf
switch WAVELET
    
    case 'gammatone'
        
        t = (0:N-1)'/sr_Hz;
        
        if ~CAUSAL
            t = t + 2/(3.5 * fc_Hz); % time * fc_Hz
            ti = t >= N/sr_Hz;
            t(ti) = mod(t(ti), N/sr_Hz);
            clear ti;
        end
        if LOWPASS
            h = (fc_Hz*t*BW).^2 .* exp(-3.5*fc_Hz*t*BW); % gammatone
        else
            h = sin(2*pi*t*fc_Hz) .* (t*fc_Hz*BW).^2 .* exp(-3.5*fc_Hz*t*BW); % gammatone
        end

    case 'morlet'
        
        
        % signal domain
        center_bin = (floor(N/2)+1);
        t = ((1:N) - center_bin)'/sr_Hz;
        
        if LOWPASS
            h = exp(-(fc_Hz*t*BW).^2);
        else
            h = cos(2*pi*t*fc_Hz) .* exp(-(fc_Hz*t*BW).^2);
        end
        
        % shift to create causal filter
        if ~CAUSAL
            h = circshift(h, [0, -(center_bin-1)]);
            t = circshift(t, [0, -(center_bin-1)]);
        end
            
    otherwise
        error('input:optargs', 'WAVELET must be gammatone or morlet, not %s', WAVELET);
end
        
% magnitude and phase of gammatone
H0 = fft(h);
th = angle(H0);
A = abs(H0);

% normalize magnitude of gammatone to have max of 1
A = A / max(A);

% create lowpass or highpass version if desired
if HIGHPASS
    
    % nyquist if present, otherwise highest pos. freq
    nyq_index = ceil((N+1)/2);
    
    % index of pos/neg frequencies with the maximum amplitude
    [~, maxA_index_pos] = max(A(1:nyq_index));
    [~, maxA_index_neg] = max(A(nyq_index+1:end));
    maxA_index_neg = maxA_index_neg + nyq_index;
    
    % check max is 1
    assert(all(A([maxA_index_pos, maxA_index_neg]) == 1))
    
    if HIGHPASS
        A(maxA_index_pos:maxA_index_neg) = 1;
    end
    
end

% recombine magnitude and phase
if RANDOM_PHASE
    max_pos_freq = ceil((N+1)/2);
    
    % negative frequencies
    if mod(N,2)==0 % nyquist present
        H_pos_freqs = A(1:max_pos_freq);
        H_pos_freqs(2:max_pos_freq-1) = H_pos_freqs(2:max_pos_freq-1).*exp(sqrt(-1)*unifrnd(0,pi,[max_pos_freq-2,1]));
        H_neg_freqs = conj(flip(H_pos_freqs(2:max_pos_freq-1)));
    else % nyquist absent
        H_pos_freqs = A(1:max_pos_freq);
        H_pos_freqs(2:max_pos_freq) = H_pos_freqs(2:max_pos_freq).*exp(sqrt(-1)*unifrnd(0,pi,[max_pos_freq-1,1]));
        H_neg_freqs = conj(flip(H_pos_freqs(2:max_pos_freq)));
    end
    H = [H_pos_freqs; H_neg_freqs];
    h = real(ifft(H));
    if N >= N_orig
        h = h(1:N_orig);
    else
        if ~CAUSAL
            N_pad = N_orig-N;
            h = [h(1:ceil(N/2)); zeros(ceil(N_pad),1); h(ceil(N/2)+1:end)]; % not exact
        else
            h = [h; zeros(N_orig-N,1)];
        end
    end
    H = fft(h);
else
    H = A .* exp(sqrt(-1)*th);
end
