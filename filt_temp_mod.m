function H = filt_temp_mod(fc_Hz, N, sr_Hz, LOWPASS, HIGHPASS, CAUSAL, BW, randomphase)

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
    randomphase = false;
end

% irf
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
H = A .* exp(sqrt(-1)*th);


