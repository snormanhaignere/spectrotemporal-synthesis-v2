function H = filt_temp_mod(fc_Hz, N, sr_Hz, LOWPASS, HIGHPASS)

% Temporal modulation filters
% derived from gen_cort.m in the nsl toolbox
% 
% -- Example --
 

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

% irf
t = (0:N-1)'/sr_Hz * fc_Hz; % time * fc_Hz
if LOWPASS
    h = t.^2 .* exp(-3.5*t); % gammatone
else
    h = sin(2*pi*t) .* t.^2 .* exp(-3.5*t); % gammatone 
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


