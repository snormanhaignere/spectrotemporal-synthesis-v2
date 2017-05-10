% Simple code illustrating how to do analysis and synthesis using the filters
% from the spectrotemporal model.
% 
% 2017-05-10: Created, Sam NH

%% Setup

% Parameters
P = synthesis_parameters_toy;

% read in waveform
[wav,sr] = audioread('speech.wav');
wav = wav(1:sr);

%% Cochleogram

% cochleogram
[coch, P, R] = wav2coch_without_filts(wav, P);

% plot cochleogram
figure;
plot_cochleogram(coch, P.f, P.t)

%% Single subbands

% single subband
spec_mod_rate = 2;
temp_mod_rate = 4;

% one-liner
coch_subband = coch2filtcoch(coch, spec_mod_rate, temp_mod_rate, P);

% break down one liner
FT_coch = fft2(coch);
Hts = filt_spectemp_mod(spec_mod_rate, temp_mod_rate, size(coch,2), size(coch,1), P);
coch_subband = real(ifft2(FT_coch .* Hts));

% plot subband
figure;
plot_cochleogram(coch_subband, P.f, P.t)

%% Multiple subbands

% many subbands
FT_coch = fft2(coch);
spec_mod_rates = P.spec_mod_rates;
temp_mod_rates = unique([-P.temp_mod_rates, P.temp_mod_rates]);
coch_subbands = nan(length(P.t), length(P.f), ...
    length(P.temp_mod_rates), length(P.spec_mod_rates));
for i = 1:length(temp_mod_rates)
    for j = 1:length(spec_mod_rates)
        
        % filter transfer function
        Hts = filt_spectemp_mod(spec_mod_rates(j), temp_mod_rates(i), ...
            size(coch,2), size(coch,1), P);    
        coch_subbands(:,:,i,j) = real(ifft2(FT_coch .* Hts));
    end
end

%% Invert back to cochleogram

% accumulate FT of matched cochleograms
accum_FT_subbands = zeros(size(coch));
accum_FT_transfer_function = zeros(size(coch));
for i = 1:length(temp_mod_rates)
    for j = 1:length(spec_mod_rates)
        
        % filter transfer function
        Hts = filt_spectemp_mod(spec_mod_rates(j), temp_mod_rates(i), ...
            size(coch,2), size(coch,1), P);        
        
        % accumulate FT of subbands
        accum_FT_subbands = ...
            accum_FT_subbands + fft2(coch_subbands(:,:,i,j)) .* conj(Hts);
        
        % accumulate FT of transfer functions
        accum_FT_transfer_function = accum_FT_transfer_function + Hts .* conj(Hts);
        
    end
end

% divide by accumulated transfer functions
coch_recon = real(ifft2(accum_FT_subbands ./ accum_FT_transfer_function));

% plot subband
figure;
subplot(2,1,1);
plot_cochleogram(coch, P.f, P.t);
subplot(2,1,2);
plot_cochleogram(coch_recon, P.f, P.t);

%% Reconstruct waveform

wav_recon = coch2wav_without_filts(coch_recon, P, R);

figure;
subplot(2,1,1);
plot(wav);
subplot(2,1,2);
plot(wav_recon);

