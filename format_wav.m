function wav = format_wav(wav, wav_sr, P)

% Format a waveform according to the parameters in struct P
% 
% 2017-07-10: Removed function that removes DC, and changed variance
% normalization to 0.01*rms(wav)

% convert to mono if stereo
wav = mean(wav,2);

% shorten
duration_sec = length(wav) / wav_sr;
if duration_sec > P.max_duration_sec
    wav = wav( 1 : wav_sr * P.max_duration_sec );
end
clear duration_sec;

% resample to desired audio rate
wav = resample(wav, P.audio_sr, wav_sr);
clear wav_sr;

% set RMS to 0.01
wav = 0.01 * wav / rms(wav(:));
