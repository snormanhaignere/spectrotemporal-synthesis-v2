
% [FILTS,HZ_CENTER_FREQS,FREQS] = MAKE_OCTAVE_COS_FILTERS2(SIGNAL_LENGTH,
% SR, LOW_LIM, HI_LIM)
%
% Returns octave-spaced filters as column vectors of FILTS 
% filters have cosine-shaped frequency responses, with center frequencies
% equally spaced on a log frequency scale from LOW_LIM to HI_LIM
%
% Adjacent filters overlap by 50%.
%
% HZ_CENTER_FREQS is a vector of the center frequencies of each filter. Because
% of the overlap arrangement, the upper cutoff of one filter is the center
% frequency of its neighbor.
%
% FREQS is a vector of frequencies the same length as FILTS, that can be
% used to plot the frequency response of the filters.
%
% FILTS does NOT contain lowpass and highpass % filters to cover the ends
% of the spectrum.
%
% The squared frequency responses of the filters sums to 1 between LOW_LIM
% and HI_LIM.
%
% Filters are to be applied multiplicatively in the frequency domain and
% thus have a length that scales with the signal length (SIGNAL_LENGTH).
%
% SR is the sampling rate
%

% Dec 2012 -- Josh McDermott <jhm@mit.edu>


function [filts,Hz_center_freqs,freqs] = make_octave_cos_filters2(signal_length, sr, low_lim, hi_lim)

if rem(signal_length,2)==0 %even length
    nfreqs = signal_length/2;%does not include DC
    max_freq = sr/2;
    freqs = [0:max_freq/nfreqs:max_freq]; %go all the way to nyquist
else %odd length
    nfreqs = (signal_length-1)/2;
    max_freq = sr*(signal_length-1)/2/signal_length; %max freq is just under nyquist
    freqs = [0:max_freq/nfreqs:max_freq];
end   

if hi_lim>sr/2
    hi_lim = max_freq;
end

%make cutoffs octave spaced
cutoffs = hi_lim./(2.^[0:20]);
cutoffs = fliplr(cutoffs(find(cutoffs>low_lim)));
center_freqs = cutoffs(1:end-1);
N=length(center_freqs);

cos_filts = zeros(nfreqs+1,N);
for k=1:N
    l = center_freqs(k)/2;
    h = cutoffs(k)*2; %adjacent filters overlap by 50%
    l_ind = min(find(freqs>l));
    h_ind = max(find(freqs<h));
    avg = (log2(l)+log2(h))/2;
    rnge = (log2(h)-log2(l));
    cos_filts(l_ind:h_ind,k) = cos((log2( freqs(l_ind:h_ind) ) - avg)/rnge*pi); %map cutoffs to -pi/2, pi/2 interval
end

filts = cos_filts;

Hz_center_freqs = center_freqs;

%subplot(2,1,1); plot(freqs,sum(filts.^2,2))
%subplot(2,1,2); semilogx(freqs,sum(filts.^2,2))
