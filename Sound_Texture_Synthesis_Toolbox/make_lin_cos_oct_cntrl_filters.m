% [FILTS,CFS,FREQS] = MAKE_LIN_COS_OCT_CNTRL_FILTERS(SIGNAL_LENGTH, SR, N,
% LOW_LIM, HI_LIM)
%
% Returns N filters as column vectors of FILTS 
% filters have cosine-shaped frequency responses on a linear frequency
% scale, with center frequencies equally spaced on a linear frequency scale
% from LOW_LIM to HI_LIM
%
% Adjacent filters overlap by 50%.
%
% CFS is a vector of the center frequencies of each filter. Because
% of the overlap arrangement, the upper cutoff of one filter is the center
% frequency of its neighbor.
%
% FREQS is a vector of frequencies the same length as FILTS, that can be
% used to plot the frequency response of the filters.
%
% FILTS does NOT contain lowpass and highpass % filters to cover the ends
% of the spectrum. The lowest filter cuts off just above 0 Hz, and thus has
% an approximately half-cosine frequency response
%
% The squared frequency responses of the filters sums to 1 between LOW_LIM
% and HI_LIM.
%
% Filters are to be applied multiplicatively in the frequency domain and
% thus have a length that scales with the signal length (SIGNAL_LENGTH).
%
% SR is the sampling rate
%
% intended to produce a linearly-spaced filter bank that can be swapped out
% for that produced by MAKE_OCTAVE_COS_FILTERS2

% Dec 2012 -- Josh McDermott <jhm@mit.edu>

function [filts,Cfs,freqs] = make_lin_cos_oct_cntrl_filters(signal_length, sr, N, low_lim, hi_lim)

if rem(signal_length,2)==0 %even length
    nfreqs = signal_length/2;%does not include DC
    max_freq = sr/2;
    freqs = [0:max_freq/nfreqs:max_freq]; %go all the way to nyquist
else %odd length
    nfreqs = (signal_length-1)/2;
    max_freq = sr*(signal_length-1)/2/signal_length; %max freq is just under nyquist
    freqs = [0:max_freq/nfreqs:max_freq];
end   
cos_filts = zeros(nfreqs+1,N);

if hi_lim>sr/2
    hi_lim = max_freq;
end

%make center frequencies evenly spaced on a linear scale
Cfs = [low_lim : (hi_lim-low_lim)/(N) : hi_lim];
cutoffs = [0 Cfs];
for k=1:N
    if k==1 %asymmetric in order to represent low frequencies
        l = cutoffs(k);
        h = cutoffs(k+2); %adjacent filters overlap by 50%
        l_ind = min(find(freqs>l));
        h_ind = max(find(freqs<h));
        center = cutoffs(k+1);
        c_ind = max(find(freqs<center));
        rnge1 = (center-l)*2;
        rnge2 = (h-center)*2;
        cos_filts(l_ind:c_ind,k) = cos((freqs(l_ind:c_ind) - center)/rnge1*pi); %map cutoffs to -pi/2, pi/2 interval
        cos_filts(c_ind:h_ind,k) = cos((freqs(c_ind:h_ind) - center)/rnge2*pi); %map cutoffs to -pi/2, pi/2 interval
    else
        l = cutoffs(k);
        h = cutoffs(k+2); %adjacent filters overlap by 50%
        l_ind = min(find(freqs>l));
        h_ind = max(find(freqs<h));
        avg = (l+h)/2;
        rnge = (h-l);
        cos_filts(l_ind:h_ind,k) = cos((freqs(l_ind:h_ind) - avg)/rnge*pi); %map cutoffs to -pi/2, pi/2 interval
    end
end

filts = cos_filts;
Cfs = Cfs(1:end-1);

%subplot(2,1,1); plot(freqs,sum(filts.^2,2))
%subplot(2,1,2); semilogx(freqs,sum(filts.^2,2))
