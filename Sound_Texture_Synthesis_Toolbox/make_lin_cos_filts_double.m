% [FILTS,HZ_CUTOFFS,FREQS] = MAKE_LIN_COS_FILTS_DOUBLE(SIGNAL_LENGTH, SR,
% N, LOW_LIM, HI_LIM)
%
% Returns 2*N+5 filters as column vectors of FILTS 
% filters have cosine-shaped frequency responses, with center frequencies
% equally spaced on a linear scale from LOW_LIM to HI_LIM
%
% This function returns a filterbank that is 2x overcomplete compared to
% MAKE_LIN_COS_FILTS (to get filterbanks that can be compared with each
% other, use the same value of N in both cases). Adjacent filters overlap
% by 75%.
%
% As in MAKE_LIN_COS_FILTS, FILTS also contains lowpass
% and highpass filters to cover the ends of the spectrum.
%
% HZ_CUTOFFS is a vector of the cutoff frequencies of each filter. 
%
% FREQS is a vector of frequencies the same length as FILTS, that can be
% used to plot the frequency response of the filters.
%
% The squared frequency responses of the filters sums to 1, so that they
% can be applied once to generate subbands and then again to collapse the
% subbands to generate a sound signal, without changing the frequency
% content of the signal.
%
% Filters are to be applied multiplicatively in the frequency domain and
% thus have a length that scales with the signal length (SIGNAL_LENGTH).
%
% SR is the sampling rate
%
% intended for use with GENERATE_SUBBANDS and COLLAPSE_SUBBANDS

function [filts,Hz_cutoffs,freqs] = make_lin_cos_filts_double(signal_length, sr, N, low_lim, hi_lim)

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

num_filters = 2*N+1;
%make cutoffs evenly spaced on a linear frequency scale
spacing = (hi_lim-low_lim)/(num_filters+1);%in Hz
center_freqs = linspace(low_lim+spacing, hi_lim-spacing, num_filters); %in Hz

for k=1:num_filters
    l = center_freqs(k)-2*spacing;
    h = center_freqs(k)+2*spacing;
    l_ind = find(freqs>l, 1 );
    h_ind = find(freqs<h, 1, 'last' );
    avg = (l+h)/2;
    rnge = (h-l);
    cos_filts(l_ind:h_ind,k) = cos((freqs(l_ind:h_ind) - avg)/rnge*pi); %map cutoffs to -pi/2, pi/2 interval
end

%add lowpass and highpass to get perfect reconstruction
filts = zeros(nfreqs+1,num_filters+4);
filts(:,3:num_filters+2) = cos_filts;
%lowpass filters go up to peaks of first, second cos filters
h_ind = find(freqs<center_freqs(1), 1, 'last' );
filts(1:h_ind,1) = sqrt(1 - filts(1:h_ind,3).^2);
h_ind = find(freqs<center_freqs(2), 1, 'last' );
filts(1:h_ind,2) = sqrt(1 - filts(1:h_ind,4).^2);
%highpass filters go down to peaks of last two cos filters
l_ind = find(freqs>center_freqs(num_filters), 1 );
filts(l_ind:nfreqs+1,num_filters+4) = sqrt(1 - filts(l_ind:nfreqs+1,num_filters+2).^2);
l_ind = find(freqs>center_freqs(num_filters-1), 1 );
filts(l_ind:nfreqs+1,num_filters+3) = sqrt(1 - filts(l_ind:nfreqs+1,num_filters+1).^2);

filts = filts/sqrt(2); %so that squared freq response adds to 1

center_freqs = [center_freqs(1)-2*spacing center_freqs(2)-2*spacing center_freqs center_freqs(num_filters-1)+2*spacing center_freqs(num_filters)+2*spacing];
Hz_cutoffs = center_freqs;
Hz_cutoffs(Hz_cutoffs<0) = 1;

%subplot(2,1,1); plot(freqs,sum(filts.^2,2))
%subplot(2,1,2); semilogx(freqs,sum(filts.^2,2))
