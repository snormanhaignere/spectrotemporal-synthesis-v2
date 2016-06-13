% [FILTS, CFS, FREQS] = MAKE_CONSTQ_COS_FILTERS(SIGNAL_LENGTH, SR, N, LOW_LIM, HI_LIM, Q)
%
% returns N filters as column vectors of FILTS 
% filters have cosine-shaped frequency responses on a linear amplitude
% scale, with center frequencies equally spaced on a log scale from
% LOW_LIM to HI_LIM, with a specified Q-value (Q)
%
% CFS is a vector of the center frequencies of each filter. 
%
% FREQS is a vector of frequencies the same length as FILTS, that can be
% used to plot the frequency response of the filters.
%
% filters are to be applied multiplicatively in the frequency domain and
% thus have a length that scales with the signal length (SIGNAL_LENGTH)
%
% SR is the sampling rate
%
% the cost of being able to separately set the spacing and Q value is that
% the filter responses do not sum to 1
%
% these filters are thus not well-suited to analysis-synthesis subband
% decompositions, and for sound texture synthesis they are used for
% measurement of modulation statistics
%

% Dec 2012 -- Josh McDermott <jhm@mit.edu>

function [filts, Cfs, freqs] = make_constQ_cos_filters(signal_length, sr, N, low_lim, hi_lim, Q)

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

%make center frequencies evenly spaced on a log scale
%want highest cos filter to go up to hi_lim
Cfs = 2.^([log2(low_lim) : (log2(hi_lim)-log2(low_lim))/(N-1) : log2(hi_lim)]);

%easy-to-implement version: filters are symmetric on linear scale
for k=1:N
    bw = Cfs(k)/Q;
    l = Cfs(k)-bw; %so that half power point is at Cf-bw/2
    h = Cfs(k)+bw;
    l_ind = find(freqs>l, 1 );
    h_ind = find(freqs<h, 1, 'last' );
    avg = Cfs(k); %(log2(l+1)+log2(h+1))/2;
    rnge = h-l;%(log2(h+1)-log2(l+1));
    cos_filts(l_ind:h_ind,k) = cos((freqs(l_ind:h_ind) - avg)/rnge*pi); %map cutoffs to -pi/2, pi/2 interval
end

temp = sum(cos_filts'.^2);
filts=cos_filts/sqrt(mean(temp(freqs>=Cfs(4) & freqs<=Cfs(end-3))));

