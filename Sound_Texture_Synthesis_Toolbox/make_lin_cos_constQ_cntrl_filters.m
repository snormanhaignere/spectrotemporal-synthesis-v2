% [FILTS, CFS, FREQS] = MAKE_LIN_COS_CONSTQ_CNTRL_FILTERS(SIGNAL_LENGTH, SR, N, LOW_LIM, HI_LIM, Q)
%
% returns N filters as column vectors of FILTS 
% filters have cosine-shaped frequency responses on a linear amplitude
% scale, with center frequencies equally spaced on a linear scale from
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
%


function [filts, Cfs, freqs] = make_lin_cos_constQ_cntrl_filters(signal_length, sr, N, low_lim, hi_lim, Q)

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

%compure amount of overlap for Q value
Q_Cfs = 2.^([log2(low_lim) : (log2(hi_lim)-log2(low_lim))/(N-1) : log2(hi_lim)]); %center_freqs for const Q
low_cut = Q_Cfs-Q_Cfs/Q/2;
hi_cut = Q_Cfs+Q_Cfs/Q/2;
low_bw = hi_cut(5)-low_cut(5);
overlap = hi_cut(5)-low_cut(5+1);
hi_bw = hi_cut(5+1)-low_cut(5+1);
prop_overlap_hi = overlap/low_bw;
prop_overlap_low = overlap/hi_bw;
avg_overlap = (prop_overlap_hi + prop_overlap_low)/2;

%make center frequencies evenly spaced on a linear scale
%want highest cos filter to go up to hi_lim
Cfs = [low_lim : (hi_lim-low_lim)/(N-1) : hi_lim];
spacing = (hi_lim-low_lim)/(N-1);
bw = spacing/(1-avg_overlap);

%easy-to-implement version: filters are symmetric on linear scale
for k=1:N
    l = Cfs(k)-bw;
    h = Cfs(k)+bw;
    if l<=0 %goes down to zero to avoid including DC
        l = 0;
        l_ind = min(find(freqs>l));
        h_ind = max(find(freqs<h));
        center = Cfs(k);
        c_ind = max(find(freqs<center));
        rnge1 = (center-l)*2;
        rnge2 = (h-center)*2;
        cos_filts(l_ind:c_ind,k) = cos((freqs(l_ind:c_ind) - center)/rnge1*pi); %map cutoffs to -pi/2, pi/2 interval
        cos_filts(c_ind:h_ind,k) = cos((freqs(c_ind:h_ind) - center)/rnge2*pi); %map cutoffs to -pi/2, pi/2 interval
    else
        l_ind = min(find(freqs>l));
        h_ind = max(find(freqs<h));
        avg = Cfs(k);
        rnge = h-l;
        cos_filts(l_ind:h_ind,k) = cos((freqs(l_ind:h_ind) - avg)/rnge*pi); %map cutoffs to -pi/2, pi/2 interval
    end
end

temp = sum(cos_filts'.^2);
filts=cos_filts/sqrt(mean(temp(find(freqs>=Cfs(4) & freqs<=Cfs(end-3)))));

