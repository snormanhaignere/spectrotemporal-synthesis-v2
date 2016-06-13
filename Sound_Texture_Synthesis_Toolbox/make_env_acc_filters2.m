% [FILTS,HZ_CUTOFFS] = MAKE_ENV_ACC_FILTERS2(SIGNAL_LENGTH, SR, LAGS_SMP)
%
% makes filters for use in calculating a version of the envelope
% autocorrelation, in which the envelope is smoothed in proportion to the
% time lag (which is a sensible thing to do is the lags at which you are
% measuring are logarithmically spaced)
%

% Josh McDermott <jhm@mit.edu>

function [filts,Hz_cutoffs] = make_env_acc_filters2(signal_length, sr, lags_smp)

if rem(signal_length,2)==0 %even length
    nfreqs = signal_length/2;%does not include DC
    max_freq = sr/2;
    freqs = [0:max_freq/nfreqs:max_freq]; %go all the way to nyquist
else %odd length
    nfreqs = (signal_length-1)/2;
    max_freq = sr*(signal_length-1)/2/signal_length; %max freq is just under nyquist
    freqs = [0:max_freq/nfreqs:max_freq];
end   

N=length(lags_smp);

cos_filts = zeros(nfreqs+1,N);

nyquist = sr/2;

for k=1:N
    if k==1
        high_cutoff = 1/(4*lags_smp(k)/1000);
        low_cutoff = .5*1/(4*lags_smp(k)/1000);
    else
        high_cutoff = 1/(4*(lags_smp(k)-lags_smp(k-1))/1000);
        low_cutoff = .5*1/(4*(lags_smp(k)-lags_smp(k-1))/1000);
    end
    
    if high_cutoff>nyquist
        cos_filts(:,k) = 1;
    else
        l_ind = min(find(freqs>low_cutoff));
        h_ind = max(find(freqs<high_cutoff));
        if l_ind<h_ind
            cos_filts(1:l_ind-1,k)=1;
            cos_filts(l_ind:h_ind,k) = cos( (freqs(l_ind:h_ind)- freqs(l_ind))/(freqs(l_ind)-freqs(h_ind))*pi/2); %map cutoffs to [0: pi/2] interval
        else
            cos_filts(1:l_ind-1,k)=1;
        end
    end
    Hz_cutoffs(k) = high_cutoff;
end

filts = cos_filts;
