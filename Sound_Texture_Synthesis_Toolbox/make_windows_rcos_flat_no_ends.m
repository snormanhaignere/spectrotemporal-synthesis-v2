
% WINDOWS = MAKE_WINDOWS_RCOS_FLAT_NO_ENDS(SIGNAL_LENGTH_SMP, N, RAMP_PROP)
%
% returns a set of N windowing vectors for a signal of length
% SIGNAL_LENGTH_SMP as columns of WINDOWS 
%
% this version makes windows symmetric - first and last are same size
%
% RAMP_PROP specifies what proportion of window should be
% occupied by ramping (raised cosine)
%

%
% this version leaves out the end windows
%
% Dec 2012 -- Josh McDermott <jhm@mit.edu>

function windows = make_windows_rcos_flat_no_ends(signal_length_smp, N, ramp_prop)

%N = floor(signal_length_smp/desired_win_half_length_smp);

N=N+2;
if N==1
    windows = 2*ones(signal_length_smp,1);
else
    if ramp_prop==.5
        ramp_length_smp = floor(signal_length_smp/ (N-1));
        flat_length_smp = 0;
    elseif ramp_prop<.5
        flat_length = signal_length_smp/ (N*(1-ramp_prop)/(1-2*ramp_prop) - ramp_prop/(1-2*ramp_prop));
        ramp_length_smp = floor(flat_length*ramp_prop/(1-2*ramp_prop));
        flat_length_smp = floor(flat_length);
    else
        error('ramp_prop must be less than .5');
    end
    win_length_smp = flat_length_smp + 2*ramp_length_smp;
    windows = zeros(signal_length_smp,N);
    windows(1:flat_length_smp,1) = 2;
    windows(flat_length_smp+1:flat_length_smp+ramp_length_smp,1) = cos([1:ramp_length_smp]/ramp_length_smp*pi)+1;
    start_pt = flat_length_smp;
    for n = 1:N-2
        windows( start_pt+1:start_pt+ramp_length_smp,n+1) = cos([-ramp_length_smp+1:0]/ramp_length_smp*pi)+1;
        windows( start_pt+ramp_length_smp+1:start_pt+ramp_length_smp+flat_length_smp,n+1) = 2;
        windows( start_pt+ramp_length_smp+flat_length_smp+1:start_pt+2*ramp_length_smp+flat_length_smp,n+1) = cos([1:ramp_length_smp]/ramp_length_smp*pi)+1;
        start_pt = start_pt + flat_length_smp+ramp_length_smp;
    end
    windows(start_pt+1 : start_pt+ramp_length_smp,N) = cos([-ramp_length_smp+1:0]/ramp_length_smp*pi)+1;
    windows(start_pt+ramp_length_smp+1 : signal_length_smp,N) = 2;
    windows = windows(:,2:end-1);
end
windows=windows/2;
