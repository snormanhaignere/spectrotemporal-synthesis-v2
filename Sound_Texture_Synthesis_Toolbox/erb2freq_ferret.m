function freq_Hz = erb2freq_ferret(n_erb)

% Inverse of freq2erb_ferret.m
% 
% 2018-03-20: Created by Sam NH

freq_kHz = (n_erb * 0.31 * (1-0.533)).^(1/(1-0.533)); 
freq_Hz = freq_kHz*1000;