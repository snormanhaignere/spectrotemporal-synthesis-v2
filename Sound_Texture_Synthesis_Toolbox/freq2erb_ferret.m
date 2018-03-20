function n_erb = freq2erb_ferret(freq_Hz)

% ERB bandwidths based on Sumner & Palmer:
% 
% ERB = 0.31 * freq_kHz^0.533
% 
% N_ERB derived by integrating the reciprocal of the equation above.
% 
% 2018-03-20: Created by Sam NH

freq_kHz = freq_Hz/1000;
n_erb = (1/0.31) * (1/(1-0.533)) * freq_kHz.^(1-0.533);