%
% assumes there is a struct SNRs from the synthesis process
%

% Jan 3 2013 -- Josh McDermott <jhm@mit.edu>

figure('Position',[400 100 800 700]);
subplot(2,1,1)
title(['SNR Plots for ' fig_title_string]);
hold on

plot(SNRs.subband_hist,'b');
plot(SNRs.subband_kurt,'b--');
plot(SNRs.env_hist,'r');
plot(SNRs.env_mean,'m');
plot(SNRs.env_var,'g');
plot(SNRs.env_skew,'c');
plot(SNRs.env_kurt,'k');
legend('Subband Hist', 'Subband Kurt','Env Hist','Env Mean','Env Var','Env Skew','Env Kurt', 'Location', 'EastOutside');
ylabel('dB');
xlabel('Iteration');

subplot(2,1,2)
hold on
plot(SNRs.mod_power,'b');
plot(SNRs.env_C,'r');
plot(SNRs.mod_C1,'m');
plot(SNRs.mod_C2,'g');
plot(SNRs.env_ac,'c');
plot(SNRs.subband_ac,'k');
legend('Mod Power', 'Env C', 'Mod C1', 'Mod C2', 'Env AC', 'Sub AC', 'Location', 'EastOutside');
ylabel('dB');
xlabel('Iteration');
