function display_before_after_error(target_S, synth_S1, synth_S2, P)

% DISPLAY_BEFORE_AFTER_ERROR(TARGET_S, SYNTH_S1, SYNTH_S2, P)
%
% computes and displays error in statistics before and after imposition
%
% rms error is used for correlation coefficients as the units are
% meaningful
%
% SNR is used for other stats
%

% This code is part of an instantiation of a sound texture synthesis
% algorithm developed with Eero Simoncelli and described in this paper:
%
% McDermott, J.H., Simoncelli, E.P. (2011) Sound texture perception via
% statistics of the auditory periphery: Evidence from sound synthesis.
% Neuron, 71, 926-940. 
%
% Dec 2012 -- Josh McDermott <jhm@mit.edu>
% May 2014 -- modified to fix incompatibility with overcomplete filter
% banks -- Josh McDermott

% compute error including only subbands with significant power
sub_power = 10*log10(max(target_S.subband_var)./target_S.subband_var);
sig_subs = find(sub_power<P.sig_sub_cutoff_db);

num_channels = length(sub_power);

env_C_to_include = zeros(num_channels);
env_C_to_include_full = zeros(num_channels);
mod_C1_to_include = zeros(num_channels);
mod_C1_to_include_full = zeros(num_channels);
for k= 1:(num_channels)
    if ismember(k,sig_subs)
        tmp=k+P.C_offsets_to_impose;
        env_C_to_include(k,tmp( tmp>0 & tmp<=(num_channels) ) )=1;
        env_C_to_include_full(k,k+1:end)=1;
        tmp=k+P.mod_C1_offsets_to_impose;
        mod_C1_to_include(k,tmp( tmp>0 & tmp<=(num_channels) ) )=1;
        mod_C1_to_include_full(k,k+1:end)=1;
    end
end

env_C_error1 = sqrt(mean(mean((target_S.env_C(env_C_to_include==1)-synth_S1.env_C(env_C_to_include==1)).^2)));
env_C_error_all1 = sqrt(mean(mean((target_S.env_C(env_C_to_include_full==1)-synth_S1.env_C(env_C_to_include_full==1)).^2)));
env_C_error2 = sqrt(mean(mean((target_S.env_C(env_C_to_include==1)-synth_S2.env_C(env_C_to_include==1)).^2)));
env_C_error_all2 = sqrt(mean(mean((target_S.env_C(env_C_to_include_full==1)-synth_S2.env_C(env_C_to_include_full==1)).^2)));

for C1_band = 1:size(target_S.mod_C1,3)
    temp1 = target_S.mod_C1(:,:,C1_band); temp2 = synth_S1.mod_C1(:,:,C1_band); temp3 = synth_S2.mod_C1(:,:,C1_band);
    mod_C1_error1(C1_band) = sqrt(mean(mean((temp1(mod_C1_to_include==1) - temp2(mod_C1_to_include==1)).^2)));
    mod_C1_error_all1(C1_band) = sqrt(mean(mean((temp1(mod_C1_to_include_full==1) - temp2(mod_C1_to_include_full==1)).^2)));
    mod_C1_error2(C1_band) = sqrt(mean(mean((temp1(mod_C1_to_include==1) - temp3(mod_C1_to_include==1)).^2)));
    mod_C1_error_all2(C1_band) = sqrt(mean(mean((temp1(mod_C1_to_include_full==1) - temp3(mod_C1_to_include_full==1)).^2)));
end
mod_C1_error1 = mean(mod_C1_error1); %average across mod bands
mod_C1_error_all1 = mean(mod_C1_error_all1); %average across mod bands
mod_C1_error2 = mean(mod_C1_error2); %average across mod bands
mod_C1_error_all2 = mean(mod_C1_error_all2); %average across mod bands

%mean is used instead of sum for stats where we want RMS error rather than
%SNR
for k= 1:(num_channels)
    subband_env_hist_error1(k) = sum((target_S.env_hist(k,:)-synth_S1.env_hist(k,:)).^2);
    subband_env_hist_error2(k) = sum((target_S.env_hist(k,:)-synth_S2.env_hist(k,:)).^2);
    env_mean_error1(k) = (target_S.env_mean(k)-synth_S1.env_mean(k))^2;
    env_var_error1(k) = (target_S.env_var(k)-synth_S1.env_var(k))^2;
    env_skew_error1(k) = (target_S.env_skew(k)-synth_S1.env_skew(k))^2;
    env_kurt_error1(k) = (target_S.env_kurt(k)-synth_S1.env_kurt(k))^2;
    env_mean_error2(k) = (target_S.env_mean(k)-synth_S2.env_mean(k))^2;
    env_var_error2(k) = (target_S.env_var(k)-synth_S2.env_var(k))^2;
    env_skew_error2(k) = (target_S.env_skew(k)-synth_S2.env_skew(k))^2;
    env_kurt_error2(k) = (target_S.env_kurt(k)-synth_S2.env_kurt(k))^2;
    
    env_ac_error1(k) = sqrt(mean((target_S.env_ac(k,:)-synth_S1.env_ac(k,:)).^2));
    env_ac_error2(k) = sqrt(mean((target_S.env_ac(k,:)-synth_S2.env_ac(k,:)).^2));
    
    mod_power_error1(k) = sum((target_S.mod_power(k,:)-synth_S1.mod_power(k,:)).^2);
    mod_power_error2(k) = sum((target_S.mod_power(k,:)-synth_S2.mod_power(k,:)).^2);
    
    mod_C2_error1(k) = sqrt(mean(sum(squeeze(target_S.mod_C2(k,:,:)-synth_S1.mod_C2(k,:,:)).^2)));
    mod_C2_error2(k) = sqrt(mean(sum(squeeze(target_S.mod_C2(k,:,:)-synth_S2.mod_C2(k,:,:)).^2)));
end

env_ac_error1 = mean(env_ac_error1); %average across bands
env_ac_error2 = mean(env_ac_error2);

mod_C2_error1 = mean(mod_C2_error1);
mod_C2_error2 = mean(mod_C2_error2);

SNRs1.env_hist = 10*log10(sum(sum(target_S.env_hist(sig_subs,:).^2))/sum(subband_env_hist_error1(sig_subs)));
SNRs1.env_mean = 10*log10(sum(target_S.env_mean(sig_subs).^2)/sum(env_mean_error1(sig_subs)));
SNRs1.env_var = 10*log10(sum(target_S.env_var(sig_subs).^2)/sum(env_var_error1(sig_subs)));
SNRs1.env_skew = 10*log10(sum(target_S.env_skew(sig_subs).^2)/sum(env_skew_error1(sig_subs)));
SNRs1.env_kurt = 10*log10(sum(target_S.env_kurt(sig_subs).^2)/sum(env_kurt_error1(sig_subs)));

SNRs2.env_hist = 10*log10(sum(sum(target_S.env_hist(sig_subs,:).^2))/sum(subband_env_hist_error2(sig_subs)));
SNRs2.env_mean = 10*log10(sum(target_S.env_mean(sig_subs).^2)/sum(env_mean_error2(sig_subs)));
SNRs2.env_var = 10*log10(sum(target_S.env_var(sig_subs).^2)/sum(env_var_error2(sig_subs)));
SNRs2.env_skew = 10*log10(sum(target_S.env_skew(sig_subs).^2)/sum(env_skew_error2(sig_subs)));
SNRs2.env_kurt = 10*log10(sum(target_S.env_kurt(sig_subs).^2)/sum(env_kurt_error2(sig_subs)));

SNRs1.mod_power = 10*log10(sum(sum(target_S.mod_power(sig_subs,:).^2))/sum(mod_power_error1(sig_subs)));
SNRs2.mod_power = 10*log10(sum(sum(target_S.mod_power(sig_subs,:).^2))/sum(mod_power_error2(sig_subs)));


fprintf('Old Mean SNR = %2.2f; New Mean SNR = %2.2f\n', SNRs1.env_mean, SNRs2.env_mean);
fprintf('Old Var SNR = %2.2f; New Var SNR = %2.2f\n', SNRs1.env_var, SNRs2.env_var);
fprintf('Old Skew SNR = %2.2f; New Skew SNR = %2.2f\n', SNRs1.env_skew, SNRs2.env_skew);
fprintf('Old Kurt SNR = %2.2f; New Kurt SNR = %2.2f\n', SNRs1.env_kurt, SNRs2.env_kurt);
fprintf('Old EnvHist SNR = %2.2f; New EnvHist SNR = %2.2f\n', SNRs1.env_hist, SNRs2.env_hist);
fprintf('Old MP SNR = %2.2f; New MP SNR = %2.2f\n', SNRs1.mod_power, SNRs2.mod_power);
fprintf('Old Env C avg error = %2.2f; New Env C avg error = %2.2f\n', env_C_error1, env_C_error2);
fprintf('Old Env AC avg error = %2.2f; New Env AC avg error = %2.2f\n', env_ac_error1, env_ac_error2);
fprintf('Old Mod C1 avg error = %2.2f; New Mod C1 avg error = %2.2f\n', mod_C1_error1, mod_C1_error2);
fprintf('Old Mod C2 avg error = %2.2f; New Mod C2 avg error = %2.2f\n', mod_C2_error1, mod_C2_error2);
fprintf('\n');


