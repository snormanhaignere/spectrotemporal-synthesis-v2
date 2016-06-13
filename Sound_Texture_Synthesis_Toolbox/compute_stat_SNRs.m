function [SNRs] = compute_stat_SNRs(target_S, synth_S, P, SNRs, it_n)

% [SNRS] = COMPUTE_STAT_SNRS(TARGET_S, SYNTH_S, P, SNRS, IT_N)
%
% computes SNRs of statistics of synthetic signal
% incorporates SNRs into existing structure SNRs in the IT_Nth position of
% each variable
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

% compute snr including only subbands with significant power
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

env_C_error = sum(sum((target_S.env_C(env_C_to_include==1)-synth_S.env_C(env_C_to_include==1)).^2));
env_C_error_all = sum(sum((target_S.env_C(env_C_to_include_full==1)-synth_S.env_C(env_C_to_include_full==1)).^2));
env_C_power = sum(sum(target_S.env_C(env_C_to_include==1).^2));
env_C_power_all = sum(sum(target_S.env_C(env_C_to_include_full==1).^2));

for C1_band = 1:size(target_S.mod_C1,3)
    temp1 = target_S.mod_C1(:,:,C1_band); temp2 = synth_S.mod_C1(:,:,C1_band);
    mod_C1_error(C1_band) = sum(sum((temp1(mod_C1_to_include==1) - temp2(mod_C1_to_include==1)).^2));
    mod_C1_error_all(C1_band) = sum(sum((temp1(mod_C1_to_include_full==1) - temp2(mod_C1_to_include_full==1)).^2));
    mod_C1_power(C1_band) = sum(sum(temp1(mod_C1_to_include==1).^2));
    mod_C1_power_all(C1_band) = sum(sum(temp1(mod_C1_to_include_full==1).^2));
end


for k= 1:num_channels
    subband_hist_error(k) = sum((target_S.subband_hist(k,:)-synth_S.subband_hist(k,:)).^2);
    subband_env_hist_error(k) = sum((target_S.env_hist(k,:)-synth_S.env_hist(k,:)).^2);
    sub_var_error(k) = (target_S.subband_var(k)-synth_S.subband_var(k))^2;
    sub_kurt_error(k) = (target_S.subband_kurt(k)-synth_S.subband_kurt(k))^2;
    env_mean_error(k) = (target_S.env_mean(k)-synth_S.env_mean(k))^2;
    env_var_error(k) = (target_S.env_var(k)-synth_S.env_var(k))^2;
    env_skew_error(k) = (target_S.env_skew(k)-synth_S.env_skew(k))^2;
    env_kurt_error(k) = (target_S.env_kurt(k)-synth_S.env_kurt(k))^2;
    
    env_ac_error(k) = sum((target_S.env_ac(k,:)-synth_S.env_ac(k,:)).^2);
    
    mod_power_error(k) = sum((target_S.mod_power(k,:)-synth_S.mod_power(k,:)).^2);
    
    mod_C2_error(k) = sum(sum(squeeze(target_S.mod_C2(k,:,:)-synth_S.mod_C2(k,:,:)).^2));
    
    sub_ac_error(k) = sum((target_S.subband_ac{k}-synth_S.subband_ac{k}).^2);
end

SNRs.subband_hist(it_n) = 10*log10(sum(sum(target_S.subband_hist(sig_subs,:).^2))/sum(subband_hist_error(sig_subs)));
SNRs.subband_var(it_n) = 10*log10(sum(target_S.subband_var(sig_subs).^2)/sum(sub_var_error(sig_subs)));
SNRs.subband_kurt(it_n) = 10*log10(sum(target_S.subband_kurt(sig_subs).^2)/sum(sub_kurt_error(sig_subs)));
SNRs.subband_ac(it_n) = 10*log10(sum(target_S.subband_ac_power(sig_subs))/sum(sub_ac_error(sig_subs)));

SNRs.env_hist(it_n) = 10*log10(sum(sum(target_S.env_hist(sig_subs,:).^2))/sum(subband_env_hist_error(sig_subs)));
SNRs.env_mean(it_n) = 10*log10(sum(target_S.env_mean(sig_subs).^2)/sum(env_mean_error(sig_subs)));
SNRs.env_var(it_n) = 10*log10(sum(target_S.env_var(sig_subs).^2)/sum(env_var_error(sig_subs)));
SNRs.env_skew(it_n) = 10*log10(sum(target_S.env_skew(sig_subs).^2)/sum(env_skew_error(sig_subs)));
SNRs.env_kurt(it_n) = 10*log10(sum(target_S.env_kurt(sig_subs).^2)/sum(env_kurt_error(sig_subs)));

SNRs.env_ac(it_n) = 10*log10(sum(sum(target_S.env_ac(sig_subs,:).^2))/sum(env_ac_error(sig_subs)));
SNRs.mod_power(it_n) = 10*log10(sum(sum(target_S.mod_power(sig_subs,:).^2))/sum(mod_power_error(sig_subs)));

SNRs.env_C(it_n) = 10*log10(env_C_power/env_C_error);
SNRs.env_C_all(it_n) = 10*log10(env_C_power_all/env_C_error_all);

SNRs.mod_C1(it_n) = 10*log10(sum(mod_C1_power)/sum(mod_C1_error));
SNRs.mod_C1_all(it_n) = 10*log10(sum(mod_C1_power_all)/sum(mod_C1_error_all));

SNRs.mod_C2(it_n) = 10*log10(sum(sum(sum(target_S.mod_C2(sig_subs,:,:).^2)))/sum(mod_C2_error(sig_subs)));


fprintf('sub hist snr = %3.2f; sub kurt snr = %3.2f; sub ac snr = %3.2f; env hist snr = %3.2f\n',...
    SNRs.subband_hist(it_n), SNRs.subband_kurt(it_n), SNRs.subband_ac(it_n), SNRs.env_hist(it_n));
fprintf('env mean snr = %3.2f; env var snr = %3.2f; env skew snr = %3.2f; env kurt snr = %3.2f\n',...
    SNRs.env_mean(it_n), SNRs.env_var(it_n), SNRs.env_skew(it_n), SNRs.env_kurt(it_n));
fprintf('env C snr = %3.2f; env C all snr = %3.2f; mod power snr = %3.2f; env ac snr = %3.2f\n',...
    SNRs.env_C(it_n), SNRs.env_C_all(it_n), SNRs.mod_power(it_n), SNRs.env_ac(it_n));
fprintf('mod C1 snr = %3.2f; mod C1 all snr = %3.2f; mod C2 snr = %3.2f\n\n',...
    SNRs.mod_C1(it_n), SNRs.mod_C1_all(it_n), SNRs.mod_C2(it_n));

