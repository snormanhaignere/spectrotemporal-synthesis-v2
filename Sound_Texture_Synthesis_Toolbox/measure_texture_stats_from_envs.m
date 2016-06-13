function S = measure_texture_stats_from_envs(subband_envs, P, measurement_win, ...
    mod_filts, env_ac_filts, mod_C1_filts, mod_C2_filts, sub_env_hist_bins)

% FUNCTION S = MEASURE_TEXTURE_STATS_FROM_ENVS(SUBBAND_ENVS, P, MEASUREMENT_WIN, ...
%   MOD_FILTS, ENV_AC_FILTS, MOD_C1_FILTS, MOD_C2_FILTS, SUB_ENV_HIST_BINS)
%
% function to measure texture statistics from envelopes of subbands of a
% sound signal contained as columns of SUBBAND_ENVS
%
% assumes that subband envelopes are compressed and downsampled as
% specificed in the parameters in P.
%
% MEASUREMENT_WIN specifies the weights on each sample of the signal in the average
% other parameters are specified in P.
%
% Uses same structure format as MEASURE_TEXTURE_STATS, but does not compute
% subband statistics, only envelope statistics.
%
% Intended for use to calculate stats following an iteration of the
% imposition procedure.
%
%

% This code is part of an instantiation of a sound texture synthesis
% algorithm developed with Eero Simoncelli and described in this paper:
%
% McDermott, J.H., Simoncelli, E.P. (2011) Sound texture perception via
% statistics of the auditory periphery: Evidence from sound synthesis.
% Neuron, 71, 926-940. 
%
% Dec 2012 -- Josh McDermott <jhm@mit.edu>


%compute marginals - of subband envelopes and modulation bands
for j=1:size(subband_envs,2) %go through subbands
    S.env_mean(j) = stat_central_moment_win(subband_envs(:,j),1,measurement_win);
    S.env_var(j) = stat_central_moment_win(subband_envs(:,j),2,measurement_win,S.env_mean(j));
    S.env_skew(j) = stat_central_moment_win(subband_envs(:,j),3,measurement_win,S.env_mean(j));
    S.env_kurt(j) = stat_central_moment_win(subband_envs(:,j),4,measurement_win,S.env_mean(j));
    [temp,bins]=hist(subband_envs(:,j),sub_env_hist_bins(j,:));
    S.env_hist(j,1:P.n_hist_bins)=temp/sum(temp);
    S.env_bins(j,1:P.n_hist_bins)=bins;

    S.env_ac(j,1:length(P.env_ac_intervals_smp)) = stat_env_ac_scaled_win(subband_envs(:,j),env_ac_filts,P.env_ac_intervals_smp,P.use_zp,measurement_win);
    S.mod_power(j,1:P.N_mod_channels) = stat_mod_power_win(subband_envs(:,j),mod_filts,P.use_zp,measurement_win);
    S.mod_C2(j,1:size(mod_C2_filts,2)-1,1:2) = stat_mod_C2_win(subband_envs(:,j),mod_C2_filts,P.use_zp,measurement_win);
end

%compute subband envelope, modulation band correlations
S.env_C = stat_corr_filt_win_full(subband_envs,ones(size(mod_C1_filts(:,1))),P.use_zp,measurement_win);
S.mod_C1 = stat_corr_filt_win_full(subband_envs,mod_C1_filts,P.use_zp,measurement_win);


