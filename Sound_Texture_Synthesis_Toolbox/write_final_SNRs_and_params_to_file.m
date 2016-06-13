% WRITE_FINAL_SNRS_AND_PARAMS_TO_FILE(P, SNRs, FILENAME)
%
% Writes the final SNRs from the synthesis process to a text file (whose
% name is generated from the FILENAME argument) along with the synthesis
% parameters.
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

function [] = write_final_SNRs_and_params_to_file(P,SNRs,filename)

all_current_snrs = [SNRs.subband_var(end) SNRs.subband_kurt(end) ...
    SNRs.env_mean(end) SNRs.env_var(end) SNRs.env_skew(end) SNRs.env_kurt(end) ...
    SNRs.env_C(end) SNRs.env_ac(end) SNRs.mod_power(end) SNRs.mod_C1(end) SNRs.mod_C2(end)];

constraints_to_impose = [P.constraint_set.sub_var P.constraint_set.sub_kurt ...
    P.constraint_set.env_mean P.constraint_set.env_var P.constraint_set.env_skew P.constraint_set.env_kurt ...
    P.constraint_set.env_C P.constraint_set.env_ac P.constraint_set.mod_pow P.constraint_set.mod_C1 P.constraint_set.mod_C2];

if length(find(all_current_snrs(constraints_to_impose>0)>=P.flag_criterion_db))~=length(find(constraints_to_impose>0))
    fid = fopen([P.output_folder filename '_NOT_CONVERGED.txt'],'wt');
    fprintf(fid,[filename '\n']);    
else
    fid = fopen([P.output_folder filename '_final_snrs.txt'],'wt');
    fprintf(fid,[filename '\n']);
end

fprintf(fid,'sub hist snr = %3.2f; sub kurt snr = %3.2f; sub ac snr = %3.2f; env hist snr = %3.2f\n',...
    SNRs.subband_hist(end), SNRs.subband_kurt(end), SNRs.subband_ac(end), SNRs.env_hist(end));
fprintf(fid,'env mean snr = %3.2f; env var snr = %3.2f; env skew snr = %3.2f; env kurt snr = %3.2f\n',...
    SNRs.env_mean(end), SNRs.env_var(end), SNRs.env_skew(end), SNRs.env_kurt(end));
fprintf(fid,'env C snr = %3.2f; env C all snr = %3.2f; mod power snr = %3.2f; env ac snr = %3.2f\n',...
    SNRs.env_C(end), SNRs.env_C_all(end), SNRs.mod_power(end), SNRs.env_ac(end));
fprintf(fid,'mod C1 snr = %3.2f; mod C1 all snr = %3.2f; mod C2 snr = %3.2f\n\n',...
    SNRs.mod_C1(end), SNRs.mod_C1_all(end), SNRs.mod_C2(end));


%write other parameters
fprintf(fid,'SYNTHESIS PARAMETERS:\n');

fprintf(fid, 'P.constraint_set.sub_var  = %d\n', P.constraint_set.sub_var); 
fprintf(fid, 'P.constraint_set.sub_kurt = %d\n', P.constraint_set.sub_kurt); 
fprintf(fid, 'P.constraint_set.env_mean = %d\n', P.constraint_set.env_mean); 
fprintf(fid, 'P.constraint_set.env_var  = %d\n', P.constraint_set.env_var); 
fprintf(fid, 'P.constraint_set.env_skew = %d\n', P.constraint_set.env_skew); 
fprintf(fid, 'P.constraint_set.env_kurt = %d\n', P.constraint_set.env_kurt); 
fprintf(fid, 'P.constraint_set.env_C    = %d\n', P.constraint_set.env_C); 
fprintf(fid, 'P.constraint_set.env_ac   = %d\n', P.constraint_set.env_ac); 
fprintf(fid, 'P.constraint_set.mod_pow  = %d\n', P.constraint_set.mod_pow); 
fprintf(fid, 'P.constraint_set.mod_C1   = %d\n', P.constraint_set.mod_C1); 
fprintf(fid, 'P.constraint_set.mod_C2   = %d\n', P.constraint_set.mod_C2); 
fprintf(fid,'\n');

fprintf(fid, 'desired_synth_dur = %3.2f; length_ratio = %3.2f; max_orig_dur = %3.2f;\n', ...
    P.desired_synth_dur_s, P.length_ratio, P.max_orig_dur_s);
fprintf(fid,'measurement_windowing = %d (1--> unwindowed, 2--> global window)\n', P.measurement_windowing);
fprintf(fid,'imposition_windowing = %d (1--> unwindowed, 2--> global window)\n', P.imposition_windowing);
fprintf(fid,'win_steepness = %3.2f\n', P.win_steepness);
fprintf(fid, 'imposition_method = %d (1--> conj grad; 2--> gauss-newton)\n', P.imposition_method);
fprintf(fid,'sub_imposition_order = %d\n', P.sub_imposition_order);
fprintf(fid,'\n');
fprintf(fid,'avg_stat_option = %d (1 = average stats across examples; 2 = take weighted average of stats for two sounds)\n', P.avg_stat_option);
fprintf(fid,'\n');
fprintf(fid,'Max number of iterations = %d\n',P.num_it); 
fprintf(fid,'SNR at which to halt imposition = %d dB\n',P.end_criterion_db);
fprintf(fid,'SNR at to flag synthesis as failed = %d dB\n', P.flag_criterion_db);
fprintf(fid,'omit_converged_stats = %d; SNR at which to omit statistics from synthesis = %d dB\n', ...
    P.omit_converged_stats, P.omission_criterion_db);
fprintf(fid,'power threshold for including subbands in snr calculation = -%d dB\n', P.sig_sub_cutoff_db);
fprintf(fid,'incremental_imposition = %d; incr. imposition iteration limit = %d\n', P.incremental_imposition, P.incr_imp_it_limit);
fprintf(fid,'initial num line searches = %d; line search num option = %d\n', P.init_LS_num, P.num_LS_option);
fprintf(fid,'\n');
fprintf(fid,'audio_sr = %d\n', P.audio_sr);
fprintf(fid,'N_audio_channels = %d; low_audio_f = %d Hz; hi_audio_f = %d Hz\n', P.N_audio_channels, P.low_audio_f, P.hi_audio_f);
fprintf(fid,'use_more_audio_filters = %d\n', P.use_more_audio_filters);
fprintf(fid,'lin_or_log_filters = %d (1--> log audio & mod; 2--> log audio, lin mod; 3--> lin audio, log mod; 4--> lin audio, lin mod)\n', P.lin_or_log_filters);
fprintf(fid,'\n');
fprintf(fid,'env_sr = %d\n', P.env_sr);
fprintf(fid,'N_mod_channels = %d; low_mod_f = %3.2f Hz; hi_mod_f = %3.2f Hz\n', P.N_mod_channels, P.low_mod_f, P.hi_mod_f);
fprintf(fid,'mod_filt_Q_value = %d\n', P.mod_filt_Q_value);
fprintf(fid,'use_more_mod_filters = %d\n', P.use_more_mod_filters);
fprintf(fid,'low_mod_f_C12 = %3.2f Hz\n', P.low_mod_f_C12);
fprintf(fid,'use_zp = %d (0 means circular convolution; 1 means zeropadding for modulation filtering)\n', P.use_zp);
fprintf(fid,'\n');
fprintf(fid,'compression option = %d; comp_exponent = %2.2f; log_constant = %3.3f\n', P.compression_option, P.comp_exponent, P.log_constant);
fprintf(fid,'desired_rms = %2.2f\n', P.desired_rms);
fprintf(fid,'\n');

fprintf(fid,'match full envelope histograms = %d\n', P.match_env_hist);
fprintf(fid,'match full subband histograms = %d\n', P.match_sub_hist);
fprintf(fid,'num histogram bins = %d\n', P.n_hist_bins);
fprintf(fid,'manual envelope mean/variance adjustment = %d\n', P.manual_mean_var_adjustment);
fprintf(fid,'\n');

%fprintf(fid,'num_C_offsets = %d\n', P.num_C_offsets);
fprintf(fid,'C_offsets_to_impose = [');
for c = 1:length(P.C_offsets_to_impose)
    fprintf(fid,'%d ', P.C_offsets_to_impose(c));
end
fprintf(fid,']\n');
fprintf(fid,'\n');
fprintf(fid,'mod_C1_offsets_to_impose = [');
for c = 1:length(P.mod_C1_offsets_to_impose)
    fprintf(fid,'%d ', P.mod_C1_offsets_to_impose(c));
end
fprintf(fid,']\n');
fprintf(fid,'\n');
fprintf(fid,'lags at which to measure env autocorr (samples) =\n');
fprintf(fid,'     [');
for c = 1:length(P.env_ac_intervals_smp)
    fprintf(fid,'%d ', P.env_ac_intervals_smp(c));
end
fprintf(fid,']\n');
fprintf(fid,'\n');
fprintf(fid,'num periods of subband cf over which to match sub ac = %d\n', P.num_sub_ac_period);
fprintf(fid,'sub_ac_undo_win = %d\n', P.sub_ac_undo_win);
fprintf(fid,'sub_ac_win_choice = %d\n', P.sub_ac_win_choice);

% use_noise_stats = [P.use_noise_stats.sub_var P.use_noise_stats.sub_kurt P.use_noise_stats.sub_ac ...
%     P.use_noise_stats.env_mean P.use_noise_stats.env_var P.use_noise_stats.env_skew P.use_noise_stats.env_kurt ...
%     P.use_noise_stats.env_C P.use_noise_stats.env_ac P.use_noise_stats.mod_pow P.use_noise_stats.mod_C1 P.use_noise_stats.mod_C2];
% fprintf(fid,'Noise Stats Used: ');
% for c=1:length(use_noise_stats)
%     fprintf(fid,'%s',num2str(use_noise_stats(c)));
% end

%impose values of noise stats for those set to 1
fprintf(fid,'\n');
fprintf(fid, 'P.use_noise_stats.sub_var  = %d\n', P.use_noise_stats.sub_var); 
fprintf(fid, 'P.use_noise_stats.sub_kurt = %d\n', P.use_noise_stats.sub_kurt); 
fprintf(fid, 'P.use_noise_stats.env_mean = %d\n', P.use_noise_stats.env_mean); 
fprintf(fid, 'P.use_noise_stats.env_var  = %d\n', P.use_noise_stats.env_var); 
fprintf(fid, 'P.use_noise_stats.env_skew = %d\n', P.use_noise_stats.env_skew); 
fprintf(fid, 'P.use_noise_stats.env_kurt = %d\n', P.use_noise_stats.env_kurt); 
fprintf(fid, 'P.use_noise_stats.env_C    = %d\n', P.use_noise_stats.env_C); 
fprintf(fid, 'P.use_noise_stats.env_ac   = %d\n', P.use_noise_stats.env_ac); 
fprintf(fid, 'P.use_noise_stats.mod_pow  = %d\n', P.use_noise_stats.mod_pow); 
fprintf(fid, 'P.use_noise_stats.mod_C1   = %d\n', P.use_noise_stats.mod_C1); 
fprintf(fid, 'P.use_noise_stats.mod_C2   = %d\n', P.use_noise_stats.mod_C2); 


fprintf(fid,'\n');
fprintf(fid,'neg_env_skew = %d; neg_mod_C2 = %d\n', P.neg_env_skew, P.neg_mod_C2);

fclose(fid);

