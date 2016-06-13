function [synth_sound] = run_synthesis(P)
%
% [SYNTH_SOUND] = RUN_SYNTHESIS(P)
%
% Generates a synthetic signal that matches the statistics of a WAV file
% whose name is specified in P.ORIG_SOUND_FILENAME.
%
% Must be passed a set of parameters in P that can be set to their default
% values by the SYNTHESIS_PARAMETERS script.
%
% Please see the README_sound_texture file for more information.
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
%
% modified Nov 2013 to allow initialization with arbitrary stored waveform.
% Nov 2013 -- Josh McDermott <jhm@mit.edu>


if P.avg_stat_option==0
    fprintf('Starting %s\n\n',P.orig_sound_filename);
elseif P.avg_stat_option==1
    fprintf('Starting %s\n\n',[P.avg_filename ' - AVERAGE']);
elseif P.avg_stat_option==2
    fprintf('Starting %s\n\n',[P.avg_filename ' - MORPH with ratio of ' num2str(P.morph_ratio)]);
end
tic

%COMPUTE STATISTICS FROM SAMPLE
if P.avg_stat_option==0
    orig_sound = format_orig_sound(P);
    measurement_win = set_measurement_window(length(orig_sound),P.measurement_windowing,P);
    target_S = measure_texture_stats(orig_sound,P,measurement_win);
    target_S = edit_measured_stats(target_S,P);
    %generate subbands and envelopes of original sound for future use
    [orig_subbands,orig_subband_envs] = generate_subbands_and_envs(orig_sound, P.audio_sr, ...
        P.env_sr, P.N_audio_channels, P.low_audio_f, P.hi_audio_f, P.lin_or_log_filters, ...
        P.use_more_audio_filters,P.compression_option,P.comp_exponent,P.log_constant);
elseif P.avg_stat_option>=1
    [target_S, waveform_avg] = measure_texture_stats_avg(P);
    target_S = edit_measured_stats(target_S,P);
    [orig_subbands,orig_subband_envs] = generate_subbands_and_envs(waveform_avg, P.audio_sr, ...
        P.env_sr, P.N_audio_channels, P.low_audio_f, P.hi_audio_f, P.lin_or_log_filters, ...
        P.use_more_audio_filters,P.compression_option,P.comp_exponent,P.log_constant);
    orig_sound = waveform_avg;
    measurement_win = set_measurement_window(length(orig_sound),P.measurement_windowing,P);
end
format_filename;

%GENERATE FILTERS FOR SYNTHESIS
ds_factor=P.audio_sr/P.env_sr; %factor by which envelopes are downsampled
synth_dur_smp = ceil(P.desired_synth_dur_s*P.audio_sr/ds_factor)*ds_factor; %ensures that length in samples is an integer multiple of envelope sr
P.length_ratio = synth_dur_smp/length(orig_sound);
%make audio filters
if P.lin_or_log_filters==1 || P.lin_or_log_filters==2
    if P.use_more_audio_filters==0
        [audio_filts, audio_cutoffs_Hz] = make_erb_cos_filters(synth_dur_smp, P.audio_sr, P.N_audio_channels, P.low_audio_f, P.hi_audio_f);
    elseif P.use_more_audio_filters==1
        [audio_filts, audio_cutoffs_Hz] = make_erb_cos_filts_double2(synth_dur_smp, P.audio_sr, P.N_audio_channels, P.low_audio_f, P.hi_audio_f);
    elseif P.use_more_audio_filters==2
        [audio_filts, audio_cutoffs_Hz] = make_erb_cos_filts_quadruple2(synth_dur_smp, P.audio_sr, P.N_audio_channels, P.low_audio_f, P.hi_audio_f);
    end
elseif P.lin_or_log_filters==3 || P.lin_or_log_filters==4
    if P.use_more_audio_filters==0
        [audio_filts, audio_cutoffs_Hz] = make_lin_cos_filters(synth_dur_smp, P.audio_sr, P.N_audio_channels, P.low_audio_f, P.hi_audio_f);
    elseif P.use_more_audio_filters==1
        [audio_filts, audio_cutoffs_Hz] = make_lin_cos_filts_double(synth_dur_smp, P.audio_sr, P.N_audio_channels, P.low_audio_f, P.hi_audio_f);
    end
end
N_sub = size(audio_filts,2);

%make modulation filters
if P.use_zp==1
    mod_filt_length = synth_dur_smp/ds_factor*2;
elseif P.use_zp==0
    mod_filt_length = synth_dur_smp/ds_factor;
end
[env_ac_filts,ac_filt_cutoffs] = make_env_acc_filters2(mod_filt_length, P.env_sr, P.env_ac_intervals_smp);
if P.lin_or_log_filters==1 || P.lin_or_log_filters==3
    [mod_filts,mod_cfreqs_Hz,mod_freqs] = make_constQ_cos_filters(mod_filt_length, P.env_sr, P.N_mod_channels, P.low_mod_f, P.hi_mod_f, P.mod_filt_Q_value);
elseif P.lin_or_log_filters==2 || P.lin_or_log_filters==4
    [mod_filts,mod_cfreqs_Hz,mod_freqs] = make_lin_cos_constQ_cntrl_filters(mod_filt_length, P.env_sr, P.N_mod_channels, P.low_mod_f, P.hi_mod_f, P.mod_filt_Q_value);
end
[mod_C2_filts,C2_cfreqs,C2_freqs] = make_octave_cos_filters2(mod_filt_length, P.env_sr, P.low_mod_f_C12, P.hi_mod_f);
if P.lin_or_log_filters==1 || P.lin_or_log_filters==3
    mod_C1_filts = mod_C2_filts(:,2:end);
    C1_cfreqs = C2_cfreqs(2:end);
elseif P.lin_or_log_filters==2 || P.lin_or_log_filters==4
    [mod_C1_filts,C1_cfreqs,C1_freqs] = make_lin_cos_oct_cntrl_filters(mod_filt_length, P.env_sr, size(mod_C2_filts,2)-1, C2_cfreqs(2), P.hi_mod_f);
end

sub_ac_N_smp = round(P.num_sub_ac_period./audio_cutoffs_Hz*P.audio_sr);
sub_ac_N_smp(sub_ac_N_smp > P.num_sub_ac_period/20*P.audio_sr)=P.num_sub_ac_period/20*P.audio_sr;

imp_win = set_measurement_window(synth_dur_smp, P.imposition_windowing, P);

%INITIALIZE SYNTH SOUND
if P.initialize_with_sound_file == 0
    synth_sound = randn(synth_dur_smp,1); %create Gaussian noise signal
    synth_sound = synth_sound/rms(synth_sound)*P.desired_rms;
elseif P.initialize_with_sound_file == 1
     [temp,sr] = wavread([P.initial_sound_folder P.initial_sound_filename]);
     if size(temp,2)==2
         temp = temp(:,1); %turn stereo files into mono
     end
     if sr ~= P.audio_sr
         temp = resample(temp, P.audio_sr, sr);
     end
     if length(temp) < synth_dur_smp
         error('File with which to initialize synthesis is not long enough!\n');
     else
         synth_sound = temp(1:synth_dur_smp);
     end
end

clear *error* *snr* *adj_mag

%SYNTHESIS LOOP
for it_n=1:P.num_it
    fprintf('starting iteration %d...\n',it_n);
    
    %save interim versions of synth signal
    if ismember(it_n,P.iteration_snapshots)
        wavwrite(synth_sound, P.audio_sr, [P.output_folder new_filename '_it' num2str(it_n) '.wav']);
    end
    
    %the incremental imposition is mostly not used in this version of the
    %code - onyl relevant for one optimization method
    if P.incremental_imposition == 1
        if it_n<P.incr_imp_it_limit
            adjustment_prop = it_n/P.incr_imp_it_limit;
        else
            adjustment_prop = 1;
        end
    else
        adjustment_prop=1;
    end
    
    %gradually increase number of line searches in optimization, or not
    if P.num_LS_option==1
        num_line_searches = P.init_LS_num+round(.2*(it_n-1));
    elseif P.num_LS_option==0
        num_line_searches = 20;
    end
    
    %GENERATE SUBBANDS AND THEIR ENVELOPES
    synth_subbands = generate_subbands(synth_sound, audio_filts);
    synth_analytic_subbands = hilbert(synth_subbands);
    synth_subband_envs = abs(synth_analytic_subbands);
    synth_subband_phases = synth_analytic_subbands./synth_subband_envs;
    %compress envelopes
    if P.compression_option==1
        synth_subband_envs = synth_subband_envs.^P.comp_exponent;
    elseif P.compression_option==2
        synth_subband_envs = log10(synth_subband_envs+log_constant);
    end
    %downsample envelopes, storing high-pass residual (important to
    %avoiding spectrum artifacts)
    env_temp1 = synth_subband_envs;
    synth_subband_envs = resample(synth_subband_envs,1,ds_factor);
    synth_subband_envs(synth_subband_envs<0)=0;
    env_temp2 = resample(synth_subband_envs,ds_factor,1);
    env_residual = env_temp1 - env_temp2;
    
    %COMPUTE DIFFERENCE BETWEEN SYNTH STATISTICS AND TARGET STATISTICS
    synth_S = measure_texture_stats(synth_sound,P,imp_win,audio_filts, ...
        audio_cutoffs_Hz, mod_filts, env_ac_filts, mod_C1_filts, mod_C2_filts, ...
        target_S.subband_bins, target_S.env_bins);
    if it_n>1
        SNRs = compute_stat_SNRs(target_S,synth_S,P,SNRs,it_n);
    else
        SNRs = compute_stat_SNRs(target_S,synth_S,P,[],it_n);
    end
    %check if snrs are sufficient
    
    all_current_snrs = [SNRs.subband_var(it_n) SNRs.subband_kurt(it_n) ...
        SNRs.env_mean(it_n) SNRs.env_var(it_n) SNRs.env_skew(it_n) SNRs.env_kurt(it_n) ...
        SNRs.env_C(it_n) SNRs.env_ac(it_n) SNRs.mod_power(it_n) SNRs.mod_C1(it_n) SNRs.mod_C2(it_n)];
    
    constraints_to_impose = [P.constraint_set.sub_var P.constraint_set.sub_kurt ...
        P.constraint_set.env_mean P.constraint_set.env_var P.constraint_set.env_skew P.constraint_set.env_kurt ...
        P.constraint_set.env_C P.constraint_set.env_ac P.constraint_set.mod_pow P.constraint_set.mod_C1 P.constraint_set.mod_C2];
    
    %halt synthesis loop if all statistics are imposed to criterion SNR
    if length(find(all_current_snrs(constraints_to_impose>0)>=P.end_criterion_db)) == length(find(constraints_to_impose>0))
        break
    else %don't impose statistics that are already sufficiently close to desired values
        temp_constraints_to_impose=constraints_to_impose;
        if P.omit_converged_stats
            temp_constraints_to_impose(all_current_snrs > P.omission_criterion_db) = 0;
        end
    end
    
    if sum(temp_constraints_to_impose(3:end))>0
        [new_envs] = impose_env_stats(synth_subband_envs, P, target_S, ...
            mod_filts, env_ac_filts, mod_C1_filts, mod_C2_filts, ...
            imp_win, num_line_searches, adjustment_prop, temp_constraints_to_impose);
    else
        new_envs = synth_subband_envs;
    end
    
    if P.check_imposition_success
        env_S = measure_texture_stats_from_envs(new_envs, P, imp_win, ...
            mod_filts, env_ac_filts, mod_C1_filts, mod_C2_filts, target_S.env_bins);
        display_before_after_error(target_S,synth_S,env_S,P);
    end
    
    if P.manual_mean_var_adjustment
        for k=1:size(filts,2)
            new_envs(:,k) = new_envs(:,k)-mean(new_envs(:,k));
            new_envs(:,k) = new_envs(:,k)/std(new_envs(:,k))*subband_mag_var(k)*subband_mag_mean(k);
            new_envs(:,k) = new_envs(:,k)+subband_mag_mean(k);
        end
    end
    
    if P.match_env_hist==1 %match full histogram
        for k=1:size(filts,2)
            new_envs(:,k) = j_histoMatch(new_envs(:,k), target_S.env_hist(k,:), target_S.env_bins(k,:));
        end
    end
    
    new_envs = real(new_envs); %make sure everything is real
    
    %get rid of nan values
    if P.compression_option==2
        new_envs(isnan(new_envs))=log10(P.log_constant);
    else
        new_envs(isnan(new_envs))=0;
    end
    
    
    %upsample envelope back to full resolution
    new_envs = resample(new_envs,ds_factor,1);
    new_envs = new_envs + env_residual;
    new_envs(new_envs<0)=0;
    if P.compression_option==1
        new_envs = new_envs.^(1/P.comp_exponent);
    elseif P.compression_option==2
        new_envs = 10.^new_envs-P.log_constant;
    end
    
    %combine new magnitudes with old phases
    new_synth_analytic_subbands = synth_subband_phases.*new_envs;
    synth_subbands = real(new_synth_analytic_subbands);
    
    %match subband marginals
    for k=1:N_sub
        if P.constraint_set.sub_var>0 || P.constraint_set.sub_kurt>0
            synth_subbands(:,k) = synth_subbands(:,k)-mean(synth_subbands(:,k));
            if P.constraint_set.sub_var==1
                if var(synth_subbands(:,k))>P.log_constant
                    synth_subbands(:,k) = synth_subbands(:,k)*sqrt(target_S.subband_var(k)/var(synth_subbands(:,k)));
                end
            end
            if P.constraint_set.sub_kurt==1
                if SNRs.subband_kurt(it_n) < P.leave_out_convergence_criterion_db
                    if exist('modkurt','file')==2
                        [synth_subbands(:,k)] = modkurt(synth_subbands(:,k), target_S.subband_kurt(k),adjustment_prop);
                    end
                end
            end
            synth_subbands(:,k) = synth_subbands(:,k)+target_S.subband_mean(k);
        end
        if P.match_sub_hist
            synth_subbands(:,k) = j_histoMatch(synth_subbands(:,k), target_S.subband_hist(k,:), target_S.subband_bins(k,:));
        end
    end
    synth_subbands(isnan(synth_subbands))=0;
    
    %resynthesize
    synth_sound = collapse_subbands(synth_subbands, audio_filts);
    
end

fprintf('Final SNRs:\n');
synth_S = measure_texture_stats(synth_sound,P,imp_win,audio_filts, ...
    audio_cutoffs_Hz, mod_filts, env_ac_filts, mod_C1_filts, mod_C2_filts, ...
    target_S.subband_bins, target_S.env_bins);
it_n=it_n+1;
SNRs = compute_stat_SNRs(target_S,synth_S,P,SNRs,it_n);
write_final_SNRs_and_params_to_file(P,SNRs,new_filename);

%for comparison, noise with spectrum of original sound
if P.write_spec_match_noise
    noise = spectrally_matched_noise(orig_sound);
    wavwrite(noise, P.audio_sr, [P.output_folder 'sm_noise_' P.orig_sound_filename]);
end
%generate copy of original with same normalization as synth sound
if P.write_norm_orig
    if P.avg_stat_option==0
        if length(orig_sound)>1.5*length(synth_sound)
            wavwrite(orig_sound(1:length(synth_sound)), P.audio_sr, [P.output_folder 'norm_' P.orig_sound_filename]);
        else
            wavwrite(orig_sound, P.audio_sr, [P.output_folder 'norm_' P.orig_sound_filename]);
        end
    elseif P.avg_stat_option>=1 %for comparison, average of waveforms
        wavwrite(waveform_avg, P.audio_sr, [P.output_folder 'wvfm_avg_' P.avg_filename]);
    end
end

%save synthetic result
if P.imposition_windowing==2 %imposition is not circular, so leave off ends where window faded to zero
    wavwrite(synth_sound(P.audio_sr:end-P.audio_sr), sr, [P.output_folder new_filename '.wav']);
else
    wavwrite(synth_sound, P.audio_sr, [P.output_folder new_filename '.wav']);
    if P.write_repeat
        wavwrite([synth_sound;synth_sound], P.audio_sr, [P.output_folder new_filename '_rep.wav']);
    end
end

toc

if P.avg_stat_option==0
    fig_file_prefix =  [P.orig_sound_filename(1:end-4) '_' constraint_text];
    fig_title_string = P.orig_sound_filename(1:end-4);
    fig_title_string(fig_title_string=='_') = ' ';
elseif P.avg_stat_option==1
    fig_file_prefix =  [P.avg_filename '_AVERAGE' '_' constraint_text];
    fig_title_string = ['Average ' P.avg_filename];
    fig_title_string(fig_title_string=='_') = ' ';
elseif P.avg_stat_option==2
    fig_file_prefix =  [P.avg_filename '_MORPH_' num2str(P.morph_ratio) '_' constraint_text];
    fig_title_string = ['Morph ' P.avg_filename];
    fig_title_string(fig_title_string=='_') = ' ';
end

if P.display_figures
    
    %compute subband envelopes once more for use in figures
    synth_subbands = generate_subbands(synth_sound, audio_filts);
    synth_subband_envs = abs(hilbert(synth_subbands));
    if P.compression_option==1
        synth_subband_envs = synth_subband_envs.^P.comp_exponent;
    elseif P.compression_option==2
        synth_subband_envs = log10(synth_subband_envs+P.log_constant);
    end
    synth_subband_envs = resample(synth_subband_envs,1,ds_factor);
    synth_subband_envs(synth_subband_envs<0)=0;
    
    figure_snr_plots;
    if P.save_figures
        saveas(gcf, [P.output_folder fig_file_prefix '_SNRs.tif']);
        close
    end
    
    figure_stat_summary(target_S, synth_S, P, orig_sound, synth_sound, ...
        audio_cutoffs_Hz, C1_cfreqs, fig_title_string);
    if P.save_figures
        saveas(gcf, [P.output_folder fig_file_prefix '_stat_summ.jpeg']);
        close
    end
    
    figure_specgrams;
    if P.save_figures
        saveas(gcf, [P.output_folder fig_file_prefix '_sgrams.jpeg']);
        close
    end
    
    figure_cochleograms;
    if P.save_figures
        saveas(gcf, [P.output_folder fig_file_prefix '_cgrams.jpeg']);
        close
    end
    
    if P.N_audio_channels>=15
        figure_env_hists;
        if P.save_figures
            saveas(gcf, [P.output_folder fig_file_prefix '_env_hists.jpeg']);
            close
        end
        
        figure_sub_hists;
        if P.save_figures
            saveas(gcf, [P.output_folder fig_file_prefix '_sub_hists.jpeg']);
            close
        end
        
        figure_sub_ac;
        if P.save_figures
            saveas(gcf, [P.output_folder fig_file_prefix '_sub_ac.jpeg']);
            close
        end
        
        figure_mod_power;
        if P.save_figures
            saveas(gcf, [P.output_folder fig_file_prefix '_MP.jpeg']);
            close
        end
        
        %%plots the full modulation spectrum of a subset of envelopes - not
        %%guaranteed to be matched between synthetic and original sound 
        %figure_mod_spectrum_full;
        %if P.save_figures
        %    saveas(gcf, [P.output_folder fig_file_prefix '_mod_spec.jpeg']);
        %    close
        %end
        
        figure_mod_C1;
        if P.save_figures
            saveas(gcf, [P.output_folder fig_file_prefix '_mod_C1.jpeg']);
            close
        end
        
        figure_mod_C2;
        if P.save_figures
            saveas(gcf, [P.output_folder fig_file_prefix '_mod_C2.jpeg']);
            close
        end
        
        figure_env_ac;
        if P.save_figures
            saveas(gcf, [P.output_folder fig_file_prefix '_env_ac.jpeg']);
            close
        end
    end
end


