% FUNCTION [TARGET_S, WAVEFORM_AVG] = MEASURE_TEXTURE_STATS_AVG(P)
%
% returns target statistics TARGET_S that are the average of the statistics
% of a set of sounds specified in P.FILES_FOR_STAT_AVG (either the first two
% files, weighted by P.MORPH_RATIO, or all the files in equal proportion)
%
% WAVEFORM_AVG is a weighted average of the sound waveforms
%

% Dec 2012 -- Josh McDermott <jhm@mit.edu>

function [target_S, waveform_avg] = measure_texture_stats_avg(P)

sub_hist_bins = ones(P.N_audio_channels+2,1)*[-.035: .07/(P.n_hist_bins-1) : .035];
env_hist_bins = ones(P.N_audio_channels+2,1)*[0: .35/(P.n_hist_bins-1) : .35];
    
if P.avg_stat_option==1
    n_files = length(P.files_for_stat_avg);
elseif P.avg_stat_option==2
    n_files=2; %take first two files in P.files_for_stat_avg
end

for f=1:n_files
    orig_sound = format_orig_sound(P,f);

    if P.write_norm_orig
        if length(orig_sound)>1.5*length(P.desired_synth_dur_s)
            wavwrite(orig_sound(1:P.desired_synth_dur_s*P.audio_sr), P.audio_sr, [P.output_folder 'norm_' P.files_for_stat_avg{f}]);
        else
            wavwrite(orig_sound, P.audio_sr, [P.output_folder 'norm_' P.files_for_stat_avg{f}]);
        end
    end
    
    measurement_win = set_measurement_window(length(orig_sound),P.measurement_windowing,P);
    
    if f==1        
        %temp_S = measure_texture_stats(orig_sound,P,measurement_win);
        temp_S = measure_texture_stats(orig_sound,P,measurement_win, ... %use fixed set of bins
            [], [], [], [], [],[], sub_hist_bins, env_hist_bins);
        if P.avg_stat_option==1
            weight = 1/n_files;
        elseif P.avg_stat_option==2
            weight = P.morph_ratio;
        end
        waveform_avg = orig_sound * weight;
        
        target_S.subband_mean = temp_S.subband_mean * weight;
        target_S.subband_var = temp_S.subband_var * weight;
        target_S.subband_skew = temp_S.subband_skew * weight;
        target_S.subband_kurt = temp_S.subband_kurt * weight;
        target_S.env_mean = temp_S.env_mean * weight;
        target_S.env_var = temp_S.env_var * weight;
        target_S.env_skew = temp_S.env_skew * weight;
        target_S.env_kurt = temp_S.env_kurt * weight;
        target_S.subband_hist = temp_S.subband_hist * weight;
        %target_S.subband_bins = temp_S.subband_bins; %use bins from first sound for subsequent sounds, in order to add histograms
        target_S.subband_bins = sub_hist_bins;
        target_S.env_hist = temp_S.env_hist * weight;
        %target_S.env_bins = temp_S.env_bins;  %use bins from first sound for subsequent sounds, in order to add histograms
        target_S.env_bins = env_hist_bins;

        target_S.env_C = temp_S.env_C * weight;
        target_S.env_ac = temp_S.env_ac * weight;
        target_S.mod_power = temp_S.mod_power * weight;
        target_S.mod_C1 = temp_S.mod_C1 * weight;
        target_S.mod_C2 = temp_S.mod_C2 * weight;
        for k=1:P.N_audio_channels+2
            target_S.subband_ac{k} = temp_S.subband_ac{k} * weight;
        end
        
    elseif f>1
        temp_S = measure_texture_stats(orig_sound,P,measurement_win, ... %use fixed set of bins
            [], [], [], [], [],[], sub_hist_bins, env_hist_bins);
        
        if P.avg_stat_option==1
            weight = 1/n_files;
        elseif P.avg_stat_option==2
            weight = 1-P.morph_ratio;
        end
        waveform_avg = waveform_avg + orig_sound * weight;
        
        target_S.subband_mean = target_S.subband_mean + temp_S.subband_mean * weight;
        target_S.subband_var = target_S.subband_var + temp_S.subband_var * weight;
        target_S.subband_skew = target_S.subband_skew + temp_S.subband_skew * weight;
        target_S.subband_kurt = target_S.subband_kurt + temp_S.subband_kurt * weight;
        target_S.env_mean = target_S.env_mean + temp_S.env_mean * weight;
        target_S.env_var = target_S.env_var + temp_S.env_var * weight;
        target_S.env_skew = target_S.env_skew + temp_S.env_skew * weight;
        target_S.env_kurt = target_S.env_kurt + temp_S.env_kurt * weight;
        target_S.subband_hist = target_S.subband_hist + temp_S.subband_hist * weight;
        target_S.env_hist = target_S.env_hist + temp_S.env_hist * weight;
        
        target_S.env_C = target_S.env_C + temp_S.env_C * weight;
        target_S.env_ac = target_S.env_ac + temp_S.env_ac * weight;
        target_S.mod_power = target_S.mod_power + temp_S.mod_power * weight;
        target_S.mod_C1 = target_S.mod_C1 + temp_S.mod_C1 * weight;
        target_S.mod_C2 = target_S.mod_C2 + temp_S.mod_C2 * weight;
        for k=1:P.N_audio_channels+2
            target_S.subband_ac{k} = target_S.subband_ac{k} + temp_S.subband_ac{k} * weight;
        end        
        
    end
end

for k=1:P.N_audio_channels+2
    target_S.subband_ac_power(k) = sum(target_S.subband_ac{k}.^2); %used in SNR calculation
end


