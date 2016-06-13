figure('Position',[5 100 1912 1025]);

log_ac_plot = 1;

%big_lag_set = [0:800]; %samples
biggest_spacing = P.env_ac_intervals_smp(end)-P.env_ac_intervals_smp(end-1);
big_lag_set = [P.env_ac_intervals_smp [P.env_ac_intervals_smp(end)+biggest_spacing : biggest_spacing : 800]];
[big_filt_set,~] = make_env_acc_filters2(length(orig_subband_envs)*2, P.env_sr, big_lag_set);
[big_filt_set2,~] = make_env_acc_filters2(length(synth_subband_envs)*2, P.env_sr, big_lag_set);

ones_win = ones(length(orig_subband_envs),1);
ones_win2 = ones(length(synth_subband_envs),1);

for k = 1:16
    ind = k*2;
    subplot(4,4,k);
    clear all_acc_values_synth all_acc_values_orig
    all_acc_values_synth = stat_env_ac_scaled_win(synth_subband_envs(:,ind),big_filt_set2,big_lag_set,ones_win2)./( (length(synth_subband_envs)-big_lag_set)/length(synth_subband_envs));
    all_acc_values_orig = stat_env_ac_scaled_win(orig_subband_envs(:,ind),big_filt_set,big_lag_set,ones_win)./( (length(orig_subband_envs)-big_lag_set)/length(orig_subband_envs));
    if log_ac_plot
        semilogx(big_lag_set/P.env_sr*1000,all_acc_values_orig,'b');hold on;semilogx(big_lag_set/P.env_sr*1000,all_acc_values_synth,'r');
        for l=1:length(P.env_ac_intervals_smp)
            semilogx(P.env_ac_intervals_smp(l)/P.env_sr*1000,all_acc_values_orig(find(big_lag_set==P.env_ac_intervals_smp(l))),'bo','MarkerSize',3);
            semilogx(P.env_ac_intervals_smp(l)/P.env_sr*1000,all_acc_values_synth(find(big_lag_set==P.env_ac_intervals_smp(l))),'ro','MarkerSize',3);
        end
        set(gca,'XLim',[big_lag_set(1)*.9 big_lag_set(end)/.9]);
    else
        plot(big_lag_set/P.env_sr*1000,all_acc_values_orig,'b');hold on;plot(big_lag_set/P.env_sr*1000,all_acc_values_synth,'r');
        for l=1:length(P.env_ac_intervals_smp)
            plot(P.env_ac_intervals_smp(l)/P.env_sr*1000,all_acc_values_orig(find(big_lag_set==P.env_ac_intervals_smp(l))),'bo','MarkerSize',3);
            plot(P.env_ac_intervals_smp(l)/P.env_sr*1000,all_acc_values_synth(find(big_lag_set==P.env_ac_intervals_smp(l))),'ro','MarkerSize',3);
        end
    end
    set(gca,'YLim',[-1 1]);
    if rem(k,4)==1
        ylabel('Autocorr','FontSize',10);
    end
    if k>12
        xlabel('Lag (ms)','FontSize',10);
    end
    if k==2
        title(['Env. Autocorr. Comparison for ' fig_title_string],'FontSize',12);
    elseif k>3
        title(['Subband ' num2str(ind) '    (' num2str(round(audio_cutoffs_Hz(ind))) ' Hz)'],'FontSize',10);
    end
end
set(gcf,'PaperOrientation', 'landscape','PaperPosition',[0.25 0.25 10.5 8]);
