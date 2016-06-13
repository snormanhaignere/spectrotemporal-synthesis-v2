log_ac_plot=0;

figure('Position',[5 100 1912 1025]);

undo_win=1;
win_choice=2;

%calculate autocorr of noise subbands for comparison
[x,n_sr]=wavread('pink_noise_20s_20kHz.wav');
if n_sr ~= P.audio_sr
    x = resample(x, P.audio_sr, n_sr);
    n_sr = P.audio_sr;
end
if length(x)/P.audio_sr > 10
    x = x(1:10*P.audio_sr);
end

if P.lin_or_log_filters==1 || P.lin_or_log_filters==2
    if P.use_more_audio_filters==0
        [n_audio_filts, n_audio_cutoffs_Hz] = make_erb_cos_filters(length(x), P.audio_sr, P.N_audio_channels, P.low_audio_f, P.hi_audio_f);
    elseif P.use_more_audio_filters==1
        [n_audio_filts, n_audio_cutoffs_Hz] = make_erb_cos_filts_double2(length(x), P.audio_sr, P.N_audio_channels, P.low_audio_f, P.hi_audio_f);
    elseif P.use_more_audio_filters==2
        [n_audio_filts, n_audio_cutoffs_Hz] = make_erb_cos_filts_quadruple2(length(x), P.audio_sr, P.N_audio_channels, P.low_audio_f, P.hi_audio_f);
    end
elseif P.lin_or_log_filters==3 || P.lin_or_log_filters==4
    if P.use_more_audio_filters==0
        [n_audio_filts, n_audio_cutoffs_Hz] = make_lin_cos_filters(length(x), P.audio_sr, P.N_audio_channels, P.low_audio_f, P.hi_audio_f);
    elseif P.use_more_audio_filters==1
        [n_audio_filts, n_audio_cutoffs_Hz] = make_lin_cos_filts_double(length(x), P.audio_sr, P.N_audio_channels, P.low_audio_f, P.hi_audio_f);
    end
end

n_subbands = generate_subbands(x, n_audio_filts);

noise_sub_ac = autocorr_mult_zp(n_subbands, P.sub_ac_win_choice, P.sub_ac_undo_win);

orig_sub_ac = autocorr_mult_zp(orig_subbands, P.sub_ac_win_choice, P.sub_ac_undo_win);

synth_sub_ac = autocorr_mult_zp(synth_subbands, P.sub_ac_win_choice, P.sub_ac_undo_win);


L2 = length(orig_subbands);
C2 = L2/2+1;
Ls = length(synth_subbands);
Cs = Ls/2+1;
Ln = length(n_subbands);
Cn = Ln/2+1;

sub_ac_N_smp = round(P.num_sub_ac_period./n_audio_cutoffs_Hz*P.audio_sr);
sub_ac_N_smp(sub_ac_N_smp > P.num_sub_ac_period/20*P.audio_sr)=P.num_sub_ac_period/20*P.audio_sr;

for k = 1:16
    ind = k*2;
    subplot(4,4,k);
    if log_ac_plot
        semilogx([0:sub_ac_N_smp(ind)*2]/P.audio_sr*1000,orig_sub_ac(C2:C2+sub_ac_N_smp(ind)*2,ind)/orig_sub_ac(C2,ind),'b');hold on;
        semilogx([0:sub_ac_N_smp(ind)*2]/P.audio_sr*1000,synth_sub_ac(Cs:Cs+sub_ac_N_smp(ind)*2,ind)/synth_sub_ac(Cs,ind),'r');
        semilogx([0:sub_ac_N_smp(ind)*2]/P.audio_sr*1000,noise_sub_ac(Cn:Cn+sub_ac_N_smp(ind)*2,ind)/noise_sub_ac(Cn,ind),'g');
    else
        plot([0:sub_ac_N_smp(ind)*2]/P.audio_sr*1000,orig_sub_ac(C2:C2+sub_ac_N_smp(ind)*2,ind)/orig_sub_ac(C2,ind),'b');hold on;
        plot([0:sub_ac_N_smp(ind)*2]/P.audio_sr*1000,synth_sub_ac(Cs:Cs+sub_ac_N_smp(ind)*2,ind)/synth_sub_ac(Cs,ind),'r');
        plot([0:sub_ac_N_smp(ind)*2]/P.audio_sr*1000,noise_sub_ac(Cn:Cn+sub_ac_N_smp(ind)*2,ind)/noise_sub_ac(Cn,ind),'g');
    end
    set(gca,'YLim',[-1 1]);
    ylabel('Autocorr.','FontSize',10);
    xlabel('Lag (ms)','FontSize',10);
    if k==2
        title(['Subband Autocorrelation Comparison for ' fig_title_string],'FontSize',12);
    elseif k>3
        title(['Subband ' num2str(ind) '    (' num2str(round(audio_cutoffs_Hz(ind))) ' Hz)'],'FontSize',10);
    end
    if k==4
        legend('Orig','Synth','Noise','Location','Best');
    end
end
set(gcf,'PaperOrientation', 'landscape','PaperPosition',[0.25 0.25 10.5 8]);
