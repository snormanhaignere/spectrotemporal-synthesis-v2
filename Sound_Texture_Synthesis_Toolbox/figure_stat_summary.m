function figure_stat_summary(target_S, synth_S, P, orig_sound, synth_sound, ...
    audio_cutoffs_Hz, C1_cfreqs, fig_title_string)
%
% function FIGURE_STAT_SUMMARY(TARGET_S, SYNTH_S, P, ORIG_SOUND, SYNTH_SOUND, ...
%    AUDIO_CUTOFFS_HZ, C1_CFREQS, FIG_TITLE_STRING)
%
% makes figure with statistical comparison of original and synth sounds
%
% written to be called following synthesis procedure, from within
% RUN_SYNTHESIS function
%

%
% Dec 2012 -- Josh McDermott

figure('Position',[5 100 1912 1025]);
subplot(4,4,1);
[Pow,Freq]=pwelch(orig_sound,[],[],[],P.audio_sr);[Pow2,Freq2]=pwelch(synth_sound,[],[],[],P.audio_sr);
plot(Freq,10*log10(Pow)); hold on; plot(Freq2,10*log10(Pow2),'r');
set(gca,'XScale','log','XLim',[audio_cutoffs_Hz(1) audio_cutoffs_Hz(end)]);
ylabel('Power (dB/Hz)','FontSize',10);
xlabel('Frequency (Hz)','FontSize',10);

subplot(4,4,2);
semilogx(audio_cutoffs_Hz,target_S.env_mean,'b');hold on;
semilogx(audio_cutoffs_Hz,synth_S.env_mean,'r');
set(gca,'XLim',[audio_cutoffs_Hz(1)/2 audio_cutoffs_Hz(end)*2]);
ylabel('Env Mean','FontSize',10);
xlabel('Center Freq (Hz)','FontSize',10);
title(['Stat Comparison for ' fig_title_string],'FontSize',12);

subplot(4,4,3);
semilogx(audio_cutoffs_Hz,target_S.env_var,'b');hold on;
semilogx(audio_cutoffs_Hz,synth_S.env_var,'r');
set(gca,'XLim',[audio_cutoffs_Hz(1)/2 audio_cutoffs_Hz(end)*2]);
ylabel('Env StdDev/Mean','FontSize',10);
xlabel('Center Freq (Hz)','FontSize',10);

subplot(4,4,4);
semilogx(audio_cutoffs_Hz,target_S.env_skew,'b');hold on;
semilogx(audio_cutoffs_Hz,synth_S.env_skew,'r');
%%leave out lowpass and highpass end subbands
%semilogx(audio_cutoffs_Hz(2:end-1),target_S.env_skew(2:end-1),'b');hold on;
%semilogx(audio_cutoffs_Hz(2:end-1),synth_S.env_skew(2:end-1),'r');
set(gca,'XLim',[audio_cutoffs_Hz(1)/2 audio_cutoffs_Hz(end)*2]);
ylabel('Env Skewness','FontSize',10);
xlabel('Center Freq (Hz)','FontSize',10);
legend('Orig','Synth','Location','Best');

subplot(4,4,5);
imagesc(flipud(target_S.mod_power)); set(gca,'CLim',[0 max(target_S.mod_power(:))]); colorbar
title('Orig. Mod. Power','FontSize',10);
ylabel('Cochlear Channel','FontSize',10);set(gca,'YTickLabel',[]);xlabel('Mod. Channel','FontSize',10);

subplot(4,4,6);
imagesc(flipud(synth_S.mod_power)); set(gca,'CLim',[0 max(target_S.mod_power(:))]); colorbar
title('Synth. Mod. Power','FontSize',10);
ylabel('Cochlear Channel','FontSize',10);set(gca,'YTickLabel',[]);xlabel('Mod. Channel','FontSize',10);

subplot(4,4,7);
imagesc(flipud([target_S.mod_C2(:,:,1) target_S.mod_C2(:,:,2)])); set(gca,'CLim',[-1 1]); colorbar
hold on; plot([6.5 6.5],[0 33],'k')
title('Orig Mod. C2');
ylabel('Cochlear Channel');set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);xlabel('Real-Real Real-Imag');

subplot(4,4,8);
imagesc(flipud([synth_S.mod_C2(:,:,1) synth_S.mod_C2(:,:,2)])); set(gca,'CLim',[-1 1]); colorbar
hold on; plot([6.5 6.5],[0 33],'k')
title('Synth Mod. C2');
ylabel('Cochlear Channel');set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);xlabel('Real-Real Real-Imag');

subplot(4,4,9);
imagesc(flipud(target_S.env_C));axis square; colormap(jet); set(gca,'CLim',[-.5 1],'XTick',[],'YTick',[]); colorbar
title('Orig Env C','FontSize',10);
ylabel('Cochlear Channel');

subplot(4,4,10);
imagesc(flipud(synth_S.env_C));axis square; colormap(jet); set(gca,'CLim',[-.5 1],'XTick',[],'YTick',[]); colorbar
title('Synth Env C','FontSize',10);
xlabel('Cochlear Channel');

subplot(4,4,11);
imagesc(flipud(target_S.mod_C1(:,:,end-3)));axis square; colormap(jet); set(gca,'CLim',[-.5 1],'XTick',[],'YTick',[]); colorbar
title(['Orig Mod C1 (' num2str(C1_cfreqs(end-3),'%3.1f') ' Hz)'],'FontSize',10);

subplot(4,4,12);
imagesc(flipud(synth_S.mod_C1(:,:,end-3)));axis square; colormap(jet); set(gca,'CLim',[-.5 1],'XTick',[],'YTick',[]); colorbar
title(['Synth Mod C1 (' num2str(C1_cfreqs(end-3),'%3.1f') ' Hz)'],'FontSize',10);
ylabel('Cochlear Channel');

subplot(4,4,13);
imagesc(flipud(target_S.mod_C1(:,:,end-2)));axis square; colormap(jet); set(gca,'CLim',[-.5 1],'XTick',[],'YTick',[]); colorbar
title(['Orig Mod C1 (' num2str(C1_cfreqs(end-3),'%3.1f') ' Hz)'],'FontSize',10);
ylabel('Cochlear Channel');

subplot(4,4,14);
imagesc(flipud(synth_S.mod_C1(:,:,end-2)));axis square; colormap(jet); set(gca,'CLim',[-.5 1],'XTick',[],'YTick',[]); colorbar
title(['Synth Mod C1 (' num2str(C1_cfreqs(end-2),'%3.1f') ' Hz)'],'FontSize',10);

subplot(4,4,15);
imagesc(flipud(target_S.mod_C1(:,:,end-1)));axis square; colormap(jet); set(gca,'CLim',[-.5 1],'XTick',[],'YTick',[]); colorbar
title(['Orig Mod C1 (' num2str(C1_cfreqs(end-1),'%3.1f') ' Hz)'],'FontSize',10);
xlabel('Cochlear Channel');

subplot(4,4,16);
imagesc(flipud(synth_S.mod_C1(:,:,end-1)));axis square; colormap(jet); set(gca,'CLim',[-.5 1],'XTick',[],'YTick',[]); colorbar
title(['Synth Mod C1 (' num2str(C1_cfreqs(end-1),'%3.1f') ' Hz)'],'FontSize',10);

set(gcf,'PaperOrientation', 'landscape','PaperPosition',[0.25 0.25 10.5 8]);
