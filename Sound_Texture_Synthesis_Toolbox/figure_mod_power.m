figure('Position',[5 100 1912 1025]);

for k = 1:16
    ind = k*2;
    subplot(4,4,k);
    semilogx(mod_cfreqs_Hz,target_S.mod_power(ind,:));hold on;
    semilogx(mod_cfreqs_Hz,synth_S.mod_power(ind,:),'r');
    set(gca,'YLim',[0 max(max(target_S.mod_power))*1.1]);
    set(gca,'XLim',[mod_cfreqs_Hz(1)/2 mod_cfreqs_Hz(end)*2]);
    if rem(k,4)==1
        ylabel('Channel Power/Total Power','FontSize',10);
    end
    if k>12
        xlabel('Mod Freq (Hz)','FontSize',10);
    end
    if k==2
        title(['Mod. Power Comparison for ' fig_title_string],'FontSize',12);
    elseif k>3
        title(['Subband ' num2str(ind) '    (' num2str(round(audio_cutoffs_Hz(ind))) ' Hz)'],'FontSize',10);
    end
end
set(gcf,'PaperOrientation', 'landscape','PaperPosition',[0.25 0.25 10.5 8]);
