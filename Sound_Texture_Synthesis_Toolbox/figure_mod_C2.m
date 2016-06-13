figure('Position',[5 100 1912 1025]);

for k = 1:16
    ind = k*2;
    subplot(4,4,k);

    semilogx(C2_cfreqs(1:end-1),target_S.mod_C2(ind,:,1),'b-o');hold on;semilogx(C2_cfreqs(1:end-1),target_S.mod_C2(ind,:,2),'b-x');
    semilogx(C2_cfreqs(1:end-1),synth_S.mod_C2(ind,:,1),'r-o');semilogx(C2_cfreqs(1:end-1),synth_S.mod_C2(ind,:,2),'r-x');
    set(gca,'YLim',[-1 1]);
    if rem(k,4)==1
        ylabel('C2 Correlation','FontSize',10);
    end
    if k>12
        xlabel('Mod Freq (Hz)','FontSize',10);
    end
    if k==2
        title(['Mod. C2 Comparison for ' fig_title_string],'FontSize',12);
    elseif k>3
        title(['Subband ' num2str(ind) '    (' num2str(round(audio_cutoffs_Hz(ind))) ' Hz)'],'FontSize',10);
    end
    if k==4
        legend('Real','Imag.','Location','Best');
    end
end
set(gcf,'PaperOrientation', 'landscape','PaperPosition',[0.25 0.25 10.5 8]);
