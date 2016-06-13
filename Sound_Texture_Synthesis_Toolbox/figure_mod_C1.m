
figure('Position',[5 100 1912 1025]);

for f = 1:size(mod_C1_filts,2)
    subplot(3,4,2*f-1); imagesc(target_S.mod_C1(:,:,f));set(gca,'CLim',[-.2 1]);colorbar;axis square
    if f>2
        title(['Orig, ' num2str(C1_cfreqs(f)) ' Hz band']);
    end
    subplot(3,4,2*f); imagesc(synth_S.mod_C1(:,:,f));set(gca,'CLim',[-.2 1]);colorbar;axis square
    if f==1
        title(['C1 Comparison for ' fig_title_string],'FontSize',12);
    else
        title(['Synth, ' num2str(C1_cfreqs(f)) ' Hz band']);
    end
end

set(gcf,'PaperOrientation', 'landscape','PaperPosition',[0.25 0.25 10.5 8]);
