% compare spectrograms computed with fft (linear frequency scale)

figure('Position',[441 85 799 1015]);

if length(orig_sound)>1.5*length(synth_sound)
    temp1=orig_sound(1:length(synth_sound));
    temp2=synth_sound;
elseif length(synth_sound)>1.5*length(orig_sound)
    temp2=synth_sound(1:length(orig_sound));
    temp1=orig_sound;
else
    temp1=orig_sound;
    temp2=synth_sound;
end

subplot(2,1,1);
j_specgram2(temp1,P.audio_sr,0);
title(['Orig ' fig_title_string],'FontSize',12);
ylabel('Frequency (Hz)','FontSize',10);
xlabel('Time (sec)','FontSize',10);
if length(temp1)>=length(temp2)
    xlim=get(gca,'XLim');
end

subplot(2,1,2);
j_specgram2(temp2,P.audio_sr,0);
title(['Synth ' fig_title_string],'FontSize',12);
ylabel('Frequency (Hz)','FontSize',10);
xlabel('Time (sec)','FontSize',10);
if length(temp1)>length(temp2)
    set(gca,'XLim',xlim);
else
    xlim=get(gca,'XLim');
    subplot(2,1,1);hold on;set(gca,'XLim',xlim);
end

