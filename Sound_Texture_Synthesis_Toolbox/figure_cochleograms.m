% compare cochleograms formed from subband envelopes

% Dec 2012 -- Josh McDermott <jhm@mit.edu>
% May 2014 -- modified to fix incompatibility with overcomplete filter
% banks -- Josh McDermott

figure('Position',[441 85 799 1015]);

N_sub = length(target_S.subband_var);
if P.use_more_audio_filters == 1
    subs_to_label = (N_sub-1):-5:1;
elseif P.use_more_audio_filters == 2
    subs_to_label = (N_sub-2):-5:1;
else
    subs_to_label = (N_sub):-5:1;
end
ytick_pos = 1:5:((length(subs_to_label)-1)*5+1);

if exist('upBlur','file')==2 %upBlur is in path

    subplot(2,1,1);
    imagesc(upBlur(flipud(orig_subband_envs'),1));
    title(['Orig ' fig_title_string],'FontSize',12);
    ylabel('Center Frequency (Hz)','FontSize',10);
    xlabel('Time (s)','FontSize',10);
    set(gca,'XTick',[1:floor(length(orig_subband_envs)/P.env_sr)]*2*P.env_sr,'XTickLabel',num2str([1:floor(length(orig_subband_envs)/P.env_sr)]'));
    set(gca,'YTick',ytick_pos*2,'YTickLabel',num2str(round(audio_cutoffs_Hz(subs_to_label))'));
    clim = get(gca,'CLim');
    if length(orig_sound)>=length(synth_sound)
        xlim = get(gca,'XLim');
    end
    
    
    subplot(2,1,2);
    imagesc(upBlur(flipud(synth_subband_envs'),1));
    title(['Synth ' fig_title_string],'FontSize',12);
    ylabel('Center Frequency (Hz)','FontSize',10);
    xlabel('Time (s)','FontSize',10);
    set(gca,'XTick',[1:floor(length(synth_subband_envs)/P.env_sr)]*2*P.env_sr,'XTickLabel',num2str([1:floor(length(synth_subband_envs)/P.env_sr)]'));
    set(gca,'YTick',ytick_pos*2,'YTickLabel',num2str(round(audio_cutoffs_Hz(subs_to_label))'));
    set(gca,'CLim',clim);
    if length(orig_sound)>=length(synth_sound) && length(orig_sound)<=1.5*length(synth_sound)
        set(gca,'XLim',xlim);
    else
        xlim=get(gca,'XLim');
        subplot(2,1,1);hold on;set(gca,'XLim',xlim);
    end
    colormap(flipud(gray));

else
    
    subplot(2,1,1);
    imagesc(flipud(orig_subband_envs'));
    title(['Orig ' fig_title_string],'FontSize',12);
    ylabel('Center Frequency (Hz)','FontSize',10);
    xlabel('Time (s)','FontSize',10);
    set(gca,'XTick',[1:floor(length(orig_subband_envs)/P.env_sr)]*P.env_sr,'XTickLabel',num2str([1:floor(length(orig_subband_envs)/P.env_sr)]'));
    set(gca,'YTick',ytick_pos,'YTickLabel',num2str(round(audio_cutoffs_Hz(subs_to_label))'));
    clim = get(gca,'CLim');
    if length(orig_sound)>=length(synth_sound)
        xlim = get(gca,'XLim');
    end
    
    
    subplot(2,1,2);
    imagesc(flipud(synth_subband_envs'));
    title(['Synth ' fig_title_string],'FontSize',12);
    ylabel('Center Frequency (Hz)','FontSize',10);
    xlabel('Time (s)','FontSize',10);
    set(gca,'XTick',[1:floor(length(synth_subband_envs)/P.env_sr)]*P.env_sr,'XTickLabel',num2str([1:floor(length(synth_subband_envs)/P.env_sr)]'));
    set(gca,'YTick',ytick_pos,'YTickLabel',num2str(round(audio_cutoffs_Hz(subs_to_label))'));
    set(gca,'CLim',clim);
    if length(orig_sound)>=length(synth_sound) && length(orig_sound)<=1.5*length(synth_sound)
        set(gca,'XLim',xlim);
    else
        xlim=get(gca,'XLim');
        subplot(2,1,1);hold on;set(gca,'XLim',xlim);
    end
    colormap(flipud(gray));
    
end

