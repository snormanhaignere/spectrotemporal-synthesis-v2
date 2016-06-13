%filenaming stuff
constraint_text=[];

stat_text = {'sub_var','sub_kurt','env_mean','env_var','env_skew','env_kurt','env_C','env_ac','mod_pow','mod_C1','mod_C2'};
for stat = 1:length(stat_text)
    if eval(['P.constraint_set.' stat_text{stat}]) && eval(['P.use_noise_stats.' stat_text{stat}])
        constraint_text = [constraint_text num2str(3)];
    else
        constraint_text = [constraint_text num2str(eval(['P.constraint_set.' stat_text{stat}]))];
    end
end

new_filename = ['_' constraint_text];

if P.avg_stat_option==0
    temp = P.orig_sound_filename;
    if strcmp(temp(end-3:end), '.wav')
        new_filename = [temp(1:end-4) new_filename];
    else
        new_filename = [temp new_filename];
    end
elseif P.avg_stat_option==1
    new_filename = ['AVG_' P.avg_filename new_filename];
elseif P.avg_stat_option==2
    new_filename = ['MORPH_' P.avg_filename '_' num2str(P.morph_ratio) new_filename];
end

