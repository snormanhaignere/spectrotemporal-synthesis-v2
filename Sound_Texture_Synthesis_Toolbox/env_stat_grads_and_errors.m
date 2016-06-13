function [all_grads, error_vect] = env_stat_grads_and_errors(current_subband_env, all_subband_envs, current_sub_num,...
    P,target_S, C_bands_to_use, mod_C1_bands_to_use,...
    mod_filts, env_ac_filts, mod_C1_filts, mod_C2_filts, imp_win, constraints_to_impose)

% [ALL_GRADS, ERROR_VECT] = ENV_STAT_GRADS_AND_ERRORS(CURRENT_SUBBAND_ENV, ALL_SUBBAND_ENVS, CURRENT_SUB_NUM,...
%    P,TARGET_S, C_BANDS_TO_USE, MOD_C1_BANDS_TO_USE,...
%    MOD_FILTS, ENV_AC_FILTS, MOD_C1_FILTS, MOD_C2_FILTS, IMP_WIN, CONSTRAINTS_TO_IMPOSE)
%
% returns a matrix of gradients and the associated errors for a set of
% statistics associated with a subband
%
% windowing is used to compute statistics and their gradients (though in
% practice IMP_WIN may be a vector of ones) 
%
%

% This code is part of an instantiation of a sound texture synthesis
% algorithm developed with Eero Simoncelli and described in this paper:
%
% McDermott, J.H., Simoncelli, E.P. (2011) Sound texture perception via
% statistics of the auditory periphery: Evidence from sound synthesis.
% Neuron, 71, 926-940. 
%
% Dec 2012 -- Josh McDermott <jhm@mit.edu>

k=current_sub_num;
all_grads=[];
error_vect=[];

%gradients are column vectors

%envelope moments
if constraints_to_impose(3) || constraints_to_impose(4) || constraints_to_impose(5) || constraints_to_impose(6)
    temp_mean = stat_central_moment_win(current_subband_env, 1, imp_win); %used in moment calculations
end
if constraints_to_impose(3)
    all_grads = [all_grads grad_central_moment_win(current_subband_env,1,imp_win,temp_mean)];
    error_vect = [error_vect; (target_S.env_mean(k) - stat_central_moment_win(current_subband_env, 1, imp_win, temp_mean))];
end
if constraints_to_impose(4)
    all_grads = [all_grads grad_central_moment_win(current_subband_env,2,imp_win,temp_mean)];
    error_vect = [error_vect; (target_S.env_var(k) - stat_central_moment_win(current_subband_env, 2, imp_win, temp_mean))];
end
if constraints_to_impose(5)
    all_grads = [all_grads grad_central_moment_win(current_subband_env,3,imp_win,temp_mean)];
    error_vect = [error_vect; (target_S.env_skew(k) - stat_central_moment_win(current_subband_env, 3, imp_win, temp_mean))];
end
if constraints_to_impose(6)
    all_grads = [all_grads grad_central_moment_win(current_subband_env,4,imp_win,temp_mean)];
    error_vect = [error_vect; (target_S.env_kurt(k) - stat_central_moment_win(current_subband_env, 4, imp_win, temp_mean))];
end

if constraints_to_impose(7) %cross-channel envelope correlation
    if length(C_bands_to_use)>0
        env_C_error_vect = [];
        env_C_grads=[];
        for f_off = 1:length(C_bands_to_use)
            env_C_value = stat_corr_win(current_subband_env, all_subband_envs(:,C_bands_to_use(f_off)), imp_win);
            env_C_error_vect = [env_C_error_vect; (target_S.env_C(k,C_bands_to_use(f_off)) - env_C_value)];
            env_C_grad = grad_corr_win(current_subband_env, all_subband_envs(:,C_bands_to_use(f_off)), env_C_value, imp_win);
            env_C_grads = [env_C_grads env_C_grad];
        end
        all_grads = [all_grads env_C_grads];
        error_vect = [error_vect; env_C_error_vect];
    end
end

if constraints_to_impose(8) %env autocorr
    ac_grads = grad_env_ac_scaled_win(current_subband_env, env_ac_filts, P.env_ac_intervals_smp, P.use_zp, imp_win);
    current_ac_values = stat_env_ac_scaled_win(current_subband_env, env_ac_filts, P.env_ac_intervals_smp, P.use_zp, imp_win);
    ac_error_vect = (target_S.env_ac(k,:) - current_ac_values)';
    all_grads = [all_grads ac_grads];
    error_vect = [error_vect; ac_error_vect];
end

if constraints_to_impose(9) %modulation power
    current_mp_values = stat_mod_power_win(current_subband_env, mod_filts, P.use_zp, imp_win);
    mp_error_vect = (target_S.mod_power(k,:) - current_mp_values)';
    mp_grads=[];
    for chan = 1:size(mod_filts,2)
        mp_grads = [mp_grads grad_mod_power_win(current_subband_env, mod_filts(:,chan), P.use_zp, imp_win)];
    end
    all_grads = [all_grads mp_grads];
    error_vect = [error_vect; mp_error_vect];
end

if constraints_to_impose(10) %C1 correlation
    if length(mod_C1_bands_to_use)>0
        mod_C1_error_vect = [];
        mod_C1_grads=[];
        mod_C1_values = stat_corr_filt_win(current_subband_env, all_subband_envs(:,mod_C1_bands_to_use),...
            mod_C1_filts, P.use_zp, imp_win);
        for f_off = 1:length(mod_C1_bands_to_use)
            n_C1_bands = size(mod_C1_filts,2);
            for b = 1:n_C1_bands
                mod_C1_value = mod_C1_values(f_off,b);
                mod_C1_error_vect = [mod_C1_error_vect; (target_S.mod_C1(k,mod_C1_bands_to_use(f_off),b) - mod_C1_value)];
                mod_C1_grad = grad_corr_filt_win(current_subband_env, all_subband_envs(:,mod_C1_bands_to_use(f_off)),...
                    mod_C1_value, mod_C1_filts(:,b), P.use_zp, imp_win);
                mod_C1_grads = [mod_C1_grads mod_C1_grad];
            end
        end
        all_grads = [all_grads mod_C1_grads];
        error_vect = [error_vect; mod_C1_error_vect];
    end
end

if constraints_to_impose(11) %C2 correlation
    current_mod_C2_values = stat_mod_C2_win(current_subband_env, mod_C2_filts, P.use_zp, imp_win);
    temp = squeeze(target_S.mod_C2(k,:,:)) - current_mod_C2_values;
    mod_C2_error_vect = temp(:);
    mod_C2_grads=[];
    for which_part=1:2 %real, imaginary
        for chan = 1:size(mod_C2_filts,2)-1
            mod_C2_grads = [mod_C2_grads grad_mod_C2_win(current_subband_env, mod_C2_filts(:,chan:chan+1), which_part, P.use_zp, imp_win)];
        end
    end
    all_grads = [all_grads mod_C2_grads];
    error_vect = [error_vect; mod_C2_error_vect];
end



