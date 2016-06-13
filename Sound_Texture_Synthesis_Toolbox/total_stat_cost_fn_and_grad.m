function [f,f_grad] = total_stat_cost_fn_and_grad(current_subband_env,...
    all_subband_envs,current_sub_num,P,target_S,C_bands_to_use,mod_C1_bands_to_use,...
    mod_filts,env_ac_filts, mod_C1_filts, mod_C2_filts, imp_win, constraints_to_impose)

% FUNCTION [F,F_GRAD] = TOTAL_STAT_COST_FN_AND_GRAD(CURRENT_SUBBAND_ENV,...
%    ALL_SUBBAND_ENVS,CURRENT_SUB_NUM,P,TARGET_S,C_BANDS_TO_USE,MOD_C1_BANDS_TO_USE,...
%    MOD_FILTS,ENV_AC_FILTS, MOD_C1_FILTS, MOD_C2_FILTS, IMP_WIN, CONSTRAINTS_TO_IMPOSE)
%
% returns the value of a cost function F based on a specified
% set of statistics, as well as the gradient F_GRAD of this cost function
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

[all_grads, error_vect] = env_stat_grads_and_errors(current_subband_env,...
    all_subband_envs,current_sub_num,P,target_S,C_bands_to_use,mod_C1_bands_to_use,...
    mod_filts,env_ac_filts, mod_C1_filts, mod_C2_filts, imp_win, constraints_to_impose);


f = sum(error_vect.^2);
f_grad = -2*sum((ones(length(current_subband_env),1)*error_vect').*all_grads,2);

