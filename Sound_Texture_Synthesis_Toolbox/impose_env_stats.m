function [new_subband_envs] = impose_env_stats(old_subband_envs, P, target_S, ...
    mod_filts, env_ac_filts, mod_C1_filts, mod_C2_filts, ...
    imp_win, num_line_searches, adjustment_prop, constraints_to_impose)
%
% [NEW_SUBBAND_ENVS] = IMPOSE_ENV_STATS(OLD_SUBBAND_ENVS, P, TARGET_S, ...
%    MOD_FILTS, ENV_AC_FILTS, MOD_C1_FILTS, MOD_C2_FILTS, ...
%    IMP_WIN, NUM_LINE_SEARCHES, ADJUSTMENT_PROP, CONSTRAINTS_TO_IMPOSE)
%
% loops through OLD_SUBBAND_ENVS, computes gradients and errors for each
% statistic flagged in CONSTRAINTS_TO_IMPOSE and adjusts each subband to
% simultaneously impose all specified statistics, either with conjugate
% gradient descent or the Gaus-Newton method
%
% needs modulation filters for calculation of statistic gradients
%
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
% May 2014 -- modified to fix incompatibility with overcomplete filter
% banks -- Josh McDermott

N_sub = length(target_S.subband_var);

if P.sub_imposition_order==1 %same order each iteration
    %find subband with most power, work outwards
    [~,starting_sub] = max(target_S.env_mean);
elseif P.sub_imposition_order==2 %somewhat stochastic
    %randomly choose subband in top 20, work outwards
    [~,subs_ordered_by_power] = sort(target_S.env_mean,1,'descend');
    starting_sub = subs_ordered_by_power(ceil(rand*20));
end    
%set order to go through subbands
sub_order = starting_sub;
low_subs = starting_sub-1:-1:1; l=1;
high_subs = starting_sub+1:N_sub; h=1;
while length(sub_order)<N_sub
    if l<=length(low_subs)
        sub_order = [sub_order low_subs(l)];
        l=l+1;
    end
    if h<=length(high_subs)
        sub_order = [sub_order high_subs(h)];
        h=h+1;
    end
end


new_subband_envs = old_subband_envs;
for k=1:length(sub_order)
    j=sub_order(k);
    if constraints_to_impose(7)
        potential_bands=P.C_offsets_to_impose+j; %correlations to be enforced
        C_bands_to_use=[];
        for b=1:length(potential_bands)
            if length(find(potential_bands(b)==sub_order(1:k-1)))>0
                C_bands_to_use = [C_bands_to_use potential_bands(b)];
            end
        end
    else
        C_bands_to_use=[];
    end
    if constraints_to_impose(10)
        potential_bands=P.mod_C1_offsets_to_impose+j; %correlations to be enforced
        mod_C1_bands_to_use=[];
        for b=1:length(potential_bands)
            if length(find(potential_bands(b)==sub_order(1:k-1)))>0
                mod_C1_bands_to_use = [mod_C1_bands_to_use potential_bands(b)];
            end
        end
    else
        mod_C1_bands_to_use=[];
    end

    if length(C_bands_to_use)>0 || length(mod_C1_bands_to_use)>0 || sum(constraints_to_impose([3:6 8 9 11]))>0
        current_sub_num=j;
                
        if P.imposition_method == 1 %conjugate gradient
            
            %checkgrad('total_stat_cost_fn_and_grad', new_subband_envs(:,current_sub_num), 10^-5,...
            %    new_subband_envs, current_sub_num, P, target_S, C_bands_to_use, mod_C1_bands_to_use,...
            %    mod_filts,env_ac_filts, mod_C1_filts, mod_C2_filts, imp_win, constraints_to_impose);
            
            [updated_subband_env, cost_fn_history, ~] = minimize(new_subband_envs(:,current_sub_num), 'total_stat_cost_fn_and_grad', num_line_searches,...
                new_subband_envs, current_sub_num, P, target_S, C_bands_to_use, mod_C1_bands_to_use,...
                mod_filts,env_ac_filts, mod_C1_filts, mod_C2_filts, imp_win, constraints_to_impose);            

            if P.compression_option==0 || P.compression_option==1
                updated_subband_env(updated_subband_env<0)=0; %force envelope to be positive
            elseif P.compression_option==2
                updated_subband_env(updated_subband_env<P.log_constant)=P.log_constant;
            end
            
            %fprintf('Subband %2.0f (%2.0f): Starting log error = %3.2f; Ending log error = %3.2f...\n',k,j,log10(cost_fn_history(1)),log10(cost_fn_history(end)));
            log_error_start(j) = log10(cost_fn_history(1)); log_error_end(j) = log10(cost_fn_history(end));
            new_subband_envs(:,j) = updated_subband_env;
        elseif P.imposition_method == 2
            for loop=1:num_line_searches %misnomer - these are not line searches but rather repeats of the alternative optimization procedure
                [all_grads, error_vect] = env_stat_grads_and_errors(new_subband_envs(:,current_sub_num),new_subband_envs,current_sub_num,...
                    P,target_S,C_bands_to_use,mod_C1_bands_to_use,...
                    mod_filts,env_ac_filts, mod_C1_filts, mod_C2_filts, imp_win);
                [U,S,V]=svd(all_grads','econ');
                if length(S)<size(all_grads,2)
                    warning(['dimension of S: ' num2str(length(S)) ' is less than number of constraints (' num2str(size(all_grads,2)) ')']);
                end
                adjustment = V*inv(S)*U'*error_vect;
                if length(adjustment)==length(new_subband_envs)
                    new_subband_envs(:,current_sub_num) = new_subband_envs(:,current_sub_num) + adjustment_prop*adjustment;
                end
            end
        end
    end
end
if P.imposition_method == 1 && P.check_imposition_success == 1
    fprintf('Mean starting log error = %3.2f; Mean ending log error = %3.2f...\n',mean(log_error_start),mean(log_error_end));
end



