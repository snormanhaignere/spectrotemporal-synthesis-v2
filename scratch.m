P = synthesis_parameters_default;

% add spectrotemporal filters
P.temp_mod_to_match = [];
P.spec_mod_to_match = [];
if P.match_spectemp_mod
    for i = 1:length(P.temp_mod_rates)
        for j = 1:length(P.spec_mod_rates);
            if P.temp_mod_rates(i) == 0
                P.temp_mod_to_match = ...
                    [P.temp_mod_to_match, 0];
                P.spec_mod_to_match = ...
                    [P.spec_mod_to_match, P.spec_mod_rates(j)];
            else
                P.temp_mod_to_match = ...
                    [P.temp_mod_to_match, ...
                    P.temp_mod_rates(i), -P.temp_mod_rates(i)];
                P.spec_mod_to_match = ...
                    [P.spec_mod_to_match, ...
                    P.spec_mod_rates(j) * ones(1,2)];
            end
        end
    end
end

% test reconstruction (i.e. match cochleogram to itself)
coch = randn(P.env_sr*4, 9*round(1/P.logf_spacing));
coch_recon = ...
    match_filtcoch_hists(coch, coch, P, 1:P.env_sr*4);
figure;
plot(coch, coch_recon);