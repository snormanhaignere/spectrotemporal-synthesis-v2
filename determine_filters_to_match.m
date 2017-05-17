function P = determine_filters_to_match(P)

P.temp_mod_to_match = [];
P.spec_mod_to_match = [];

% add temporal modulation filters
if P.match_temp_mod
    P.temp_mod_to_match = ...
        [P.temp_mod_to_match, P.temp_mod_rates];
    P.spec_mod_to_match = ...
        [P.spec_mod_to_match, nan(1,length(P.temp_mod_rates))];
end

% add spectral modulation filters
if P.match_spec_mod
    P.temp_mod_to_match = ...
        [P.temp_mod_to_match, nan(1,length( P.spec_mod_rates))];
    P.spec_mod_to_match = ...
        [P.spec_mod_to_match, P.spec_mod_rates];
end

% add spectrotemporal filters
if P.match_spectemp_mod
    for i = 1:length(P.temp_mod_rates)
        for j = 1:length(P.spec_mod_rates);
            if P.temp_mod_rates(i) == 0 || P.spec_mod_rates(j) == 0
                P.temp_mod_to_match = ...
                    [P.temp_mod_to_match, P.temp_mod_rates(i)];
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

% add separate low rate filters
% modulated in time but with a flat spectrum
if (P.match_temp_mod || P.match_spectemp_mod) ...
        && ~isempty(P.lowrate_tempfilts_flat_spec);
    P.temp_mod_to_match = ...
        [P.temp_mod_to_match, P.lowrate_tempfilts_flat_spec];
    P.spec_mod_to_match = ...
        [P.spec_mod_to_match, zeros(1,length(P.lowrate_tempfilts_flat_spec))];
end

% add separate low rate filters
% modulated in time but with an impulse spectrum
if (P.match_temp_mod || P.match_spectemp_mod) ...
        && ~isempty(P.lowrate_tempfilts_impulse_spec);
    P.temp_mod_to_match = ...
        [P.temp_mod_to_match, P.lowrate_tempfilts_impulse_spec];
    P.spec_mod_to_match = ...
        [P.spec_mod_to_match, nan(1,length(P.lowrate_tempfilts_impulse_spec))];
end

% check that the two vectors match in size
assert(length(P.temp_mod_to_match) == length(P.spec_mod_to_match));