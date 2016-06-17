function coch_synth = ...
    match_coch_hists_every_Nsec(coch_orig, coch_synth, ti, P)

% coch_synth = match_coch_hists_every_Nsec(coch_orig, coch_synth, ti)
%
% histogram matches cochlear envelopes such that hist(coch_orig(:, i)) =
% hist(coch_synth(:, i)). the histogram matching is once for the entire
% cochleogram, and then a second time for short Nsec blocks of the cochleogram

assert_equal_dimensions(coch_orig, coch_synth);

% samples per segment and total number of segments
smps_per_seg = round(P.match_coch_every_Nsec * P.env_sr);
n_seg = ceil(length(ti) / smps_per_seg);

for j = 1:size(coch_synth,2)
    
    % match all temporal indices
    [~,xi] = sort(coch_synth(:,j));
    coch_synth(xi,j) = sort(coch_orig(:,j));
    
    % match segments
    for k = 1:n_seg
        
        % indices to match
        xi = (1:smps_per_seg) + smps_per_seg * (k-1);
        ti_for_this_seg = ti(xi);
        clear xi;
        
        % matching
        [~,xi] = sort(coch_synth(ti_for_this_seg,j));
        coch_synth(ti_for_this_seg(xi),j) = sort(coch_orig(ti_for_this_seg,j));
    end
    
end

