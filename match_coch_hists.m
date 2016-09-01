function coch_synth = match_coch_hists(coch_orig, coch_synth, ti)

% coch_synth = match_coch_hists(coch_orig, coch_synth, ti)
%
% histogram matches cochlear envelopes such that hist(coch_orig(:, i)) =
% hist(coch_synth(:, i)). the histogram matching is performed twice, once you
% all time points, and once you just the temporal indices specified by ti such
% that hist(coch_orig(ti, i)) = hist(coch_synth(ti, i))
%
% -- Example --
%
% coch_orig = randn(1000,10);
% coch_synth = randn(1000,10)*3;
% ti = 200:800;
%
% figure;
% hist([coch_orig(ti,5), coch_synth(ti,5)]);
%
% figure;
% coch_synth = match_coch_hists(coch_orig, coch_synth, ti);
% hist([coch_orig(ti,5), coch_synth(ti,5)]);

assert_equal_dimensions(coch_orig, coch_synth);

% match subset of teporal indices, if they're not equal to the full set
all_ti = 1:size(coch_orig,1);
match_subset = ~isempty(setdiff(all_ti, ti));

for j = 1:size(coch_synth,2)
    
    % match all temporal indices
    [~,xi] = sort(coch_synth(:,j));
    coch_synth(xi,j) = sort(coch_orig(:,j));
    
    % match subset of temporal indices
    if match_subset
        [~,xi] = sort(coch_synth(ti,j));
        coch_synth(ti(xi),j) = sort(coch_orig(ti,j));
    end
    
    %     % match outside indices
    %     ti_other = setdiff(1:size(coch_synth,1), ti);
    %     [~,xi] = sort(coch_synth(ti_other,j));
    %     coch_synth(ti_other(xi),j) = sort(coch_orig(ti_other,j));
    
end