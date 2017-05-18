function coch_synth = match_coch_mean_and_var(coch_orig, coch_synth, ti)

% coch_synth = match_coch_mean_and_var(coch_orig, coch_synth, ti)
%
% matches mean and variance of the cochlear envelopes. the matching is performed
% twice, once using all time points, and once you just the temporal indices
% specified by ti such that moment(coch_orig(ti, i)) = moment(coch_synth(ti, i))
%
% -- Example --
%
% coch_orig = randn(1000,10)-1;
% coch_synth = randn(1000,10)*3+2;
% ti = 200:800;
%
% figure;
% hist([coch_orig(ti,5), coch_synth(ti,5)]);
% moment_measures(coch_orig, 1)
% moment_measures(coch_synth, 1)
%
% figure;
% coch_synth = match_coch_mean_and_var(coch_orig, coch_synth, ti);
% hist([coch_orig(ti,5), coch_synth(ti,5)]);
% moment_measures(coch_orig, 1)
% moment_measures(coch_synth, 1)
% 
% 2017-05-17: Created, Sam NH

assert_equal_dimensions(coch_orig, coch_synth);

% match all time points
coch_synth = match(coch_orig, coch_synth);

% match subset of indices
all_ti = 1:size(coch_orig,1);
if ~isempty(setdiff(all_ti, ti))
    coch_synth(ti, :) = match(coch_orig(ti, :), coch_synth(ti, :));
end

function coch_synth = match(coch_orig, coch_synth)

% demean
coch_synth = bsxfun(@minus, coch_synth, mean(coch_synth));

% normalize std/variance, leave channels with zero variance
std_synth = std(coch_synth);
std_synth(std_synth==0) = 1;
coch_synth = bsxfun(@times, coch_synth, 1./std_synth);

% match std/variance
std_orig = std(coch_orig);
coch_synth = bsxfun(@times, coch_synth, std_orig);

% match mean
coch_synth = bsxfun(@plus, coch_synth, mean(coch_orig));

