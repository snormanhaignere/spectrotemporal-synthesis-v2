function coch_synth_matched = match_filtcoch_covariance(coch_orig, coch_synth, P, ti)

%% Subbands

% filtered cochleograms
filtcoch_orig = coch2filtcoch_allsubbands(coch_orig, P);
filtcoch_synth = coch2filtcoch_allsubbands(coch_synth, P);

filtcoch_orig_tomatch = filtcoch_orig(ti,:,:);
filtcoch_synth_tomatch = filtcoch_synth(ti,:,:);

%% Match correlation structure

% unwrap into matrix
assert(all(size(filtcoch_orig_tomatch) == size(filtcoch_synth_tomatch)));
dims = size(filtcoch_orig_tomatch);
filtcoch_orig_tomatch = reshape(filtcoch_orig_tomatch, [dims(1), prod(dims(2:end))]);
filtcoch_synth_tomatch = reshape(filtcoch_synth_tomatch, [dims(1), prod(dims(2:end))]);

% compute svd of each filtered cochleogram
[U_orig, S_orig, V_orig] = svd(filtcoch_orig_tomatch, 'econ'); %#ok<ASGLU>
[U_synth, S_synth, V_synth] = svd(filtcoch_synth_tomatch, 'econ'); %#ok<ASGLU>

% impose covariance
% whiten -> U_synth * V_synth'
% project onto PCA bases of original -> * V_orig
% expand bases with singular values -> * S_orig
% project back -> * V_orig'
filtcoch_synth_matched = (U_synth * (V_synth' * V_orig) * S_orig) * V_orig';

% reshape, separating out cochleograms of different filters
filtcoch_synth_matched = reshape(filtcoch_synth_matched, dims);

% collapse subbands
filtcoch_synth_matched_withpadding = filtcoch_synth;
filtcoch_synth_matched_withpadding(ti,:,:) = filtcoch_synth_matched;
coch_synth_matched = filtcoch2coch(filtcoch_synth_matched_withpadding, P);

% figure;
% subplot(3,1,1);
% plot_cochleogram(coch_orig, P.f, P.t)
% subplot(3,1,2);
% plot_cochleogram(coch_synth, P.f, P.t)
% subplot(3,1,3);
% plot_cochleogram(coch_synth_match, P.f, P.t)
