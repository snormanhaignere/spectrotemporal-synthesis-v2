P = synthesis_parameters_default;
coch = randn(P.env_sr*P.max_duration_sec, 9*round(1/P.logf_spacing));
padded_and_depadded_coch = remove_pad(pad_coch(coch, P), P);
figure;
plot(coch(:), padded_and_depadded_coch(:));