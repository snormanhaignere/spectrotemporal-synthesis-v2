function X = squeeze_dim(X, DIM)

% Remove selected singleton dimension from array X. Safer than squeeze because
% only one dimension is removed and the function checks that the dimension is a
% singleton.
% 
% -- Example --
% X = randn(1,1,4);
% size(squeeze(X))
% size(squeeze_dim(X,2))
% squeeze_dim(X,3)

dims = [size(X),1];
assert(dims(DIM) == 1);
X = reshape(X, dims([1:DIM-1, DIM+1:end]));
