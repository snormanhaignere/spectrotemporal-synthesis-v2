function M = moment_measures(X, DIM)

% M = moment_measures(X, DIM)
% 
% calculates the mean, variance, skew, and kurtosis of matrix X along dimension
% DIM. 
% 
% -- Example --
% X = randn(3,10000);
% moment_measures(X, 2)

dims = size(X);
dims = [dims,1];

% reshaping function
dims_reduced = dims( [1:DIM-1, DIM+1:end] );
fn_reshape = @(X)reshape(X, dims_reduced);

% concatenate reshaped moment arrays
M = cat(...
    ndims(X), ...
    fn_reshape( mean(X, DIM) ), ...
    fn_reshape( var(X, [], DIM) ), ...
    fn_reshape( skewness(X, [], DIM) ), ...
    fn_reshape( kurtosis(X, [], DIM) ));