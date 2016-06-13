function X = unwrap_allbutlast(X)

% unwrap all of the dimensions but the last in array X

dims = size(X);
X = reshape(X, [prod(dims(1:end-1)), dims(end)]);

