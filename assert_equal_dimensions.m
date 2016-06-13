function assert_equal_dimensions(X,Y)

% tests whether too matrices have the same size
% 
% assert_equal_dimensions(rand(4,3), rand(4,3))
% assert_equal_dimensions(rand(4,3), rand(4,2))

Z = size(X) == size(Y);
assert(all(Z(:)));