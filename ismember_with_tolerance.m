function z = ismember_with_tolerance(x, y, tol)

% z = ismember_with_tolerance(x, y, tol)
% 
% like the built-in function ismember, but can be used with numbers that are the
% same up to some small tolerance
% 
% ismember([2 3 4 2], [2,3]+1e-15)
% ismember_with_tolerance([2 3 4 2], [2,3]+1e-15)
% ismember_with_tolerance([2 3 4 2], [2,3]+1e-6)
% ismember_with_tolerance([2 3 4 2], [2,3]+1e-6, 1e-3)

if nargin < 3
    tol = 1e-10;
end

z = false(size(x));

for i = 1:length(x)
    z(i) = any(abs(x(i) - y) < tol);
end