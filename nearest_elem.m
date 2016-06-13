function [y, ai] = nearest_elem(a, x)

% [y, ai] = nearest_elem(a, x)
% 
% Finds the element in array a that is nearest in euclidian distance to scalar x
% 
% -- Example --
% 
% [y,ai] = nearest_elem([3 2 4.5 9 10 7.03], 4.4)

[~, ai] = min(abs(a-x).^2);
y = a(ai);