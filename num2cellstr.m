function s = num2cellstr(x, varargin)

% s = num2cellstr(x, varargin)
% 
% applies num2str to each number of an array x, and returns a cell of the same
% dimension as x
% 
% -- Example --
% 
% x = randn(4,3)
% num2cellstr( x, '%.2f')

s = cell(size(x));
for i = 1:numel(x)
    s{i} = num2str(x(i), varargin{:});
end