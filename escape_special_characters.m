function str = escape_special_characters(str)

% place a backslash before all special characters in a string
% special characters are ' ', '(' and ')'

special_chars = {' ', '(', ')'};
for i = 1:length(special_chars);
    str = strrep(str, special_chars{i}, ['\' special_chars{i}]);
end