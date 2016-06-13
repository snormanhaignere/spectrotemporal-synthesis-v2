function cellstr = append_str_to_each_cellstr(str, cellstr, front_or_back)
% appends a string to the front or back of each string in a cell array

for i = 1:length(cellstr)
    switch front_or_back
        case 'front'
            cellstr{i} = [str cellstr{i}];
        case 'back'
            cellstr{i} = [cellstr{i} str];
        otherwise
            error('front_or_back cannot be %s\n', front_or_back);
    end
end