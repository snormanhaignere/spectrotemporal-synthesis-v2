function index = optInputs(argList, argTest)
% index = optInputs(argList, argTest)

index = 0;
for i = 1:length(argList)
    if isequal(argList{i}, argTest)
        if index == 0;
            index = i;
        else
            index = [index i]; %#ok<AGROW>
        end
    end
end