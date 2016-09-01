function I = parse_optInputs_keyvalue(optargs, I)

% Parse optional inputs specified as a key-value pair. Key must be a string.
% 
% I is an optional input argument with a structure containing the default values
% of the parameters, which are overwritten by the optional arguments. If I is
% specified then every key must match one of the fields of the structure I. 
% 
% -- Example --
% 
% % parse key value pairs
% optargs = {'key1', {1 2 3}, 'key2', 1:3, 'key3', 'abc'};
% I = parse_optInputs_keyvalue(optargs)
% 
% % specifying default values
% I.key1 = {1 2 3};
% I.key2 = 1:3;
% I.key3 = 'abc';
% I = parse_optInputs_keyvalue({'key1', {4,5,6}}, I)
% 
% % use defaults to catch a spelling mistake
% I = parse_optInputs_keyvalue({'keys1', {4,5,6}}, I)
% 
% 2016-08-27: Created, Sam NH

% should be an event number of arguments
n_optargs = length(optargs);
assert(mod(n_optargs,2)==0);

% initialize with empty structure if not specified
if nargin < 2
    I = struct;
else % extract list of possible keys from the values of I if specified
    possible_keys = fieldnames(I);
end

% immediately return if there are no optional arguments
if n_optargs == 0
    return;
end


% assign keys and values
i_key = 1:2:n_optargs;
i_val = 2:2:n_optargs;
for j = 1:n_optargs/2
    key = optargs{i_key(j)};
    value = optargs{i_val(j)};

    % check key is a string
    if ~ischar(key)
        error('Optional arguments not formatted propertly\nAll keys must be strings\n');
    end
    
    % check key is one of the possible keys, and if it has the same type class
    % type
    if exist('possible_keys', 'var')
        if ~any(strcmp(key, possible_keys))
            error(['Optional arguments not formatted propertly\n' ...
                '''%s'' not a valid key\n'], key);
        end
        
        if ~isequal(class(I.(key)), class(value))
            allowed_class_swaps = {'double', 'int32'};
            allowed_swap = false;
            for k = 1:size(allowed_class_swaps,1)
                if strcmp(class(I.(key)), allowed_class_swaps{k,1}) ...
                        && strcmp(class(value), allowed_class_swaps{k,2})
                    allowed_swap = true;
                end
                
                if strcmp(class(value), allowed_class_swaps{k,1}) ...
                        && strcmp(class(I.(key)), allowed_class_swaps{k,2})
                    allowed_swap = true;
                end                
            end
            if ~allowed_swap
                error(['Optional arguments not formatted propertly\n' ...
                    'Value of ''%s'' should be of type %s\n'], key, class(I.(key)));
            end
        end
    end
    
    % assign
    I.(key) = value;
    
end