function outputfiles = mydir(inputdir, regexp_to_select, regexp_to_remove, show_hidden)

% Wrapper for the built-in function dir. Reads files from within a directory
% removes hidden files by default (that begin with dot).  
% 
% You can filter for files with a matching regular expression (regexp_to_select)
% or remove files with a matching regular expression (regexp_to_remove). You can
% optionally show the hidden files as well (show_hidden = true).
% 
% Last modified 2016-07-16 by Sam NH

if nargin < 2
    regexp_to_select = [];
end

if nargin < 3
    regexp_to_remove = [];
end

if nargin < 4
    show_hidden = false;
end

% read in file names
filenames = dir(inputdir);

% cell array of file names
filecell = cell(length(filenames),1);
for fileindex = 1:length(filenames)
    filecell{fileindex} = filenames(fileindex).name;
end

% remove hidden files
if ~show_hidden
    outputfiles = filecell(setxor(strmatch('.',filecell), (1:length(filecell))));
else
    % discard just . and .. files
    outputfiles = outputfiles(3:end);
end

% select files with desired regular expression
if ~isempty(regexp_to_select)
    xi = ~cellfun(@isempty, regexp(outputfiles, regexp_to_select));
    outputfiles = outputfiles(xi);
end

% select files with desired regular expression
if ~isempty(regexp_to_remove)
    xi = ~cellfun(@isempty, regexp(outputfiles, regexp_to_remove));
    outputfiles = outputfiles(~xi);
end




