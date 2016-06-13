function random_seed = ResetRandStream2(random_seed)
% Resets the random stream with a given seed.
%
% -- Example --
% ResetRandStream2(1); 
% rand
% ResetRandStream2(2); 
% rand
% ResetRandStream2(1);
% rand
% 
% Also allows arbitrary inputs, via a hashing function:
% 
% ResetRandStream2({'Version 1', 3}); 
% rand
% ResetRandStream2({'Version 2', 3}); 
% rand
% ResetRandStream2({'Version 1', 3}); 
% rand
% 
% Updated 15-07-21 by Sam NH to allow arbitrary inputs.

% Use data hashing to produce a random seed
if length(random_seed) > 1 || any(~isnumeric(random_seed))
    opt.Format = 'double';
    hashvalues = DataHash(random_seed,opt);
    random_seed = sum(hashvalues(1:3).*256.^(0:2));
end

% Reinitializing random seed
random_stream = RandStream('mt19937ar','Seed', random_seed);
try 
    RandStream.setDefaultStream(random_stream); %#ok<SETRS>
catch
    RandStream.setGlobalStream(random_stream);
end
