function flag = stash_opts(stash_dir,opts)
% A quick function for caching the configuration files used in a given data
% run for structured recollection by get_stashed_opts() later on.
% Required inputs: A structure s with mandatory fields:
    % stash_dir         The destination for storage (full form)
    % label             Identifier to be inserted as part of standard fname
% Outputs:      
    % flag              Binary, 0 = write failed, 1 = write success
    
% Example:
% opts.stash_dir = '.\stash_test';
% opts.label = 'this_test';
% stash_opts(opts)
% retrieved_opts = get_stashed_opts(opts.stash_dir,opts.label);

% Possible improvements:
% Check to see whether an identical option-set exists.
% Simple version: Check latest saved option.
% Or: Hash the input structure and search for matching hash before
    % writing
    
flag = 0;
if ~isfield(opts,'label')
   error('You forgot the stash identifier!'); 
end
if ~exist(stash_dir,'dir')
    mkdir(stash_dir);
end
timestamp = sprintf('%.0f',posixtime(datetime()));
% Test to see whether a stash already exists
    % If is it the same as the present input, write nothing
fname = sprintf('stash__%s__%s.mat',opts.label,timestamp);
savedir=fullfile(stash_dir,fname);

try
    save(savedir,'opts','-nocompression', '-v7.3')
    flag = 1;
catch
    warning('Structure write failed');
end

if isfield(opts,'verbose')
    if opts.verbose>0
        cli_header(2,'Options stashed');
    end
end


end