function retrieved_opts = get_stashed_opts(stash_dir,label)

% A quick function for retrieving the latest stashed 
% structures with a standard filename stash__label__timestamp.mat
% Input:
    % stash_dir         Full form of directory where structures are kept
    % label             Unique identifier for the option set
% Output:
    % retrieved_opts    The last-written option stash from stash_opts()

% Example:
% opts.stash_dir = '.\stash_test';
% opts.label = 'this_test';
% stash_opts(opts)
% retrieved_opts = get_stashed_opts(opts.stash_dir,opts.label);


%Default output in case of failure
retrieved_opts = nan;

if ~exist(stash_dir,'dir')
    error("Your stash doesn't exist!")
end
dir_contents = dir(stash_dir);
dir_names = {dir_contents.name};
dir_names=dir_names(~ismember(dir_names,{'.','..'}));
if numel(dir_names) == 0
    warning('Your stash is empty!')
else
    data_dirs = cell(numel(dir_names),1);
    delim_idx = zeros(numel(dir_names),2);
    
    ctr = 0;
    for ii = 1:numel(dir_names)
        this_fname = dir_names{ii};
        [start_idx,end_idx] = regexp(this_fname,sprintf('__%s__',label));
        if length([start_idx,end_idx])==2
            ctr = ctr + 1;
            data_dirs{ctr}=this_fname;
            delim_idx(ctr,:) = [start_idx,end_idx];
        end
    end
    if ctr == 0
        warning("Couldn't find any matches in stash")
    else
        data_dirs = data_dirs(1:ctr);
        delim_idx = delim_idx(1:ctr,:);
        timestamps = zeros(numel(data_dirs),1);
        for idx = 1:numel(data_dirs)
            timestamps(idx) = str2double(data_dirs{idx}(delim_idx(idx,2)+1:delim_idx(idx,2)+10));
        end
        [~,place] = max(timestamps);
        latest_dir = data_dirs{place};
        retrieved_opts = open(fullfile(stash_dir,latest_dir));
    end
end
end