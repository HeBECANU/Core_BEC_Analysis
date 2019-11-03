cli_header(0,'TESTING stash_opts and get_stashed_opts:');
% pass_labels = {'PASSED','FAILED'};
opts.stash_dir = '.\stash_test';
opts.label = 'this_test';
opts.verbose=0;
% fprintf('Here! %s',pass_labels{1})
% Write one instance
stash_opts(opts);
% Write another instance after changing a field
opts.newfield = 'abcdefg';
stash_opts(opts);
retrieved=get_stashed_opts(opts.stash_dir,opts.label);
test1 = isequal(opts,retrieved.opts);
cli_header(2,'Write, update, and retrieve test %s',pass_labels{2-test1});


% Clean up the side effects...
cli_header(1,'Cleaning up after test...');
all_delete = fullfile(opts.stash_dir,'*.mat');
delete(all_delete)
rem_dir = dir(opts.stash_dir);
rem_dir = {rem_dir.name};
rem_dir = rem_dir(~ismember(rem_dir,{'.','..'}));
if isempty(rem_dir)
    cli_header(2,'Test directory cleared');
else
    cli_header(2,'Test directory NOT cleared');
end
rmdir(opts.stash_dir)
if exist(opts.stash_dir,'dir')
    cli_header(2,'Test directory NOT removed.');
else
    cli_header(1,'Test cleanup complete.');
end
