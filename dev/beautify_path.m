function signal = beautify_path(dir_name)
% Runs MBeautifier on all .m files throughout a directory tree

allfiles = dir(dir_name);
filenames = {allfiles.name};
are_dirs = [allfiles.isdir];
is_dir = filenames(are_dirs);
is_dir = is_dir(~ismember(is_dir,{'.','..'}));
is_mfile = filenames(endsWith(filenames,'.m'));

[~,short_name]=fileparts(dir_name);
cli_header(0,'Beautifying directory %s...',short_name);

for thisfile = is_mfile
    cli_header(2,'Formatting %s...',thisfile{1});
    MBeautify.formatFile(fullfile(dir_name,thisfile{1}),fullfile(dir_name,thisfile{1}))
end


for subdir = is_dir
   beautify_path(fullfile(dir_name,subdir{1}));
   cli_header(1,'Done.');
end

cli_header(0,' All done!');
end