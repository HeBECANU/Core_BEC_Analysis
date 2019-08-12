function set_up_project_path(path_to_project_root,folders_to_add,cd_proj_root)
%this function is used to set up the path and current directory
%for a function in a project

% inputs
% path_to_project_root  - path to the root folder of the project
% folders_to_add        - cell array of char vectors, eg {'dev','lib','bin'}
% cd_proj_root          - should the current directory be changed to the project root

if nargin<1
    path_to_project_root='.';
end

if nargin<2
    folders_to_add={'dev','lib','bin','test'};
end

if nargin<3
    cd_proj_root=true;
end

% if the path to root is relative then se the current directory to the calling function
if strcmp(path_to_project_root(1),'.')
    % find the path of the script/function that called this function 
    fun_stack_obj=dbstack;

    calling_function_name=fun_stack_obj(2).file;
    %fprintf('name of file what path will be relative to "%s" \n',calling_function_name)
    calling_function_path = which(calling_function_name);
    calling_function_path=fileparts(calling_function_path); %remove the function name
    path_to_project_root=fullfile(calling_function_path,path_to_project_root);
end

%change to the project root if requested and if the current folder is not that already
if cd_proj_root && ~strcmp(pwd,cd_proj_root)
    cd(path_to_project_root)
end

possible_genpath_places={'./lib/Core_BEC_Analysis/bin/genpath_exclude/',...
    './bin/genpath_exclude/'};
%try to find genpath_exclude without adding everything
look_for_genpath=1;
found_genpath=0;
search_index=1;
while look_for_genpath
    if search_index>numel(possible_genpath_places)
        look_for_genpath=false;
    elseif exist(possible_genpath_places{search_index},'dir')==7
        look_for_genpath=false;
        found_genpath=true;
        path_to_genpath=possible_genpath_places{search_index};
    end
end

%if it has not been found in these locations try the brute option
if ~found_genpath
    %add eveything in all the subfolders so that we can find genpath_exclude
    % Add that folder plus all subfolders to the path.
    addpath(genpath(path_to_project_root));%add all subfolders to the path to find genpath_exclude
    if exist('genpath_exclude','file')~=2
        error('could not find genpath_exclude check that it exists in your project folder')
    end
    path_to_genpath=fileparts(which('genpath_exclude'));
end

path(pathdef) %clean up the path back to the default state to remove all the .git that were added
addpath(path_to_project_root)
addpath(path_to_genpath)
addpath(genpath_exclude(fullfile(path_to_project_root,'lib'),'\.')) %dont add hidden folders
addpath(genpath_exclude(fullfile(path_to_project_root,'dev'),'\.'))
addpath(genpath_exclude(fullfile(path_to_project_root,'bin'),'\.'))




end