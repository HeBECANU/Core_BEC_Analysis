function out_status=mkdir_if_not_exist(dirpath)
%mkdir_if_not_exist - check if a directory path exists and if it does not, make as many folders as needed
% built on code initialy posted in https://gist.github.com/ferryzhou/2269380
%
% Syntax:         status=mkdir_if_not_exist(dirpath)
%
% Inputs:
%    dirpath    - path to folder
%
% Outputs:
%    status     - 1 it does exist
%                 0 it did not but now it does
%                 
% Example: 
%     out_dir='./out'
%     mkdir_if_not_exist(out_dir);
%
% Other m-files required: none
% Also See: none
% Subfunctions: none
% MAT-files required: none
%
% Known BUGS/ Possible Improvements
%  - none 
%
% Author: Bryce Henson
% email: Bryce.Henson@live.com
% Last revision:2020-10-16

%------------- BEGIN CODE --------------

if dirpath(end) ~= filesep
    dirpath = strcat(dirpath,filesep); 
end

out_status=1;
if (exist(dirpath, 'dir') == 0) 
    mkdir(dirpath)
    out_status=0;
end

end
