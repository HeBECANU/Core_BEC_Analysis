function new_num=change_num_workers_local(new_num)
% no input  sets to number of logical processors (with warning) 
% 0 sets to number of logical processors (without warning)
% -1 sets to number of logical processors minus 1 (so computer is usable when running)
% -2 sets to number of physical processors



if nargin<1 || isempty(new_num)
    new_num=str2double(getenv('NUMBER_OF_PROCESSORS'));
    warning('setting number of workers to number of logical processors depending on the application you may be better setting it to the number of physical processors')
    %new_num = feature('numcores')
end
    

if isequal(new_num,0)
    new_num=str2double(getenv('NUMBER_OF_PROCESSORS'));
elseif isequal(new_num,-1)
    new_num=str2double(getenv('NUMBER_OF_PROCESSORS'))-1;    
elseif isequal(new_num,-2)
    new_num = feature('numcores');
end


myCluster = parcluster('local');
myCluster.NumWorkers = new_num;  % 'Modified' property now TRUE
saveProfile(myCluster);    % 'local' profile now updated,
                           % 'Modified' property now FALSE
                           
end