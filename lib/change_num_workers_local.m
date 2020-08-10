function change_num_workers_local(new_num)

myCluster = parcluster('local');
myCluster.NumWorkers = new_num;  % 'Modified' property now TRUE
saveProfile(myCluster);    % 'local' profile now updated,
                           % 'Modified' property now FALSE
                           
end