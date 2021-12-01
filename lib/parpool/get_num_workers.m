function num_workers=get_num_workers()
    myCluster = parcluster('local');
    num_workers=myCluster.NumWorkers;

end