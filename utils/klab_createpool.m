function klab_createpool(CPUs)

try
    myCluster = gcp('nocreate');
    if isempty(myCluster)
        myCluster = parcluster('local');
        if ~isempty(CPUs)
            myCluster.NumWorkers=CPUs;
        end
        parpool(myCluster);
    else
        if myCluster.NumWorkers>CPUs;
            delete(myCluster);
            myCluster = parcluster('local');
            if ~isempty(CPUs)
                myCluster.NumWorkers=CPUs;
            end
            parpool(myCluster);
        end
    end
catch err
    warning('FAILED TO CREATE POOL: %s',err.message);
end