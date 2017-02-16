function add_clean_signals(ops, db)
ops = build_ops3(db, ops);

fprintf('\nAdding cleaned signals into output\n');
for i = 1:length(ops.planesToProcess)
    iplane  = ops.planesToProcess(i);
    
    fpath = sprintf('%s/F_%s_%s_plane%d_proc.mat', ops.ResultsSavePath, ...
        ops.mouse_name, ops.date, iplane);
    if exist(fpath, 'file')
        A = load(fpath);
    else
        fpath = sprintf('%s/F_%s_%s_plane%d.mat', ops.ResultsSavePath, ...
            ops.mouse_name, ops.date, iplane);
        A = load(fpath);
    end
    
    if isfield(A, 'dat')
        A = A.dat; % just in case...
    end    
    
    if min([A.stat.neuropilCoefficient])<0
       warning('Total %i neuropil coefficients were less than zero, all set to zero!',sum([A.stat.neuropilCoefficient]<0));
    end
    coefNeu = max(0,[A.stat.neuropilCoefficient]');
    
    fprintf('..plane %i: %i segments with neuropil coefficients mean %3.2f with 95%% interval %3.2f - %3.2f\n',iplane,length(coefNeu),mean(coefNeu),prctile(coefNeu,2.5),prctile(coefNeu,97.5));
    for k = 1:length(A.Fcell)
        F_clean{k} = A.Fcell{k} -  bsxfun(@times,A.FcellNeu{k},coefNeu);
        F_clean{k} = bsxfun(@times,F_clean{k},A.scalefactors);
        F_clean{k} = bsxfun(@plus,F_clean{k},A.baselines);
    end
    
    save(fpath,'F_clean','-append');
    
end
%
