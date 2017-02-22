function cfg = klab_suite2pConverter(cfg)
% Converter SBX to TIFF

if ischar(cfg)
    load(cfg);
end

if ~exist(cfg.outpath,'dir')
    try
        mkdir(cfg.outpath);
    catch me
        error('Failed to create output folder! (%s)',me.message);
    end
end

if cfg.channels>2
    error('This converter does not currently support more than 2 channels (1=Green and 2=Red) !!')
end

for i=1:length(cfg.sbxfiles)
    [~,b]=fileparts(cfg.sbxfiles{i});
    if isempty(strfind(b,cfg.experiment_ID))
        error('!!! Experiment ID %s does not match with filename %s !!!',cfg.experiment_ID,cfg.sbxfiles{i});
    end
    if ~exist([cfg.sbxfiles{i}],'file')
        error('SBX file %s not found !!!',cfg.sbxfiles{i});
    end
    if ~exist([cfg.sbxfiles{i}(1:(end-4)),'.mat'],'file')
        error('.mat file for %s not found !!!',cfg.sbxfiles{i});
    end
end

cfg.processing_stage = 0;

if cfg.use_cluster_computing
        
    cfg = doClusterConversion(cfg);
    
else
    
    cfg = klab_sbx2tiff(cfg);
    
end

cfg.expred(isnan(cfg.expred))=[];
cfg.expts(isnan(cfg.expts))=[];
if isempty(cfg.expred)
    cfg.expred=[];
end

if nnz(isnan(cfg.expts))>0
    error('expts vector is bad!')
end
if nnz(isnan(cfg.expred))>0
    error('expred vector is bad!')
end

% if not set (overridden), set framerate here
if ~isfield(cfg,'framerate') || (isfield(cfg,'framerate') && isempty(cfg.framerate))
    framerate = [];
    for i=1:length(cfg.fileinfo)
        framerate(i)=cfg.fileinfo{i}.info.framerate;
    end
    
    if any(abs((framerate-median(framerate))/median(framerate)))>0.10
        for i=1:length(framerate)
            fprintf('file %i: %f\n',framerate(i)) ;
        end
        error('!!!!! Framerates do not match between files (over 10%% deviations found) !!!!!');
    end
    
    cfg.framerate = mean(framerate); % use same mean framerate for all files
end

fprintf('\nGreen channel subfolders are: ');
for i=1:length(cfg.expts)
   fprintf('%i ',cfg.expts(i));
end
fprintf('\n');
fprintf('Red channel subfolders are: ');
for i=1:length(cfg.expred)
   fprintf('%i ',cfg.expred(i));
end
fprintf('\n\n');

cfg.processing_stage = 1;

fprintf('Conversion finished!\n\n');

end

function cfg = doClusterConversion(cfg)

    JOBNAME_PREFIX = 'job_output_conversion_file';
    jobfiles_stage0 = [];
    jobnames = [];    
    
    CONFIGFILE = cfg.CONFIGFILE;
    
    JOBDIR = cfg.JOBDIR;
    
    N_files = length(cfg.sbxfiles);
    
    for i = 1:N_files
        CONFIGFILES{i} = [JOBDIR,filesep,cfg.experiment_ID,'_CONFIGFILE_file',num2str(i)];
        save(CONFIGFILES{i},'cfg');
    end
    
    for i = 1:N_files
        
        CFG = CONFIGFILES{i};
        
        filename = fullfile(JOBDIR,[JOBNAME_PREFIX,num2str(i),'.pbs']);
        
        dlmwrite(filename, '#!/bin/bash', '');
        dlmwrite(filename, '#PBS -S /bin/bash','-append','delimiter','');
        dlmwrite(filename, '#PBS -j oe','-append','delimiter','');
        dlmwrite(filename, ['#PBS -o ',JOBDIR,filesep,JOBNAME_PREFIX,num2str(i),'.out'],'-append','delimiter','');
        dlmwrite(filename, '#PBS -m n','-append','delimiter','');
        dlmwrite(filename, '#PBS -l nodes=1:ppn=1,mem=20gb','-append','delimiter','');
        dlmwrite(filename, '#PBS -l walltime=0:40:00','-append','delimiter','');
        dlmwrite(filename, '#PBS -q default','-append','delimiter','');
        dlmwrite(filename, 'hostname','-append','delimiter','');
        dlmwrite(filename, 'echo "job starting"','-append','delimiter','');
        dlmwrite(filename, 'module load matlab-8.6','-append','delimiter','');
        dlmwrite(filename, ['cd ',cfg.programpath],'-append','delimiter','');
        dlmwrite(filename, 'echo "RUNNING MATLAB"','-append','delimiter','');
        dlmwrite(filename, sprintf('matlab -nosplash -nodisplay -nodesktop -r "add_suite2p_paths;klab_sbx2tiff(''%s'',%i,''%s'');exit;"',CONFIGFILE,i,CFG),'-append','delimiter','');
        dlmwrite(filename, 'module unload matlab-8.6','-append','delimiter','');
        dlmwrite(filename, 'echo "job finished"','-append','delimiter','');
        [notused,jobnames{i}] = system(['qsub ' filename]);
        
        outfiles_stage0{i} = [JOBDIR,filesep,JOBNAME_PREFIX,num2str(i),'.out'];
        jobfiles_stage0{i} = filename;
        
    end
    
    fprintf('\nWaiting for all jobs to finish (preprocessing stage 1)\n')
    wait_for_jobs(CONFIGFILES,1,outfiles_stage0,jobfiles_stage0,jobnames);
        
    expred = nan(1,N_files);
    expts = nan(1,N_files);
    
%     for i = 1:N_files
%         klab_sbx2tiff(CONFIGFILE,i,CONFIGFILES{i});
%     end
    
    for i = 1:N_files
                
        A=load(CONFIGFILES{i});
        
        ind = find(~isnan(A.cfg.expts));
        expts(ind) = A.cfg.expts(ind);
        
        ind = find(~isnan(A.cfg.expred));
        expred(ind) = A.cfg.expred(ind);
                
        cfg.fileinfo{i}.info = A.cfg.fileinfo{i}.info;
        cfg.fileinfo{i}.sbxinfo = A.cfg.fileinfo{i}.sbxinfo;
        cfg.fileinfo{i}.folders = A.cfg.fileinfo{i}.folders;
        
    end   
    
    cfg.expred=expred;
    cfg.expts=expts;

end




