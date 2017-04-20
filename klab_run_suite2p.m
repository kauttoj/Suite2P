function cfg = klab_run_suite2p(cfg)

if exist(cfg.programpath, 'dir')
	addpath(genpath(cfg.programpath)); % add local path to the toolbox
else
	error('toolbox_path does not exist, please change toolbox_path');
end

diaryfile = [cfg.outpath,filesep,cfg.experiment_ID,'_suite2p_diary.txt'];

if ~exist(cfg.outpath,'dir')
    try
        mkdir(cfg.outpath);
    catch me
        error('Failed to create output folder! (%s)',me.message);
    end
end

cfg.tempdata_folder = [cfg.outpath,filesep,'raw_data'];

diary off
pause(2); % allow some time to write buffer
if exist(diaryfile,'file')
    try
        delete(diaryfile);
    catch
        warning('old diary exists, but it could not be removed, appending to old');
    end
end
try
    diary(diaryfile);
catch err
    warning('Failed to write to diary: %s',err.message);
end

% run converter to create cropped TIFF files and suite2p structure from SBX data

starttime = tic;

fprintf('\n-------- Starting Suite2P pipeline (%s) ----------\n\n',char(datetime('now')));

if ~isfield(cfg,'num_cores')
    cfg.num_cores=[];
end

if ~isfield(cfg,'tempfile_folder') || ( isfield(cfg,'tempfile_folder') && isempty(cfg.tempfolder) )
    cfg.tempfile_folder = [cfg.tempdata_folder,filesep];
end

if ~isfield(cfg,'write_red_bin')
   cfg.write_red_bin = 0;
end 

cfg.mouse_name = cfg.experiment_ID; % subfolder, level1
cfg.date = 'suite2p'; % subfolder, level 2
cfg.processing_stage = 0;

if ~isfield(cfg,'use_cluster_computing')
    cfg.use_cluster_computing = 1;
end

if cfg.use_cluster_computing
    fprintf('\nNote: Pipeline is running in CLUSTER mode using PBS scheduler (not for desktop use!)\n\n');
    
    JOBDIR = [cfg.outpath,filesep,'jobfiles'];
    if ~exist(JOBDIR,'dir')
        mkdir(JOBDIR);
    else
        for i=1:20
            JOBDIR = [cfg.outpath,filesep,'jobfiles',num2str(i)];
            if ~exist(JOBDIR,'dir')
                mkdir(JOBDIR);
                break;
            end
        end
        if i==20
            error('you already have 20 job folders!')
        end
    end
    
    cfg.JOBDIR = JOBDIR;
    
else
    klab_createpool(cfg.num_cores);
end

cfg.CONFIGFILE = [cfg.outpath,filesep,cfg.experiment_ID,'_suite2p_CONFIGFILE.mat'];

if exist(cfg.CONFIGFILE,'file')
    
    fprintf('!! Old CFG file found, trying to use that...');
    old_cfg = load([cfg.outpath,filesep,cfg.experiment_ID,'_suite2p_CONFIGFILE.mat']);
    
    try 
        assert(length(old_cfg.cfg.sbxfiles)<=length(cfg.sbxfiles));
        for i=1:length(old_cfg.cfg.sbxfiles)
            if ~strcmp(old_cfg.cfg.sbxfiles{i},cfg.sbxfiles{i})
               assert(0);
            end
        end
        for i=1:length(cfg.sbxfiles)
            if ~exist(old_cfg.cfg.fileinfo{i}.folders{1},'dir')
                assert(0);
            end
        end
        cfg.mouse_name=old_cfg.cfg.mouse_name;
        cfg.date=old_cfg.cfg.date;
        cfg.expts=old_cfg.cfg.expts;
        cfg.expred=old_cfg.cfg.expred;
        cfg.framerate = old_cfg.cfg.framerate;
        cfg.fileinfo=old_cfg.cfg.fileinfo;
        cfg.tempdata_folder=old_cfg.cfg.tempdata_folder; 
        fprintf(' success!! Skipping SBX conversion\n\n');        
    catch err
       fprintf(' FAILED (reason: %s)! Running full SBX conversion\n',err.message);
       save(cfg.CONFIGFILE,'cfg');
       cfg = klab_suite2pConverter(cfg); 
    end
else
    save(cfg.CONFIGFILE,'cfg');
    cfg = klab_suite2pConverter(cfg);
end

cfg.processing_stage = 1;
save(cfg.CONFIGFILE,'cfg');

ops0.klab_info = cfg;

%%------------------------------------
db(1).mouse_name    = cfg.mouse_name;
db(1).date          = cfg.date;
db(1).expts         = cfg.expts;
db(1).expred        = cfg.expred;
db(1).nchannels = 1;
db(1).gchannel = 1;
db(1).nplanes = cfg.planes;
db(1).diameter = cfg.cell_diameter;
db(1).nonrigid = cfg.do_nonridig_correction;
db(1).nchannels_red = cfg.channels - 1;

ops0.toolbox_path = cfg.programpath;

% mex -largeArrayDims SpikeDetection/deconvL0.c (or .cpp) % MAKE SURE YOU COMPILE THIS FIRST FOR DECONVOLUTION

ops0.useGPU                 = cfg.use_GPU; % if you can use an Nvidia GPU in matlab this accelerates registration approx 3 times. You only need the Nvidia drivers installed (not CUDA).
ops0.fig                    = 0; % turn off figure generation with 0

% root paths for files and temporary storage (ideally an SSD drive. my SSD is C:/)
ops0.RootStorage            = cfg.tempdata_folder; % Suite2P assumes a folder structure, check out README file
ops0.temp_tiff              = [cfg.tempfile_folder,'suite2p_tempfile.tiff']; % copies each remote tiff locally first, into this file
ops0.RegFileRoot            = cfg.outpath;  % location for binary file
ops0.DeleteBin              = 0; % set to 1 for batch processing on a limited hard drive
ops0.ResultsSavePath        = cfg.outpath; % a folder structure is created inside
if cfg.write_aligned_tiffs
    ops0.RegFileTiffLocation = cfg.outpath;
else
    ops0.RegFileTiffLocation    = []; %'D:/DATA/'; % leave empty to NOT save registered tiffs (slow)
end
ops0.fix_baseline 		 = cfg.fix_baseline;

% registration options
ops0.doRegistration         = 1; % skip (0) if data is already registered
ops0.showTargetRegistration = 0; % shows the image targets for all planes to be registered
if ~isfield(cfg,'use_phase_correlation')
    ops0.PhaseCorrelation       = 1; % set to 0 for non-whitened cross-correlation
else
    ops0.PhaseCorrelation       = cfg.use_phase_correlation;
end
ops0.SubPixel               = Inf; % 2 is alignment by 0.5 pixel, Inf is the exact number from phase correlation
ops0.NimgFirstRegistration  = 500; % number of images to include in the first registration pass 
ops0.nimgbegend             = 250; % frames to average at beginning and end of blocks

ops0.dobidi                 = 0; % is bidirectional scanning?
ops0.showTargetRegistration = 0;

% cell detection options
ops0.ShowCellMap            = 0; % during optimization, show a figure of the clusters
ops0.sig                    = 0.5;  % spatial smoothing length in pixels; encourages localized clusters
ops0.nSVDforROI             = 1300; % was 1000, how many SVD components for cell clustering
ops0.NavgFramesSVD          = 6000; % was 5000, how many (binned) timepoints to do the SVD based on
ops0.saveNeuropil           = 0;

% spike deconvolution options
ops0.imageRate              = cfg.framerate;   % imaging rate (cumulative over planes!). Approximate, for initialization of deconvolution kernel.
ops0.sensorTau              = 0.4; % decay half-life (or timescale). Approximate, for initialization of deconvolution kernel.
ops0.maxNeurop              = Inf; % for the neuropil contamination to be less than this (sometimes good, i.e. for interneurons)
ops0.recomputeKernel        = 1; % whether to re-estimate kernel during optimization (default kernel is "reasonable", if you give good timescales)
ops0.sameKernel             = 1; % whether the same kernel should be estimated for all neurons (robust, only set to 0 if SNR is high and recordings are long)
ops0.detrend_window         = 100;   % in secods, 100s <--> 1/100 = 0.01Hz

cfg.do_deconvolution = 1; % always need run this to get refined neuropil coefficients!

if ~isfield(cfg,'max_shift_limit')
    ops0.max_shift_limit = inf;
else
    ops0.max_shift_limit = cfg.max_shift_limit; 
    fprintf('\nNote: Setting maximum motion correction shift radius to %i pixels\n\n',ops0.max_shift_limit);
end


    
% red channel options
% redratio = red pixels inside / red pixels outside
% redcell = redratio > mean(redratio) + redthres*std(redratio)
% notred = redratio < mean(redratio) + redmax*std(redratio)
ops0.redthres               = 1.5; % the higher the thres the less red cells
ops0.redmax                 = 1; % the higher the max the more NON-red cells

%% RUN THE PIPELINE HERE
db0 = db;

if cfg.use_cluster_computing
    run_pipeline_cluster(db(1), ops0,cfg);
else    
    run_pipeline(db(1), ops0);
    add_deconvolution(ops0, db0(1));
end

try
    diary(diaryfile);
catch err
    warning('Failed to write to diary: %s',err.message);
end

% add red channel information (if it exists)
if isfield(db0,'expred') && ~isempty(db0(1).expred)
    
    fprintf('\n-------- adding red channel (%s) ----------\n\n',char(datetime('now')));
    
    ops0.nchannels_red = db0(1).nchannels_red;
    ops0.write_red_bin = cfg.write_red_bin;
            
    if cfg.use_cluster_computing
        
        DetectRedCells_cluster(db(1),ops0,cfg); % fills dat.cl.redcell and dat.cl.notred
        
    else        
                      
        DetectRedCells(db(1),ops0); % fills dat.cl.redcell and dat.cl.notred
        
    end
end

save(cfg.CONFIGFILE,'cfg');

try
    if cfg.use_cluster_computing
        create_diagnostic_figures_cluster(ops0, db0(1),cfg);
    else
        create_diagnostic_figures(ops0, db0(1));
    end
catch err
    warning('Failed to create diagnostic figures!, reason: %s',err.message);
end


if exist(ops0.temp_tiff,'file')
    delete(ops0.temp_tiff);
end

root = [ops0.ResultsSavePath,filesep,db(1).mouse_name,filesep,db(1).date];
d = dir(root);
for i=1:length(d)
    if d(i).isdir && ~(strcmp(d(i).name,'.') || strcmp(d(i).name,'..'))
        outpath = [root,filesep,d(i).name];
    else
       continue; 
    end
    dd = dir(outpath);
    for j=1:length(dd) 
        if ~(strcmp(dd(j).name,'.') || strcmp(dd(j).name,'..'))
            movefile([outpath,filesep,dd(j).name],ops0.ResultsSavePath);
        end
    end
    rmdir(outpath);
end

if cfg.delete_raw_tiffs==1    
    fprintf('\nRemoving RAW tiff files...');
    count = 0;
    for i=1:length(cfg.fileinfo)
        for j=1:length(cfg.fileinfo{i}.folders)
            d = dir(cfg.fileinfo{i}.folders{j});
            for k=1:length(d)
                if ~d(k).isdir && ~(strcmp(d(k).name,'.') || strcmp(d(k).name,'..'))
                    [~,~,type]=fileparts(d(k).name);
                    if strcmp(type,'.tiff')
                        delete([cfg.fileinfo{i}.folders{j},filesep,d(k).name]);
                        count = count + 1;
                    end
                end
            end
            try
                rmdir(cfg.fileinfo{i}.folders{j});
            catch
                
            end
        end
    end    
    fprintf(' done (total %i TIFF files removed)\n',count);    
    try
        rmdir([cfg.tempdata_folder,filesep,cfg.mouse_name,filesep,cfg.date]);
        rmdir([cfg.tempdata_folder,filesep,cfg.mouse_name]);
        rmdir([cfg.tempdata_folder]);
    catch err
        warning('Failed to clean folders: %s',err.message);
    end
end

try
    rmdir(root);
    rmdir([ops0.ResultsSavePath,filesep,db(1).mouse_name]);
catch
    
end

try
    if cfg.write_diagnostic_videos
        if cfg.use_cluster_computing
            create_diagnostic_movies_cluster(ops0, db0(1),cfg);
        else
            create_diagnostic_movies(ops0, db0(1));
        end
    else
       fprintf('\nSkipping diagnostic videos\n');
    end
catch err
    warning('Failed to create diagnostic movies!, reason: %s',err.message);
end

fprintf('\n\n-------- All done! (%s, took %imin) ----------\n\n',char(datetime('now')),round(toc(starttime)/60));

try
    diary(diaryfile);
catch err
    warning('Failed to write to diary: %s',err.message);
end

diary off;

%% STRUCTURE OF RESULTS FILE
% cell traces are in dat.Fcell
% neuropil traces are in dat.FcellNeu
% manual, GUI overwritten "iscell" labels are in dat.cl.iscell

% stat(icell) contains all other information
% autoamted iscell label, based on anatomy
% neuropil subtraction coefficient 
% st are the deconvolved spike times (in frames)
% c  are the deconvolved amplitudes
% kernel is the estimated kernel

%     stat(icell).c                      = dcell{icell}.c;     % spike amplitudes
%     stat(icell).st                     = dcell{icell}.st; % spike times    
%     stat(icell).neuropilCoefficient    = coefNeu(icell); 
%     stat(icell).noiseLevel             = sn(icell);     % noise level
%     stat(icell).kernel                 = dcell{icell}.kernel;     % noise level

% f = trace
% s = spike
% b = baseline
% p = neuropil
% beta = neuropil contamination coefficient
% 
% model:
% f = conv(s,kernel) + beta*p + b
% i.e., 

end

%% ALL STUFF BELOW ARE ONLY FOR CLUSTER COMPUTING (writing and submitting jobs)

function create_diagnostic_movies_cluster(ops0,db,cfg)

db_FILE = [cfg.JOBDIR,filesep,'db_file.mat'];
ops_FILE = [cfg.JOBDIR,filesep,'ops_file.mat'];

save(db_FILE,'db');
save(ops_FILE,'ops0');

JOBNAME_PREFIX = 'job_output_diagmovies_file';

JOBDIR = cfg.JOBDIR;

fold_root = [cfg.outpath,filesep,'diagnostics'];
fold = [fold_root,filesep,'videos'];
if ~exist(fold,'dir')
    mkdir(fold);
end

%%----
N_files = length(cfg.sbxfiles);

for i = 1:N_files
    CONFIGFILES{i} = [JOBDIR,filesep,cfg.experiment_ID,'_CONFIGFILE_file',num2str(i),'.mat'];
    cfg.files_to_process = i;
    save(CONFIGFILES{i},'cfg');
end


for i = 1:N_files
        
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
    dlmwrite(filename, sprintf('matlab -nosplash -nodisplay -nodesktop -r "add_suite2p_paths;create_diagnostic_movies(''%s'',''%s'',''%s'');exit;"',ops_FILE,db_FILE,CONFIGFILES{i}),'-append','delimiter','');
    dlmwrite(filename, 'module unload matlab-8.6','-append','delimiter','');
    dlmwrite(filename, 'echo "job finished"','-append','delimiter','');
    [notused,jobname{i}] = system(['qsub ' filename]);
    
    outfiles_stage1{i} = [JOBDIR,filesep,JOBNAME_PREFIX,num2str(i),'.out'];
    jobfiles_stage1{i} = filename;
    
end

fprintf('\nWaiting for job to finish (diagnostic movies)\n')
wait_for_jobs(CONFIGFILES,6,outfiles_stage1,jobfiles_stage1,jobname);

%create_diagnostic_movies(ops_FILE,db_FILE,CONFIGFILES);

end

function create_diagnostic_figures_cluster(ops0,db,cfg)

db_FILE = [cfg.JOBDIR,filesep,'db_file.mat'];
ops_FILE = [cfg.JOBDIR,filesep,'ops_file.mat'];

save(db_FILE,'db');
save(ops_FILE,'ops0');

CONFIGFILES = cfg.CONFIGFILE;

JOBDIR = cfg.JOBDIR;

%%----

JOBNAME_PREFIX = 'job_output_diagfigures_file';

filename = fullfile(JOBDIR,[JOBNAME_PREFIX,'.pbs']);

dlmwrite(filename, '#!/bin/bash', '');
dlmwrite(filename, '#PBS -S /bin/bash','-append','delimiter','');
dlmwrite(filename, '#PBS -j oe','-append','delimiter','');
dlmwrite(filename, ['#PBS -o ',JOBDIR,filesep,JOBNAME_PREFIX,'.out'],'-append','delimiter','');
dlmwrite(filename, '#PBS -m n','-append','delimiter','');
dlmwrite(filename, '#PBS -l nodes=1:ppn=1,mem=16gb','-append','delimiter','');
dlmwrite(filename, '#PBS -l walltime=0:40:00','-append','delimiter','');
dlmwrite(filename, '#PBS -q default','-append','delimiter','');
dlmwrite(filename, 'hostname','-append','delimiter','');
dlmwrite(filename, 'echo "job starting"','-append','delimiter','');
dlmwrite(filename, 'module load matlab-8.6','-append','delimiter','');
dlmwrite(filename, ['cd ',cfg.programpath],'-append','delimiter','');
dlmwrite(filename, 'echo "RUNNING MATLAB"','-append','delimiter','');
dlmwrite(filename, sprintf('matlab -nosplash -nodisplay -nodesktop -r "add_suite2p_paths;create_diagnostic_figures(''%s'',''%s'',''%s'');exit;"',ops_FILE,db_FILE,CONFIGFILES),'-append','delimiter','');
dlmwrite(filename, 'module unload matlab-8.6','-append','delimiter','');
dlmwrite(filename, 'echo "job finished"','-append','delimiter','');
[notused,jobname{1}] = system(['qsub ' filename]);

outfiles_stage1{1} = [JOBDIR,filesep,JOBNAME_PREFIX,'.out'];
jobfiles_stage1{1} = filename;

fprintf('\nWaiting for job to finish (diagnostic figures)\n')
wait_for_jobs(CONFIGFILES,5,outfiles_stage1,jobfiles_stage1,jobname);

%create_diagnostic_figures(ops_FILE,db_FILE,CONFIGFILES);

end


function DetectRedCells_cluster(db, ops0,cfg)

db_FILE = [cfg.JOBDIR,filesep,'db_file.mat'];
ops_FILE = [cfg.JOBDIR,filesep,'ops_file.mat'];

save(db_FILE,'db');
save(ops_FILE,'ops0');

CONFIGFILES = cfg.CONFIGFILE;

JOBDIR = cfg.JOBDIR;

%%----

JOBNAME_PREFIX = 'job_output_redcell_file';

filename = fullfile(JOBDIR,[JOBNAME_PREFIX,'.pbs']);

dlmwrite(filename, '#!/bin/bash', '');
dlmwrite(filename, '#PBS -S /bin/bash','-append','delimiter','');
dlmwrite(filename, '#PBS -j oe','-append','delimiter','');
dlmwrite(filename, ['#PBS -o ',JOBDIR,filesep,JOBNAME_PREFIX,'.out'],'-append','delimiter','');
dlmwrite(filename, '#PBS -m n','-append','delimiter','');
dlmwrite(filename, '#PBS -l nodes=1:ppn=1,mem=20gb','-append','delimiter','');
dlmwrite(filename, '#PBS -l walltime=0:40:00','-append','delimiter','');
dlmwrite(filename, '#PBS -q default','-append','delimiter','');
dlmwrite(filename, 'hostname','-append','delimiter','');
dlmwrite(filename, 'echo "job starting"','-append','delimiter','');
dlmwrite(filename, 'module load matlab-8.6','-append','delimiter','');
dlmwrite(filename, ['cd ',cfg.programpath],'-append','delimiter','');
dlmwrite(filename, 'echo "RUNNING MATLAB"','-append','delimiter','');
dlmwrite(filename, sprintf('matlab -nosplash -nodisplay -nodesktop -r "add_suite2p_paths;DetectRedCells(''%s'',''%s'',''%s'');exit;"',db_FILE,ops_FILE,CONFIGFILES),'-append','delimiter','');
dlmwrite(filename, 'module unload matlab-8.6','-append','delimiter','');
dlmwrite(filename, 'echo "job finished"','-append','delimiter','');
[notused,jobname{1}] = system(['qsub ' filename]);

outfiles_stage1{1} = [JOBDIR,filesep,JOBNAME_PREFIX,'.out'];
jobfiles_stage1{1} = filename;

fprintf('\nWaiting for job to finish (red cells)\n')
wait_for_jobs(CONFIGFILES,4,outfiles_stage1,jobfiles_stage1,jobname);

%DetectRedCells(db_FILE,ops_FILE,CONFIGFILES);

end


function run_pipeline_cluster(db, ops0,cfg)

db_FILE = [cfg.JOBDIR,filesep,'db_file.mat'];
ops_FILE = [cfg.JOBDIR,filesep,'ops_file.mat'];

save(db_FILE,'db');
save(ops_FILE,'ops0');

CONFIGFILES = cfg.CONFIGFILE;

JOBDIR = cfg.JOBDIR;

%%----

JOBNAME_PREFIX = 'job_output_pipeline_file';

filename = fullfile(JOBDIR,[JOBNAME_PREFIX,'.pbs']);

dlmwrite(filename, '#!/bin/bash', '');
dlmwrite(filename, '#PBS -S /bin/bash','-append','delimiter','');
dlmwrite(filename, '#PBS -j oe','-append','delimiter','');
dlmwrite(filename, ['#PBS -o ',JOBDIR,filesep,JOBNAME_PREFIX,'.out'],'-append','delimiter','');
dlmwrite(filename, '#PBS -m n','-append','delimiter','');
dlmwrite(filename, '#PBS -l nodes=1:ppn=1,mem=26gb','-append','delimiter','');
dlmwrite(filename, '#PBS -l walltime=5:00:00','-append','delimiter','');
dlmwrite(filename, '#PBS -q default','-append','delimiter','');
dlmwrite(filename, 'hostname','-append','delimiter','');
dlmwrite(filename, 'echo "job starting"','-append','delimiter','');
dlmwrite(filename, 'module load matlab-8.6','-append','delimiter','');
dlmwrite(filename, ['cd ',cfg.programpath],'-append','delimiter','');
dlmwrite(filename, 'echo "RUNNING MATLAB"','-append','delimiter','');
dlmwrite(filename, sprintf('matlab -nosplash -nodisplay -nodesktop -r "add_suite2p_paths;run_pipeline(''%s'',''%s'',''%s'');exit;"',db_FILE,ops_FILE,CONFIGFILES),'-append','delimiter','');
dlmwrite(filename, 'module unload matlab-8.6','-append','delimiter','');
dlmwrite(filename, 'echo "job finished"','-append','delimiter','');
[notused,jobname{1}] = system(['qsub ' filename]);

outfiles_stage1{1} = [JOBDIR,filesep,JOBNAME_PREFIX,'.out'];
jobfiles_stage1{1} = filename;

fprintf('\nWaiting for job to finish (main pipeline)\n')
wait_for_jobs(CONFIGFILES,2,outfiles_stage1,jobfiles_stage1,jobname);

%run_pipeline(db_FILE,ops_FILE,CONFIGFILES);

%%------

JOBNAME_PREFIX = 'job_output_deconv_file';

filename = fullfile(JOBDIR,[JOBNAME_PREFIX,'.pbs']);

dlmwrite(filename, '#!/bin/bash', '');
dlmwrite(filename, '#PBS -S /bin/bash','-append','delimiter','');
dlmwrite(filename, '#PBS -j oe','-append','delimiter','');
dlmwrite(filename, ['#PBS -o ',JOBDIR,filesep,JOBNAME_PREFIX,'.out'],'-append','delimiter','');
dlmwrite(filename, '#PBS -m n','-append','delimiter','');
dlmwrite(filename, '#PBS -l nodes=1:ppn=1,mem=20gb','-append','delimiter','');
dlmwrite(filename, '#PBS -l walltime=1:00:00','-append','delimiter','');
dlmwrite(filename, '#PBS -q default','-append','delimiter','');
dlmwrite(filename, 'hostname','-append','delimiter','');
dlmwrite(filename, 'echo "job starting"','-append','delimiter','');
dlmwrite(filename, 'module load matlab-8.6','-append','delimiter','');
dlmwrite(filename, ['cd ',cfg.programpath],'-append','delimiter','');
dlmwrite(filename, 'echo "RUNNING MATLAB"','-append','delimiter','');
dlmwrite(filename, sprintf('matlab -nosplash -nodisplay -nodesktop -r "add_suite2p_paths;add_deconvolution(''%s'',''%s'',''%s'');exit;"',ops_FILE,db_FILE,CONFIGFILES),'-append','delimiter','');
dlmwrite(filename, 'module unload matlab-8.6','-append','delimiter','');
dlmwrite(filename, 'echo "job finished"','-append','delimiter','');
[notused,jobname{1}] = system(['qsub ' filename]);

outfiles_stage2{1} = [JOBDIR,filesep,JOBNAME_PREFIX,'.out'];
jobfiles_stage2{1} = filename;

fprintf('\nWaiting for job to finish (deconvolution)\n')
wait_for_jobs(CONFIGFILES,3,outfiles_stage2,jobfiles_stage2,jobname);

%add_deconvolution(ops_FILE,db_FILE,CONFIGFILES);

end