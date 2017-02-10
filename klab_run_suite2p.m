function cfg = klab_run_suite2p(cfg)

if exist(cfg.programpath, 'dir')
	addpath(genpath(cfg.programpath)) % add local path to the toolbox
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

diary off
pause(0.5);
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
    
klab_createpool(cfg.num_cores);

if exist([cfg.outpath,filesep,cfg.experiment_ID,'_suite2p_CONFIGFILE.mat'],'file')
    
    fprintf('!! Old CFG file found, trying to use that...\n');
    old_cfg = load([cfg.outpath,filesep,cfg.experiment_ID,'_suite2p_CONFIGFILE.mat']);
    
    try 
        for i=1:length(old_cfg.cfg.sbxfiles)
            if ~strcmp(old_cfg.cfg.sbxfiles{i},cfg.sbxfiles{i})
               assert(0);
            end
        end
        cfg.mouse_name=old_cfg.cfg.mouse_name;
        cfg.date=old_cfg.cfg.date;
        cfg.expts=old_cfg.cfg.expts;
        cfg.expred=old_cfg.cfg.expred;
        cfg.fileinfo=old_cfg.cfg.fileinfo;
        fprintf(' success!! Skipping SBX conversion...\n');
    catch
       fprintf(' FAILED. Running full SBX conversion again.\n');
       cfg = klab_suite2pConverter(cfg); 
    end
else
    cfg = klab_suite2pConverter(cfg);
end

save([cfg.outpath,filesep,cfg.experiment_ID,'_suite2p_CONFIGFILE.mat'],'cfg');

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
ops0.temp_tiff              = [tempdir,filesep,'suite2p_tempfile.tiff']; % copies each remote tiff locally first, into this file
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
ops0.PhaseCorrelation       = 1; % set to 0 for non-whitened cross-correlation
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
ops0.saveNeuropil           = 1;

% spike deconvolution options
ops0.imageRate              = cfg.framerate;   % imaging rate (cumulative over planes!). Approximate, for initialization of deconvolution kernel.
ops0.sensorTau              = 0.4; % decay half-life (or timescale). Approximate, for initialization of deconvolution kernel.
ops0.maxNeurop              = Inf; % for the neuropil contamination to be less than this (sometimes good, i.e. for interneurons)
ops0.recomputeKernel        = 1; % whether to re-estimate kernel during optimization (default kernel is "reasonable", if you give good timescales)
ops0.sameKernel             = 1; % whether the same kernel should be estimated for all neurons (robust, only set to 0 if SNR is high and recordings are long)
ops0.detrend_window = 100;   % in secods, 50s = 1/50 = 0.02Hz

% red channel options
% redratio = red pixels inside / red pixels outside
% redcell = redratio > mean(redratio) + redthres*std(redratio)
% notred = redratio < mean(redratio) + redmax*std(redratio)
ops0.redthres               = 1.5; % the higher the thres the less red cells
ops0.redmax                 = 1; % the higher the max the more NON-red cells

%% RUN THE PIPELINE HERE
db0 = db;

run_pipeline(db(1), ops0);

try
    diary(diaryfile);
catch err
    warning('Failed to write to diary: %s',err.message);
end
    
% deconvolved data into (dat.)cl.dcell, and neuropil subtraction coef
if cfg.do_deconvolution
    add_deconvolution(ops0, db0(1));
end

% add red channel information (if it exists)
if isfield(db0,'expred') && ~isempty(db0(1).expred)
    
    fprintf('\n-------- adding red channel (%s) ----------\n\n',char(datetime('now')));
    
    ops0.nchannels_red = db0(1).nchannels_red;
    
    %run_REDaddon(1, db0, ops0) ; % create redcell array
    
    ops1 = build_ops3(db(1), ops0);
    
    % custom function to u
    mimgR = klab_red_channel_mean(ops1);
    add_red_channel(ops1, mimgR, []);
    
    iexp=1;
    DetectRedCells; % fills dat.cl.redcell and dat.cl.notred
end

save([cfg.outpath,filesep,cfg.experiment_ID,'_suite2p_CONFIGFILE.mat'],'cfg');

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
    fprintf(' done (total % TIFF files removed)\n',count);    
    try
        rmdir([cfg.tempdata_folder,filesep,cfg.mouse_name,filesep,cfg.date]);
        rmdir([cfg.tempdata_folder,filesep,cfg.mouse_name]);
        rmdir([cfg.tempdata_folder]);
    catch
end

try
    rmdir(root);
    rmdir([ops0.ResultsSavePath,filesep,db(1).mouse_name]);
catch
    
end

fprintf('\n-------- All done! (%s, took %imin) ----------\n\n',char(datetime('now')),round(toc(starttime)/60));

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

