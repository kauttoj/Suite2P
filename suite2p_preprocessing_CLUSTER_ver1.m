clearvars
close all
clc 

SKIP_FILE_TRANSFER = 0; % set to 1 if you want to skip initial data copying (you already have data on Cluster)


% Template for Ca data preprocessing, including data copying with rsync
% Updated 1/27/17 to new room 191 FOV settings after scanbox card was replaced.  Will also copy .txt file into /processed folder

%% Initial step to copy data from Sparrowhawk

% Raw Calcium data folder at Sparrowhawk storage, this is permanent, ALL FILES AND SUBFOLDERS ARE COPIED
RAW_DATA_PATH='/datasmart/kuhlmanlab/RawCaDataArchive/Pati/2036_2R_170223/';%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Processed data (after processing is finished) Calcium data folder at Sparrowhawk storage, this is permanent (leave empty to skip step)
PROCESSED_DATA_PATH = '/datasmart/kuhlmanlab/ProcessedDataArchive/WAVE_512/2036_2R_170223/';%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Processing folder at CNBC cluster, this is a temporary storage
CLUSTER_PATH = '/data2/kuhlmanlab/CurrentProcessingSKa/Pati/2036_2R_170223/';%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% here we create the command string
COMMAND = sprintf('rsync -av --chmod=+rw --stats --inplace -L -e "ssh -l kuhlmanlab" sparrowhawk.cnbc.cmu.edu:%s %s',RAW_DATA_PATH,CLUSTER_PATH);

% execute command
if not(SKIP_FILE_TRANSFER)
   system(COMMAND);
end



%%----------------- processing settings start here -------------------

cfg.experiment_ID = '2036_2R_170223';   % experiment ID, must match your SBX files!

cfg.do_nonridig_correction = 0;		% set to 1 if there is lots of fast motion or you have a large FOV
cfg.channels = 1; 			% 1 (only green) or 2 (green+red). If less than real number, remaining channels are omitted.
cfg.planes = 1; 			% how many planes
cfg.image_FOV = [17,778,3,512]; 	% [left,right,up,down]
cfg.grating_size = [0,0]; 		% [left,right] in pixels AFTER cropping operation with selected FOV (interpolation to remove empty lines)
cfg.cell_diameter = 14; 		% in pixels! (14 is optimal for WAVE stuff)
cfg.delete_raw_tiffs = 1; 		% if 1, all raw tiffs are kept
cfg.write_aligned_tiffs = 0; 		% if 1, all final tiffs are saved (.bin is kept anyway)
cfg.write_diagnostic_videos = 1; 	% make video of aligned data
cfg.use_GPU = 0; %			% for systems with good nvidia GPU (sometimes unstable, does not work with cluster!)
cfg.fix_baseline = 0; 			% 1 = normalize intensity, 2 = normalize intensity + detrending (for each file, slow!). Set to 0 if data is good as is.

DATAPATH = CLUSTER_PATH;
cd(DATAPATH);

% grab all sbx files
allfiles = dir('*.sbx');
for i = 1:length(allfiles)
    cfg.sbxfiles{i} = [pwd,filesep,allfiles(i).name];
end

SUITE2P_PATH = '/data2/kuhlmanlab/code/klab_Suite2P';
cfg.programpath = SUITE2P_PATH;

cfg.outpath = [DATAPATH,filesep,'processed_suite2p'];	% change this if you want to re-run processing with different settings


% run code
addpath(SUITE2P_PATH);
klab_run_suite2p(cfg);  % start processing




%%--------------------------------------------------------
%% copy processed files to Sparrowhawk

fprintf('Copying data to Sparrowhawk\n');

CLUSTER_PATH_PROCESSED = cfg.outpath;
if ~isempty(PROCESSED_DATA_PATH)
    COMMAND = sprintf('rsync -av --stats --inplace -L -e "ssh -l kuhlmanlab" %s sparrowhawk.cnbc.cmu.edu:%s',CLUSTER_PATH_PROCESSED,PROCESSED_DATA_PATH);
    system(COMMAND);

    % copy over the log .txt file
    CLUSTER_PATH_TEXTFILE = [CLUSTER_PATH cfg.experiment_ID '.txt'];
    PROCESSED_DATA_PATH_TEXTFILE = [PROCESSED_DATA_PATH '/processed'];
    COMMAND = sprintf('rsync -av --stats --inplace -L -e "ssh -l kuhlmanlab" %s sparrowhawk.cnbc.cmu.edu:%s',CLUSTER_PATH_TEXTFILE,PROCESSED_DATA_PATH_TEXTFILE);
    system(COMMAND);
end
