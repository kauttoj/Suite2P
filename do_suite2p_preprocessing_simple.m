clc;
clearvars;
close all;

SUITE2P_PATH = 'C:\Users\skOpti_6\Documents\MATLAB\klab_Suite2P';

cfg.do_nonridig_correction = 0;
cfg.channels = 1; 				% 1 (only green) or 2 (green+red). If less than real number, remaining channels omitted.
cfg.planes = 1; 				% how many planes
cfg.do_deconvolution = 1; 			% should be always 1
cfg.image_FOV = [17,778,3,512]; 	% [left,right,up,down]
cfg.grating_size = [0,0]; 			% [left,right] in pixels AFTER cropping operation with selected FOV (interpolation to remove empty lines)
cfg.cell_diameter = 14; 			% in pixels! (14 is optimal for WAVE stuff)
%cfg.framerate = 15.4883;   		% acquisition! Automatically converted for planes. Keep empty to read from SBX files (recommended).
cfg.delete_raw_tiffs = 1; 			% if 0, all raw tiffs are kept (not recommended)
cfg.write_aligned_tiffs = 1; 		% if 1, all final tiffs are saved (.bin is kept always)
cfg.write_diagnostic_videos = 1; 	% make video of aligned data (NOT YET WORKING!)
cfg.use_GPU = 0; %				% for systems with good nvidia GPU (sometimes unstable)
cfg.fix_baseline = 0; 			% 1 = normalize intensity, 2 = normalize intensity + detrending (for each file). If no water loss occurs, can use 0
%cfg.num_cores = 4;				% set maximum number of cores to use (reduce if memory runs out)
cfg.use_cluster_computing = 0; % use 1 only with CNBC cluster
cfg.max_shift_limit = 30;      % max total shift in PIXELS, all above this are marked as bad (default = inf, no limits) and do not limit FOV. Does not work with non-ridig registration (disabled by authors)

cfg.programpath = SUITE2P_PATH;

cfg.experiment_ID = '2093_NC_170331';
cfg.sbxfiles{1}='D:\TEMP\testdata\2093_NC_170331_000_004.sbx';
cfg.sbxfiles{2}='D:\TEMP\testdata\2093_NC_170331_000_005.sbx';

cfg.outpath = 'D:\TEMP\testdata\output';

% run pipeline
addpath(SUITE2P_PATH);
klab_run_suite2p(cfg);
