clc;
clearvars;
close all;

SUITE2P_PATH = 'C:\Users\skOpti_6\Documents\MATLAB\klab_Suite2P';

cfg.do_nonridig_correction = 0; % set to 1 for large FOV and/or with lots of fast motion (grooming). Note: Will disable cfg.max_shift_limit and bad frame detection
cfg.channels = 2; 				% 1 (only green) or 2 (green+red). If less than real number, remaining channels omitted.
cfg.planes = 1; 				% how many planes
cfg.image_FOV = [17,778,3,512]; 	% [left,right,up,down]
cfg.grating_size = [0,0]; 			% [left,right] in pixels AFTER cropping operation with selected FOV (interpolation to remove empty lines)
cfg.cell_diameter = 14; 			% in pixels! (14 is optimal for WAVE stuff)
%cfg.framerate = 15.4883;   		% acquisition! Automatically converted for planes. Keep empty to read from SBX files (recommended).
%cfg.tempfile_folder = []; 			% default is subfolder in output
cfg.delete_raw_tiffs = 1; 			% if 0, all raw tiffs are kept
cfg.write_aligned_tiffs = 0; 		% if 1, all final tiffs are saved (.bin is kept anyway)
cfg.write_diagnostic_videos = 0; 	% make video of aligned data
cfg.use_GPU = 1; %				% for systems with good nvidia GPU (can fail for some hardware/version combinations)
cfg.fix_baseline = 0; 			% 1 = normalize intensity, 2 = normalize intensity + detrending (for each file). If no water loss occurs, can use 0
cfg.num_cores = 2;				% set maximum number of cores to use (reduce if memory runs out)
cfg.use_cluster_computing = 0; % use 1 only with CNBC cluster
cfg.write_red_bin = 1; % write aligned red channel

% If you notice spurious motion correction spikes (say >30 pixels), use these options
cfg.max_shift_limit = 30;      % max total shift in PIXELS, all above this are omitted. Does not currently work with non-ridig registration (disabled by authors).
%cfg.use_phase_correlation = 1; % default is 1, but you can try 0 for classic (standard) correlation

%%% INPUT FILES AND OUTPUT FOLDER

cfg.experiment_ID = 'T111_1R_170407';
cfg.sbxfiles{1}='C:\Users\skOpti_6\Documents\MATLAB\klab_Suite2P\TESTDATA\T111_1R_170407_000_014.sbx';
cfg.sbxfiles{2}='C:\Users\skOpti_6\Documents\MATLAB\klab_Suite2P\TESTDATA\T111_1R_170407_000_015.sbx';

cfg.outpath = 'C:\Users\skOpti_6\Documents\MATLAB\klab_Suite2P\TESTDATA\output_single';

%%% run pipeline
cfg.programpath = SUITE2P_PATH;
addpath(SUITE2P_PATH);
klab_run_suite2p(cfg);
