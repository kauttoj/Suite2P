clc;
clearvars;
close all;

SUITE2P_PATH = 'C:\Users\skOpti_6\Documents\MATLAB\klab_Suite2P';

cfg.do_nonridig_correction = 0;  % set to 1 if you have large FOV and/or lots of fast motion
cfg.channels = 2;           % 1 (only green) or 2 (green+red). If less than real number, remaining channels omitted.
cfg.planes = 6;             % how many planes
cfg.image_FOV = [17,778,3,512];  % [left,right,up,down]
%cfg.grating_size = [0,0]; 	% [left,right] in pixels AFTER cropping operation with selected FOV (interpolation to remove empty lines)
cfg.grating_size = [76,12];
cfg.cell_diameter = 10; 	% in pixels! (14 is optimal for WAVE stuff)
%cfg.framerate = 15.4883;   % acquisition! Automatically converted for planes. Keep empty to read from SBX files (recommended).
%cfg.tempfile_folder = [];   % path to tempfile, leave empty to use local path (recommended)
cfg.tempdata_folder = []; 	
cfg.delete_raw_tiffs = 1; 	% if 1, all raw tiffs are kept (not recommended!)
cfg.write_aligned_tiffs = 0; % if 1, all final tiffs are saved (.bin is kept always)
cfg.write_diagnostic_videos = 0; % make videos of aligned data
cfg.use_GPU = 0; %          % only for systems with good nvidia GPU (can be unstable)
cfg.fix_baseline = 0; 		% 1 = normalize intensity, 2 = normalize intensity + detrending (for each file). If no water loss occurs, set to 0!
cfg.num_cores = 1;          % set maximum number of cores to use (need to limit this if memory runs out)
cfg.use_cluster_computing = 0;   % must have this 0 when not using cluster!
cfg.write_red_bin = 1; % write aligned red channel
% If you notice spurious motion correction spikes (say >30 pixels), use these options
cfg.max_shift_limit = 30;      % max total shift in PIXELS, all above this are omitted. Does not currently work with non-ridig registration (disabled by authors).
%cfg.use_phase_correlation = 1; % default is 1, but you can try 0 for classic (standard) correlation

%%% INPUT FILES AND OUTPUT FOLDER

cfg.experiment_ID = 'T060_1L_10252016';
cfg.sbxfiles{1}='C:\Users\skOpti_6\Documents\MATLAB\klab_Suite2P\TESTDATA\T060_1L_10252016_000_026.sbx';
cfg.sbxfiles{2}='C:\Users\skOpti_6\Documents\MATLAB\klab_Suite2P\TESTDATA\T060_1L_10252016_000_027.sbx';

cfg.outpath = 'C:\Users\skOpti_6\Documents\MATLAB\klab_Suite2P\TESTDATA\output';

%%% run pipeline
cfg.programpath = SUITE2P_PATH;
addpath(SUITE2P_PATH);
klab_run_suite2p(cfg);
