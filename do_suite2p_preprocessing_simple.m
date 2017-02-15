clc;
clearvars;
close all;

SUITE2P_PATH = 'Z:\commonKlabData\klab_Suite2P_2.14.2017';

cfg.do_nonridig_correction = 0;
cfg.channels = 1; 				% 1 (only green) or 2 (green+red). If less than real number, remaining channels omitted.
cfg.planes = 1; 				% how many planes
cfg.do_deconvolution = 1; 			% should be always 1
cfg.image_FOV = [40,720,8,512]; 	% [left,right,up,down]
cfg.grating_size = [0,0]; 			% [left,right] in pixels AFTER cropping operation with selected FOV (interpolation to remove empty lines)
cfg.cell_diameter = 14; 			% in pixels! (14 is optimal for WAVE stuff)
%cfg.framerate = 15.4883;   		% acquisition! Automatically converted for planes. Keep empty to read from SBX files (recommended).
cfg.tempdata_folder = []; 			% if 0, using system tempdir()
cfg.delete_raw_tiffs = 1; 			% if 1, all raw tiffs are kept
cfg.write_aligned_tiffs = 0; 		% if 1, all final tiffs are saved (.bin is kept anyway)
cfg.write_diagnostic_videos = 0; 	% make video of aligned data (NOT YET WORKING!)
cfg.use_GPU = 0; %				% for systems with good nvidia GPU (sometimes unstable)
cfg.fix_baseline = 0; 			% 1 = normalize intensity, 2 = normalize intensity + detrending (for each file). If no water loss occurs, can use 0
%cfg.num_cores = 4;				% set maximum number of cores to use (reduce if memory runs out)

cfg.programpath = SUITE2P_PATH;

cfg.experiment_ID = '2068_2R_170124';
cfg.sbxfiles{1}='z:/RawCaDataArchive/Pati/2068_2R_170124/2068_2R_170124_000_018.sbx';
cfg.sbxfiles{2}='z:/RawCaDataArchive/Pati/2068_2R_170124/2068_2R_170124_000_020.sbx';
cfg.sbxfiles{3}='z:/RawCaDataArchive/Pati/2068_2R_170124/2068_2R_170124_000_026.sbx';

cfg.outpath = 'Z:\commonKlabData\klab_Suite2P_testdata\output';

% run pipeline
addpath(SUITE2P_PATH);
klab_run_suite2p(cfg);
