clc;
clearvars;
close all;

SUITE2P_PATH = 'Z:\commonKlabData\klab_Suite2P_2.14.2017';

cfg.do_nonridig_correction = 0;
cfg.channels = 1; 				% 1 (only green) or 2 (green+red). If less than real number, remaining channels omitted.
cfg.planes = 1; 				% how many planes
cfg.do_deconvolution = 1; 			% should be always 1
cfg.image_FOV = [17,778,3,512]; 	% [left,right,up,down]
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
cfg.use_cluster_computing = 0; % use 1 only with CNBC cluster

cfg.programpath = SUITE2P_PATH;
addpath(SUITE2P_PATH);

DATAFOLDER{1} = '/mnt/Sparrowhawk/RawCaDataArchive/Pati/2036_2R_12102016/';
DATAFOLDER{2} = '/mnt/Sparrowhawk/RawCaDataArchive/Pati/2068_2R_12192016/';
DATAFOLDER{3} = '/mnt/Sparrowhawk/RawCaDataArchive/Pati/2068_2L_161222/';
DATAFOLDER{4} = '/mnt/Sparrowhawk/RawCaDataArchive/Pati/2069_2R_161229/';
DATAFOLDER{5} = '/mnt/Sparrowhawk/RawCaDataArchive/Pati/2077_1R_170105/';
DATAFOLDER{6} = '/mnt/Sparrowhawk/RawCaDataArchive/Pati/2077_1R_170112/';
DATAFOLDER{7} = '/mnt/Sparrowhawk/RawCaDataArchive/Pati/2068_2R_170124/';

sbxfiles = [];
experiment_IDs = [];
outpaths = [];

for dummy = 1:length(DATAFOLDER)
    k=0;
    files = dir([DATAFOLDER{dummy},'*.sbx']);
    for i=1:length(files)
        if files(i).isdir==0
            [~,a,b] = fileparts(files(i).name);
            if strcmp(b,'.sbx')
                k=k+1;
                sbxfiles{dummy}{k} = [DATAFOLDER{dummy},a,b];
                if ~exist(sbxfiles{dummy}{k},'file')
                    error('File ''%s'' does not found!',sbxfiles{dummy}{k});
                end
            end
        end
    end
    experiment_IDs{dummy} = a(1:end-8);
    
    % all old files will be overwritten in outpath
    outpaths{dummy} = ['/mnt/CurrentProcessing/WAVE_512_suite2p/',experiment_IDs{dummy}];
    
end

for dummy = 1:length(DATAFOLDER)
    cfg.sbxfiles = sbxfiles{dummy};
    cfg.outpath = outpaths{dummy};
    cfg.experiment_ID = experiment_IDs{dummy};
    klab_run_suite2p(cfg);
end
