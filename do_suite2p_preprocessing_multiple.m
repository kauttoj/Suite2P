clc;
clearvars;
close all;

SUITE2P_PATH = '/home/kuhlmanlab/code/klab_Suite2P';

cfg.do_nonridig_correction = 0;
cfg.channels = 1;
cfg.planes = 1;
cfg.do_deconvolution = 1;
cfg.image_FOV = [25,758,4,511]; % [left,right,up,down]
cfg.grating_size = [0,0]; % [left,right] in pixels AFTER cropping operation with selected FOV (interpolation to remove empty lines)
cfg.cell_diameter = 15; % in pixels!
cfg.framerate = 15.4883;   % For each pixel. Approximate, for initialization of deconvolution kernel.
cfg.tempdata_folder = [];%;tempdir;
cfg.delete_raw_tiffs = 1;
cfg.write_aligned_tiffs = 0;
cfg.write_diagnostic_videos = 0;
cfg.use_GPU = 0; %
cfg.fix_baseline = 1; % do piecewise linear detrending of aligned data and fix baselines between files (slow, can take ~2 hours)

addpath(SUITE2P_PATH);
cfg.programpath = SUITE2P_PATH;

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
