function cfg = klab_sbx2tiff(cfg,INDEX,CONFIGFILE)

N_max_frames = 2000; % split TIFF's (single TIFF should never reach 4GB limit!)

HOME = pwd;

%klab_createpool( cfg.processes );

fprintf('\n-------- Klab SBX --> TIFF (suite2p format) converter module ----------\n');

if ischar(cfg)
    load(cfg);
end
   
N_files = length(cfg.sbxfiles);
if nargin>1   
   INDICES = INDEX;
else
   INDICES = 1:N_files;
end

ROOT_FOLDER = [cfg.tempdata_folder,filesep,cfg.mouse_name,filesep,cfg.date];
if ~exist(ROOT_FOLDER,'dir')
    try
        mkdir(ROOT_FOLDER);
    catch me
        error('Failed to create temp TIFF folder! (%s)',me.message);
    end
end

infiles = cfg.sbxfiles;

expts_all = cell(1,length(cfg.sbxfiles));
expred_all = cell(1,length(cfg.sbxfiles));
fileinfo = cell(1,length(cfg.sbxfiles));

for ii = INDICES
    % make sure we are in data folder (required by sbx)
    filename = infiles{ii};
    
    fprintf('Processing file %s\n',filename);
    
    [a,b,~]=fileparts(filename);
    if ~isempty(a)
        cd(a);
        filename=b;
    end
    
    fname = [filename,'_raw'];
    
    info = [];
    [data,info,orig_FOV] = readsbx(filename,cfg.image_FOV);
    fprintf(' done\n');
    
    rem = mod(size(data,4),cfg.planes);
    if rem>0
        data(:,:,:,size(data,4)+1-rem)=[];
        warning('Frame count is not multiply of planes, dropping last %i frames!',rem);
    end
    
    if nnz(cfg.grating_size)>0
        fprintf('Interpolating interleaved gratings with %i (left) and %i (right)...',cfg.grating_size(1),cfg.grating_size(2));
        tic
        data = round(do_grating_interpolation(data,cfg.grating_size,mod(cfg.image_FOV(3),2)));
        fprintf(' done (%is)\n',round(toc));
    end
    
    sbxinfo = info;
    info=[];
        
    if ndims(data)==4
        info.channels = size(data,1);
    else
        error('invalid data size (should be 4D)');
    end
    
    if cfg.channels < info.channels
        warning('!!! More channels available than specified in CFG, dropping all extra channels !!!')
        data((cfg.channels+1):end,:,:,:)=[];
    elseif cfg.channels > info.channels
        warning('!!! More channels requested than available in data !!!')
    end
    
    info.channels = size(data,1);
    info.framerate = (sbxinfo.resfreq/sbxinfo.config.lines*(2-sbxinfo.scanmode));
    
    % test if the last GREEN frame contains good data, if not, drop it
    last = squeeze(data(1,:,:,end));
    last = last(:);
    if var(double(last))<1e-5 || nnz(isnan(last))>0
        data(:,:,:,end)=[];
        rem = mod(size(data,4),cfg.planes);
        if rem>0
            data(:,:,:,size(data,4)+1-rem)=[];
            warning('Frame count is not multiply of planes, dropping last %i frames!',rem);
        end
    end
    
    bad = data>floor(intmax('uint16')*0.99);
    if nnz(bad)>0
       warning('!!! Total %i datapoints (time and space) were over 99%% of uint16 maximum, data might be clipped !!! ',nnz(bad));
    end
    bad=[];
    
    data = uint16(data);
    
    N_frames = size(data,4);
    
    info.frames = N_frames;
    
    info.orig_FOV = orig_FOV;
    
    info.planes = cfg.planes;
    
    info.frames_per_plane = N_frames/cfg.planes;
    
    if info.frames ~= info.frames_per_plane*info.planes
        error('Number of frames does not match (total %i frames, %i planes, %i frames per plane)',info.frames,info.planes,info.frames_per_plane);
    end
    
    sz = size(data);
    info.dims = [sz(2),sz(3)];
    
    sbxinfo.dims = info.dims;
    
    info.fname = filename;
    
    file_subfolder=[];
    
    all_channels(ii) = info.channels;
    
    save_parts = round(linspace(1,N_frames+1,1+ceil(N_frames/N_max_frames)));
    
    info.tiff_parts = length(save_parts)-1;
    
    expts_all{ii}  = nan(1,info.channels*N_files);
    expred_all{ii} = nan(1,info.channels*N_files);
    
    for channel = 1:info.channels
        
        offset = N_files*(channel-1);
        fold_id = ii + offset;
        
        file_subfolder{channel} = [ROOT_FOLDER,filesep,num2str(fold_id)];
        
        expts_all{ii}(ii) = ii;
        
        if channel > 1
            expred_all{ii}(offset+ii) = offset+ii;
        end
        
        if ~exist(file_subfolder{channel} ,'dir')
            try
                mkdir(file_subfolder{channel});
            catch me
                error('Failed to create temp TIFF folder! (%s)',me.message);
            end
        end
        
        for i=1:(length(save_parts)-1)
            fname_part = sprintf('%s_PART%i',fname,i);
            inds = save_parts(i):(save_parts(i+1)-1);
            if ~exist([file_subfolder{channel},filesep,fname_part,'.tiff'],'file')
                fprintf('saving TIFF data (channel %i, part %i)...',channel,i);
                savedata([file_subfolder{channel},filesep,fname_part,'.tiff'],squeeze(data(channel,:,:,inds)));
                fprintf(' done\n');
            else
                fprintf('TIFF data (channel %i, part %i) already exists! Skipping save\n',channel,i);
            end
        end
        
    end
    
    fileinfo{ii}.info = info;
    fileinfo{ii}.sbxinfo = sbxinfo;
    fileinfo{ii}.folders = file_subfolder;
    
end

expred = nan(1,cfg.channels*N_files);
expts = nan(1,cfg.channels*N_files);

for ii = INDICES
    
    ind = find(~isnan(expts_all{ii}));
    expts(ind) = expts_all{ii}(ind);
    ind = find(~isnan(expred_all{ii}));
    expred(ind) = expred_all{ii}(ind);
    
end

cfg.fileinfo = fileinfo;

%expts(isnan(expts))=[];
%expred(isnan(expred))=[];

cfg.expts = expts;
cfg.expred = expred;

cfg.processing_stage = 1;

if nargin==3
    save(CONFIGFILE,'cfg');
end

cd(HOME);

end

function savedata(fname,data)

%save(fname,'data','info','sbxinfo','-v7.3');
TiffWriter(uint16(data),fname,16);

end

function [data,info,orig_FOV] = readsbx(filename,FOV)

% get info
clearvars -global
global info;
img = sbxread(filename,0,1);

fprintf('reading data (RAW)...');

%data(:,cfg.image_FOV(3):cfg.image_FOV(4),cfg.image_FOV(1):cfg.image_FOV(2),:)
sy = length(FOV(3):FOV(4));
sx = length(FOV(1):FOV(2));

if size(img,2)>sy || size(img,3)>sx
    fprintf(' original FOV [y,x]=[%i,%i] cropped to [%i,%i] ',size(img,2),size(img,3),sy,sx);
end

for i=0:(info.max_idx-1)
    a = sbxread(filename,i,1);
    %a = sbxread(filename,0,i);
    if i==0
        data = zeros(size(a,1),sy,sx,info.max_idx,'uint16');
        orig_FOV = [size(a,2),size(a,3)];
    end
    data(:,:,:,i+1) = a(:,FOV(3):FOV(4),FOV(1):FOV(2));
end

end

function data = do_grating_interpolation(data,extend,startind)

% startind = {0,1};

sz = size(data);

bad_left = (startind+1):2:sz(2);
if bad_left(1)==1
    bad_right = (startind+1):2:sz(2);
else
    bad_right = (startind+2):2:sz(2);
end

cols_left = 1:extend(1);
cols_right = (sz(3)-extend(2)+1):sz(3);

N_bad_left = length(bad_left);
N_bad_right = length(bad_right);

data = single(data);

for i = 1:size(data,4)
    frame = (data(:,:,:,i));
    
    for r = 1:N_bad_left
        
        row = bad_left(r);
                
        if ~isempty(cols_left)
            if row == 1
                frame(:,row,cols_left) = frame(:,row+1,cols_left);
            elseif row == sz(2)
                frame(:,row,cols_left) = frame(:,row-1,cols_left);
            else
                frame(:,row,cols_left) = (frame(:,row-1,cols_left) + frame(:,row+1,cols_left))/2;
            end
        end
        
    end
    
    for r = 1:N_bad_right
        
        row = bad_right(r);
        
        if ~isempty(cols_right)
            if row == 1
                frame(:,row,cols_right) = frame(:,row+1,cols_right);
            elseif row == sz(2)
                frame(:,row,cols_right) = frame(:,row-1,cols_right);
            else
                frame(:,row,cols_right) = (frame(:,row-1,cols_right) + frame(:,row+1,cols_right))/2;
            end
        end
                
    end
    
    data(:,:,:,i) = frame;
end

end
