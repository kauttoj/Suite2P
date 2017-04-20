function mimgR = klab_red_channel_mean(ops)

% numPlanes = length(ops.planesToProcess);
if ~isfield(ops,'write_red_bin')
    warning('Field write_red_bin not set, setting default 0 (no binaries will be written for red)\n');
    ops.write_red_bin = 0;
end
%% build file list with red channel

if (isfield(ops, 'SubDirsRed') && ~isempty(ops.SubDirsRed))
   subDirsRed = ops.SubDirsRed;
else
    if (isfield(ops, 'expred') && ~isempty(ops.expred))
        for i = 1:length(ops.expred)
            subDirsRed{i} = sprintf('%d', ops.expred(i));
        end
    else
        warning('could not find red channel info, returning...')
        return;
    end 
end
%%
% build file list
fsrootRED = [];

for j = 1:length(subDirsRed)
    fsrootRED{j} = dir(fullfile(ops.RootDir, subDirsRed{j}, '*.tiff'));
    for k = 1:length(fsrootRED{j})
         fsrootRED{j}(k).name = fullfile(ops.RootDir, subDirsRed{j}, fsrootRED{j}(k).name);
    end
end

%%
if (isfield(ops, 'SubDirs') && ~isempty(ops.SubDirs))
   subDirsGreen = ops.SubDirs;
else
    if (isfield(ops, 'expred') && ~isempty(ops.expts))
        for i = 1:length(ops.expts)
            subDirsGreen{i} = sprintf('%d', ops.expts(i));
        end
    else
        warning('could not find green channel info to match files, returning...')
        return;
    end 
end
% build file list
fsrootGREEN = [];

for j = 1:length(subDirsRed)
    fsrootGREEN{j} = dir(fullfile(ops.RootDir, subDirsGreen{j}, '*.tiff'));
    for k = 1:length(fsrootGREEN{j})
         fsrootGREEN{j}(k).name = fullfile(ops.RootDir, subDirsGreen{j}, fsrootGREEN{j}(k).name);
    end
end

%% check if the RED channel file has already been registered, find indices to shift by
fsRED   = [];
fs      = [];
nimg    = [];
nimgall = [];

nimgFirst = zeros(length(ops.fsroot), ops.nplanes);

for j = 1:length(ops.fsroot)
    if ~isempty(ops.fsroot{j})
        ops.fsroot{j}(1).nFrames = nFrames(ops.fsroot{j}(1).name);
        
        for k = 1:length(ops.fsroot{j})
            if abs(ops.fsroot{j}(1).bytes - ops.fsroot{j}(k).bytes)<10
                ops.fsroot{j}(k).nFrames = ops.fsroot{j}(1).nFrames;
            else
                ops.fsroot{j}(k).nFrames = nFrames(ops.fsroot{j}(k).name);
            end
        end
        
        fs = [fs {ops.fsroot{j}.name}];
        nimg{j} = [ops.fsroot{j}.nFrames];
        nimgall = [nimgall ops.fsroot{j}.nFrames];
    end
end
%
nindx = zeros(1, ops.nplanes);
for k = 1:length(fs)
    nimgFirst(k, :) = nindx;
    nindx = nindx + ceil((-[0:1:(ops.nplanes-1)] + nimgall(k)/ops.nchannels_red)/ops.nplanes);
end
%%
for j = 1:length(fsrootRED)
    fsRED = [fsRED {fsrootRED{j}.name}];
end
indx = 1:length(fsRED);

fsGREEN   = [];
for j = 1:length(fsrootGREEN)
    fsGREEN = [fsGREEN {fsrootGREEN{j}.name}];
end

for k = 1:length(fsRED) 
    [~,a]=fileparts(fsRED{k});
    [~,b]=fileparts(fsGREEN{k});   
    if isempty(strcmp(a,b))
       error('Red and Green filenames mismatch!!') 
    end  
end

% sum(nimgall(indx))/ops.nchannels/ops.nplanes
ntifs = ceil(1000/length(fsRED));

%
DS = cell(ops.nplanes, 1);
try
    root = ops.ResultsSavePath;
    fname = sprintf('regops_%s_%s.mat', ops.mouse_name, ops.date);
    load(fullfile(root, fname))
    for j = 1:ops.nplanes
        DS{j} = ops1{j}.DS;
    end
catch err
    error(' Error loading motion corrections params from regops (to align Red): %s!',err.message);
end
%
    
fprintf('\nCreating mean red image from subset of red frames\n');
ntf0 = 0;
numPlanes = ops.nplanes;
iplane0 = 1:1:ops.nplanes;
for k = 1:length(fsRED)
    %if nimgall(indx(k))>=median(nimgall(indx))        
        iplane0 = mod(iplane0-1, numPlanes) + 1;
         
        ichanset = [ops.nchannels_red*ops.nplanes + [ops.nchannels_red (ntifs*ops.nchannels_red*ops.nplanes)] ...
            ops.nchannels_red]; 
        data = loadFrames(fsRED{k}, ichanset(1),ichanset(2), ichanset(3));        
        if ~exist('mimgR')
            [Ly, Lx, ~] = size(data);
            mimgR = zeros(Ly, Lx, ops.nplanes);
        end
        data = reshape(data, Ly, Lx, ops.nplanes, ntifs);         

        for j = 1:size(data,3)
            dsall = DS{j}(nimgFirst(indx(k), j)+1 + [1:size(data,4)], :);
            data(:, :, j,:)        = ...
                    register_movie(data(:, :, j, :), ops, dsall);
        end
        mimgR = mimgR + mean(data(:,:,iplane0,:), 4);           
        ntf0 = ntf0 + 1;                
    %end
end

mimgR = mimgR/ntf0;

if ops.write_red_bin
    
    try
        tic;
        fprintf('\nPreparing aligned red channel binaries\n');
        old = 0;
        for file = 1:length(fsrootRED)
            markers = [0,cumsum(nimg{file})];
            data = zeros(Ly,Lx,markers(end),'uint16');
            for set = 1:length(fsrootRED{file})
                
                a = loadFrames(fsrootRED{file}(set).name,1,nimg{file}(set),1);
                data(:,:,(markers(set)+1):markers(set+1)) = a;
                
            end
            data = reshape(data, Ly, Lx, ops.nplanes,markers(end)/ops.nplanes);
            for j = 1:size(data,3)
                dsall = DS{j}(old+(1:size(data,4)), :);
                data(:, :, j,:)        = ...
                    register_movie(data(:, :, j, :), ops, dsall);
            end
            fprintf('...writing aligned file %i (%i planes, %i frames)\n',file,ops.nplanes,size(data,4));
            write_red_bin(ops,data,file);
            old = old + size(data,4);
        end
        fprintf('all done! (took %is)\n\n',round(toc));
        
    catch err
        
        warning('\nFailed to write red binary files!!!: %s\n',err.message);
        
    end
end

end

function write_red_bin(ops,data,file)

data = uint16(data);

numPlanes = size(data,3);

for i = 1:numPlanes
    
    % open bin file for writing
    RegFile = fullfile(ops.RegFileRoot, ...
        sprintf('%s_%s_%s_plane%d_RED.bin', ops.mouse_name, ops.date, ...
        ops.CharSubDirs,i));
    regdir = fileparts(RegFile);
    if ~exist(regdir, 'dir')
        mkdir(regdir);
    end
    
    if file>1        
        fid = fopen(RegFile, 'a');        
    else
        fid = fopen(RegFile, 'w');
    end
    
    fwrite(fid,squeeze(data(:,:,i,:)), class(data));
    
    fclose(fid);
    
end

end