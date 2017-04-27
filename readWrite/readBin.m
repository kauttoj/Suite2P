function data = readBin(regopsfile,frame_start,framecount,plane,isred,coord)
% Function to read suite2P binary (uint16) into 3d matrix
%
% regopsfile = file generated after registering, e.g., regops_2080_1R_170131_suite2p.mat. Binary must be in the same folder!
% frame_start = starting frame starting from index 0 (i.e., first frame = 0)
% framecount = how many frame to return starting from frame_start
% plane = which plane to read (default = 1)
% isred = read red (1) or green (0) binary

% Note: uint16 data type is assumed (2 bytes)!!

if nargin<4
    plane = 1;
end
if nargin<5
    isred = 0;
end
if nargin<6
    coord = [];
end

HOME = pwd;
path = pwd; 

if ischar(regopsfile)
    [path,regopsfile] = fileparts(regopsfile);
end

data = [];

fclose('all');

try
    cd(path);
    
    if ischar(regopsfile)
        fprintf('--- loading RegOps file...')
        load(regopsfile);
        fprintf(' done\n');
        if nargin<4 && length(ops1)>1
            warning('Plane not defined for %i-plane data, getting data from first plane!',length(ops1));
        end
        ops = ops1{plane};
        if isred
            ops.RegFile=[ops.RegFile(1:end-4),'_RED.bin'];
        end
        [~,binfile,id] = fileparts(ops.RegFile);
    else
        ops = regopsfile;
        binfile = ops.RegFile;
        id=[];
    end
    
    [Ly, Lx] = size(ops.mimg);
    ntotframes = ceil(sum(ops.Nframes));
    
    if ischar(framecount) && strcmp(framecount,'all')
        framecount = ntotframes;
    end
    
    if frame_start+framecount>ntotframes
        warning('Requested %i frames, only %i available from startpoint %i!',framecount-frame_start+1,ntotframes-frame_start,frame_start);
        framecount = ntotframes-frame_start;
    end        
              
    fid = fopen([binfile,id]);
    try
        fprintf('--- reading frames %i to %i from plane %i...',frame_start,frame_start+framecount-1,plane);
        fseek(fid,Ly*Lx*frame_start*2,'bof');        
        tic;
        if isempty(coord)
            data = zeros(Ly,Lx,framecount,'single');
            data = single(fread(fid,Ly*Lx*framecount, '*uint16'));
            data = reshape(data, Ly, Lx, []);
        else            
            BATCH = 1000;
            bad = coord(:,1)<1 | coord(:,2)<1 | coord(:,1)>Lx | coord(:,2)>Ly;
            coord(bad,:)=[];
            ind = sub2ind([Ly,Lx],coord(:,2),coord(:,1));
            remaining = framecount;
            position = 1;
            data = nan(1,remaining);
            PARTS = ceil(remaining/BATCH);
            PART = 1;
            fprintf(' using mask with %i pixels and %i parts\n',length(ind),PARTS);
            while remaining>0                
                if mod(PART,10)==0
                    fprintf('.. part %i of %i\n',PART,PARTS);
                end
                f = min(BATCH,remaining);
                temp_data = single(fread(fid,Ly*Lx*f, '*uint16'));
                temp_data = reshape(temp_data,Ly,Lx,f);
                temp_data = reshape(temp_data,Ly*Lx,f);
                data(position + (1:f) - 1) = mean(temp_data(ind,:));
                remaining = max(0,remaining - BATCH);
                position = position+BATCH;
                PART = PART + 1;
            end             
        end
        fprintf('success! (took %is)\n',round(toc));
    catch err1
        fprintf(' FAILED (%s)!',err1.message);
    end
    fclose(fid);
catch err2    
    warning('\n!!!!!!! Failed to read data, reason: %s\n',err2.message);
    if exist('fid','var')
        fclose(fid);
    end
end

cd(HOME);

end

