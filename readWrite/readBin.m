function data = readBin(regopsfile,frame_start,framecount,plane,isred)
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

HOME = pwd;
[path,regopsfile] = fileparts(regopsfile);
if isempty(path)
   path = pwd; 
end

data = [];

fclose('all');

try
    cd(path);
    fprintf('--- loading RegOps file...')
    load(regopsfile);
    fprintf(' done\n');
    
    if nargin<4 && length(ops1)>1
        warning('Plane not defined for %i-plane data, getting data from first plane!',length(ops1));
    end
    ops = ops1{plane};
    
    [Ly, Lx] = size(ops.mimg);
    ntotframes = ceil(sum(ops.Nframes));
    
    if frame_start+framecount>ntotframes
        warning('Requested %i frames, only %i available from startpoint %i!',framecount-frame_start+1,ntotframes-frame_start,frame_start);
        framecount = ntotframes-frame_start;
    end
    
    data = zeros(Ly,Lx,framecount,'single');
    
    if isred
        ops.RegFile=[ops.RegFile(1:end-4),'_RED.bin'];
    end
    [~,binfile,id] = fileparts(ops.RegFile);
    
    fid = fopen([binfile,id]);
    try
        fprintf('--- reading frames %i to %i from plane %i...',frame_start,frame_start+framecount-1,plane);
        fseek(fid,Ly*Lx*frame_start*2,'bof');
        data = single(fread(fid,Ly*Lx*framecount, '*uint16'));
        data = reshape(data, Ly, Lx, []);
        fprintf('success!\n');
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

