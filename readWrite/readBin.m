function data = readBin(regopsfile,frame_start,framecount,plane)
% Function to read suite2P binary (uint16) into 3d matrix
%
% regopsfile = file generated after registering, e.g., regops_2080_1R_170131_suite2p.mat. Binary must be in the same folder!
% frame_start = starting frame starting from index 0 (i.e., first frame = 0)
% framecount = how many frame to return starting from frame_start
% plane = which plane to read (default = 1)


if nargin<4
    plane = 1;
end

HOME = pwd;
[path,regopsfile] = fileparts(regopsfile);
if isempty(path)
   path = pwd; 
end
data = [];

try
    cd(path);
    fprintf('--- loading RegOps file...')
    load(regopsfile);
    fprintf(' done\n');
    
    if nargin<4 && length(ops1)>1
        warning('Plane not defined for %i-plane data, gettin data from first plane!',length(ops1));
    end
    ops = ops1{plane};
    
    [Ly, Lx] = size(ops.mimg);
    ntotframes = ceil(sum(ops.Nframes));
    
    if frame_start+framecount>ntotframes
        warning('Requested %i frames, only %i available from startpoint %i!',framecount,ntotframes-frame_start,frame_start);
        framecount = ntotframes-frame_start;
    end
    
    data = zeros(Ly,Lx,framecount,'single');
    
    [~,binfile,id] = fileparts(ops.RegFile);
    
    fid = fopen([binfile,id]);
    try
        fprintf('--- reading frames %i to %i from plane %i...',frame_start,frame_start+framecount-1,plane);
        fseek(fid,Ly*Lx*frame_start,'bof');
        data = fread(fid,Ly*Lx*framecount, '*uint16');
        data = reshape(data, Ly, Lx, []);
        fprintf('success!\n');
    catch err1
        fprintf(' FAILED (%s)!',err1.message);
    end
    fclose(fid);
catch err2    
    warning('!!!!!!! Failed to read data, reason: %s',err2.message);
    if exist('fid','var')
        fclose(fid);
    end
end

cd(HOME);

end

