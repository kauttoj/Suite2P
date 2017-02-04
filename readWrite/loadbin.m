function data = loadbin(opsfile,binfile)

[~,a] = fileparts(binfile);
plane = str2double(a(end));
try
    .Nframes;
    load(opsfile);

[Ly, Lx] = size(ops1{plane}.mimg);
ntotframes = ceil(sum(ops1{plane}.Nframes));

fid = fopen(binfile);

try
    fprintf('reading %i frames from plane %i...',ntotframes,plane);
    data = fread(fid,  Ly*Lx*ntotframes, '*uint16');
    data = reshape(data, Ly, Lx, []);
    fprintf('success!');
catch err
    fprintf(' FAILED (%s)!',err.message);
end
fclose(fid);

if min(data(:))<0
   warning('Data has total %i negative values (possible overflow)!',nnz(data<0));
end

data = single(data);

end

