function create_diagnostic_movies(ops0, db,cfg)

if ischar(db) && ischar(ops0)
    load(db);
    load(ops0);
end

BATCH_SIZE = 500; % how many frames to read at once

ops = build_ops3(db, ops0);

fregops =  sprintf('regops_%s_%s.mat', ops.mouse_name, ops.date);
ops1 = load(fullfile(ops.RegFileRoot, fregops));
ops1 = ops1.ops1;

fold_root = [ops.klab_info.outpath,filesep,'diagnostics'];
if ~exist(fold_root,'dir')
    mkdir(fold_root);
end

fold = [fold_root,filesep,'videos'];
if ~exist(fold,'dir')
    mkdir(fold);
end

tic;

fprintf('\n\n--- Creating diagnostic movies (total %i)\n\n',numel(ops1));

parfor file = 1:numel(ops1)
    
    fprintf('..making video for plane %i\n',file);
    
    ops = ops1{file};
       
    ops.iplane  = file;
    [Ly, Lx] = size(ops.mimg);
    
    fid = fopen(ops.RegFile, 'r');
    
    fprintf('... %i blocks (files) and %i frames\n',length(ops.Nframes),sum(ops.Nframes));
    
    outfile = [fold,filesep,sprintf('aligned_data_plane%i_%s',file,ops.klab_info.experiment_ID)];
       
    limits = get_color_limits(fid,sum(ops.Nframes),300,Ly*Lx); % estimate intensity limits from 300 frames
    
    fprintf('... setting intensity limits to %i - %i\n',limits(1),limits(2));
    
    map = linspace(0,1,2^12);
        
    try
        writerObj = VideoWriter([outfile,'.m4v'],'MPEG-4');
        writerObj.Quality = 95;
    catch err
        writerObj = VideoWriter([outfile,'.avi'],'Motion JPEG AVI');
    end
    
    writerObj.FrameRate = 30;
    open(writerObj);
    
    frewind(fid);
                
    inds = cell(1,3);
    f = zeros(Ly,Lx,3);
    for j=1:3    
        f(:,:,j)=j;
        inds{j} = find(f==j);
    end
    ind = inds{1};
    
    while 1
        
        data = fread(fid,  Ly*Lx*BATCH_SIZE, '*uint16');
        
        if isempty(data)
           break; 
        end
        
        data = double(reshape(data,Ly,Lx,[]));
        
        for j=1:size(data,3)  
            data(:,:,j) = medfilt2(data(:,:,j), [2,2]);
        end
        
        data(data<limits(1))=limits(1);
        data(data>limits(2))=limits(2);
        data = data - min(data(:)) + 1e-4;
        data =data/(max(data(:)));
        data(data>0.50)=0.5 - 1e-4;
        data = 2*data;
        data = floor((length(map)-2)*data)+1;
        
        data(data>length(map))=length(map);
        data(data<1)=1;
        
        for j=1:size(data,3)                   
            a = data(:,:,j);
            f(inds{1}) = map(a(ind));
            f(inds{2}) = map(a(ind));
            f(inds{3}) = map(a(ind));
            frame = im2frame(f);
            writeVideo(writerObj,frame);     
        end        
        
    end
    
    close(writerObj);
    fclose(fid);    
    
end

fprintf('all done (took %is)\n\n',round(toc));

if nargin==3
    load(cfg);    
    cfg.processing_stage = 6;
    save(cfg.CONFIGFILE,'cfg');
end

end

function limits = get_color_limits(fid,frames,steps,S)

pos = round(linspace(1,frames,min(frames,steps+3)));

pos = pos(2:(end-1));
mi=0;
ma=0;

for i=1:length(pos)
   fseek(fid,S*(pos(i)-1),'bof');
   data = single(fread(fid,S,'*uint16'));
   mi = mi + prctile(data(:),2.5)/length(pos);
   ma = ma + prctile(data(:),97.5)/length(pos);
end

limits = [ceil(mi),floor(ma)];

end