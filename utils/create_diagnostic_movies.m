function create_diagnostic_movies(ops0, db,cfg_file)

if ischar(db) && ischar(ops0)
    load(db);
    load(ops0);
end

CUT_PERCENTAGE = 0.75;

ops = build_ops3(db, ops0);

fregops =  sprintf('regops_%s_%s.mat', ops.mouse_name, ops.date);
ops1 = load(fullfile(ops.RegFileRoot, fregops));
ops1 = ops1.ops1;

frame_limits = [0,cumsum(ops1{1}.Nframes)]; 
frame_counts = ops1{1}.Nframes;
files_to_process = 1:length(frame_counts);
if nargin==3
   A=load(cfg_file);
   files_to_process = A.cfg.files_to_process;
   A=[];
end   

fold_root = [ops.klab_info.outpath,filesep,'diagnostics'];
if ~exist(fold_root,'dir')
    mkdir(fold_root);
end

fold = [fold_root,filesep,'videos'];
if ~exist(fold,'dir')
    mkdir(fold);
end

tic;

map = linspace(0,1,2^12);                

fprintf('\n\n--- Creating diagnostic movies (total %i)\n\n',numel(ops1));

for plane = 1:numel(ops1)
    
    fprintf('..making video for plane %i\n',plane);
    
    ops = ops1{plane};
    
    ops.iplane  = plane;
    [Ly, Lx] = size(ops.mimg);
    
    inds = cell(1,3);
    f = zeros(Ly,Lx,3);
    for j=1:3
        f(:,:,j)=j;
        inds{j} = find(f==j);
    end
    ind = inds{1};
    
    fid = fopen(ops.RegFile, 'r');
    
    for block = files_to_process
        
        [~,str]=fileparts(ops.klab_info.sbxfiles{block});
        
        fprintf('... block %i with %i frames\n',block,frame_counts(block));
        
        outfile = [fold,filesep,sprintf('aligned_data_plane%i_%s',plane,str)];
        
        fseek(fid,Ly*Lx*frame_limits(block)*2,'bof');
        limits = get_color_limits(fid,frame_limits(block)+1,frame_counts(block),200,Ly*Lx); % estimate intensity limits from 300 frames
        
        fprintf('... intensity limits are [%i,%i] (cutted to [%i,%i])\n',limits(1),limits(2),limits(1),round(limits(1)+range(limits)*CUT_PERCENTAGE));        
        
        fseek(fid,Ly*Lx*frame_limits(block)*2,'bof');
        data = fread(fid,Ly*Lx*frame_counts(block), '*uint16');
        data = single(reshape(data,Ly,Lx,[]));
        
        if isempty(data)
            warning('Data for block %i was empty!',block);
            continue;
        end
        
        try
            writerObj = VideoWriter([outfile,'.m4v'],'MPEG-4');
            writerObj.Quality = 95;
        catch err
            writerObj = VideoWriter([outfile,'.avi'],'Motion JPEG AVI');
        end
        
        writerObj.FrameRate = 30;
        open(writerObj);                
        
        for j=1:size(data,3)
            data(:,:,j) = medfilt2(data(:,:,j), [2,2]);
        end
        
        data(data<limits(1))=limits(1);
        data(data>limits(2))=limits(2);
        data = data - min(data(:)) + 1e-4;
        data =data/(max(data(:)));
        data(data>CUT_PERCENTAGE)=CUT_PERCENTAGE - 1e-4;
        data = data/CUT_PERCENTAGE;
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
        
        data = [];
        
        close(writerObj);                
        
    end
    
    fclose(fid);
    
end

if nargin==3
    load(cfg_file);    
    cfg.processing_stage = 6;
    save(cfg_file,'cfg');
end

fprintf('all done (took %is)\n\n',round(toc));

end

function limits = get_color_limits(fid,startind,frames,steps,S)

pos = round(linspace(startind,startind+frames,min(frames,steps+3)));

pos = pos(2:(end-1));
mi=0;
ma=0;

for i=1:length(pos)
   fseek(fid,2*S*(pos(i)-1),'bof');
   data = single(fread(fid,S,'*uint16'));
   mi = mi + prctile(data(:),1)/length(pos);
   ma = ma + prctile(data(:),99)/length(pos);
end

limits = [ceil(mi),floor(ma)];

end