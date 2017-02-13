function  run_pipeline(db, ops0)

% ops0.TileFactor (or db(iexp).TileFactor) can be set to multiply the number of default tiles for the neuropil

ops0.nimgbegend                     = getOr(ops0, {'nimgbegend'}, 0);
ops0.splitROIs                      = getOr(ops0, {'splitROIs'}, 1);
ops0.LoadRegMean                    = getOr(ops0, {'LoadRegMean'}, 0);
ops0.NiterPrealign                  = getOr(ops0, {'NiterPrealign'}, 10);
ops0.registrationUpsample           = getOr(ops0, {'registrationUpsample'}, 1);  % upsampling factor during registration, 1 for no upsampling is much faster, 2 may give better subpixel accuracy
ops0.getROIs                        = getOr(ops0, {'getROIs'}, 1);   % whether to run the optimization
ops0.getSVDcomps                    = getOr(ops0, {'getSVDcomps'}, 0);   % whether to save SVD components to disk for later processing
ops0.nSVD                           = getOr(ops0, {'nSVD'}, 1000);   % how many SVD components to save to disk
if isfield(ops0, 'numBlocks') && ~isempty(ops0.numBlocks) && ops0.numBlocks> 1
    ops0.nonrigid                   = 1;
end
ops0.nonrigid                       = getOr(ops0, 'nonrigid', 0);   
ops0.kriging                        = getOr(ops0, 'kriging', 1);  

ops                                 = build_ops3(db, ops0);

if ~isfield(ops, 'diameter') || isempty(ops.diameter)
    warning('you have not specified mean diameter of your ROIs')
    warning('for best performance, please set db(iexp).diameter for each experiment')
end
ops.diameter                        = getOr(ops, 'diameter', 10);
ops.clustrules.diameter             = ops.diameter;
ops.clustrules                      = get_clustrules(ops.clustrules);
%%
% this loads ops1 and checks if processed binary files exist
opath = sprintf('%s/regops_%s_%s.mat', ops.ResultsSavePath, ops.mouse_name, ops.date);
processed = 1;
if exist(opath, 'file')
    load(opath);
    for j = 1:numel(ops1)       
       if ~exist(ops1{j}.RegFile, 'file') % check if the registered binary file exists
          processed = 0; 
       end
    end
else
    processed = 0;
end
%%
% do registration if the processed binaries do not exist
if processed==0
    if ops.nonrigid
        ops1 = blockReg2P(ops);  % do non-rigid registration
    else
        ops1 = reg2P(ops);  % do registration
    end
else
    disp('already registered binary found');
end

if ops.fix_baseline && processed==0
    
    fprintf('\n--- Starting file-wise baseline fix ---\n\n')
    
    for i = 1:numel(ops1)
        
        ops         = ops1{i};
        ops.iplane  = i;
        [Ly, Lx] = size(ops.mimg);
        
        fid = fopen(ops.RegFile, 'r+');
        
        mean_images = zeros(Ly,Lx,length(ops.Nframes),'single');
        corr_images = zeros(Ly,Lx,length(ops.Nframes),'single');
        count_images = zeros(Ly,Lx,length(ops.Nframes),'single');
        
        MAXINT = single(intmax('uint16'));
        
        fprintf('plane %i with %i blocks and %i frames\n',i,length(ops.Nframes),sum(ops.Nframes));
        
        fprintf(' computing mean...\n');
        for block = 1:length(ops.Nframes);
            
            data = fread(fid,  Ly*Lx*ops.Nframes(block), '*uint16');
            data = single(reshape(data, Ly, Lx, []));
            %
            if nnz(data<0)>0
                warning('Total %i pixels are negative (binary %s, block %i) !!!',nnz(data<0),ops.RegFile,block);
            end
            if nnz(data>MAXINT)>0
                warning('Total %i pixels are over uint16 limit (binary %s, block %i) !!!',nnz(data>MAXINT),ops.RegFile,block);
            end
            %
            nans = isnan(data);
            data(nans) = 0;
            % Count up non-NaNs.
            n = sum(~nans,3);
            nans = [];
            n(n==0) = NaN; % prevent divideByZero warnings
            % Sum up non-NaNs, and divide by the number of non-NaNs.
            m = sum(data,3)./n;
            count_images(:,:,block)  = n;
            mean_images(:,:,block) = m;            
            corr_images(:,:,block) = correlation_image(data,12,Ly,Lx);                       
        end
        data=[];
                
        background = 0;
        corrmap = 0;
        coef = sum(count_images,3);
        for block = 1:length(ops.Nframes);
            background = background + mean_images(:,:,block).*count_images(:,:,block)./coef;
            corrmap = corrmap + corr_images(:,:,block).*count_images(:,:,block)./coef;
        end
        count_images=[];
        corr_images=[];
        
        if nnz(isnan(background))>0
           warning('Background has %i NaN''s !!!!',nnz(isnan(background)));
        end
        
        mimg1 = zeros(Ly,Lx);
       % tempfile = [ops.RegFile(1:end-4),'_TEMP.bin'];
        %fid_new = fopen(tempfile, 'w');
        frewind(fid);
        
        fprintf(' rescaling with common mean...\n');
        median_multip = nan(length(ops.Nframes),3);
        for block = 1:length(ops.Nframes)
            
            pos = ftell(fid);
            data = fread(fid,  Ly*Lx*ops.Nframes(block), '*uint16');
            data = single(reshape(data, Ly, Lx, []));
            
            multimat = background./mean_images(:,:,block);
            data = bsxfun(@times,data,multimat);
            median_multip(block,1) = median(multimat(:));
            median_multip(block,2) = prctile(multimat(:),2.5);
            median_multip(block,3) = prctile(multimat(:),97.5);

            if ops.fix_baseline==2
                fprintf('...running mean detrend removal (block %i)\n',block);
                data = klab_detrend(data,ops.imageRate,ops.detrend_window,'MEAN');
            end
            %data = data + shift_levels(block);
            
            data = max(data,0);
            if nnz(data>MAXINT)>0
                warning('Total %i pixels are over uint16 limit after level fix (binary %s, block %i) !!!',nnz(data>MAXINT),ops.RegFile,block);
            end
            
            mimg1 = mimg1 + sum(data,3);
                        
            %count = fwrite(fid_new,uint16(data),'*uint16');
            fseek(fid,pos,'bof');
            count = fwrite(fid,uint16(data),'*uint16');
            
            if count ~= Ly*Lx*ops.Nframes(block)
                error('Number of written elements mismatch! (BUG)');
            end
            
        end
        data=[];
        mean_images=[];
        
        ops.mimg1 = mimg1/sum(ops.Nframes);
        
        fprintf('level-fix multipliers for plane %i:\n',ops.iplane);
        for block = 1:length(ops.Nframes)
            fprintf(' median %3.2f, 95%% interval %3.2f - %3.2f (file %i with %i frames)\n',median_multip(block,1),median_multip(block,2),median_multip(block,3),block,ops.Nframes(block));
        end
        fprintf('\n');
        
        fclose(fid);     
        
        ops.common_mean = background;
        ops.common_corrmap = corrmap;
        ops.median_multip = median_multip;
        
        ops1{i} = ops;
                       
    end    
    
    fprintf('--- signal level fixing finished ---\n\n')
end

%%
for i = 1:numel(ops1)
    ops         = ops1{i};    
    ops.iplane  = i;
    
    if numel(ops.yrange)<10 || numel(ops.xrange)<10
        warning('valid range after registration very small, continuing to next plane')
        continue;
    end
    
    if getOr(ops, {'getSVDcomps'}, 0)
        % extract and write to disk SVD comps (raw data)
        ops    = get_svdcomps(ops);
    end
    
    if ops.getROIs || getOr(ops, {'writeSVDroi'}, 0)
        % extract and/or write to disk SVD comps (normalized data)
        [ops, U, model]    = get_svdForROI(ops);
    end
        
    if ops.getROIs
        % get sources in stat, and clustering images in res
        [ops, stat, model]           = sourcery(ops,U, model);
        
        % extract dF
        [ops, stat, Fcell, FcellNeu,baselines,scalefactors] = extractSignals(ops, model, stat);

        % apply user-specific clustrules to infer stat.iscell
        stat                         = classifyROI(stat, ops.clustrules);
                
%         mimg1 = ops.mimg1(ops.yrange, ops.xrange);
%         baselines1 = nan(1,numel(stat),1);
%         scalefactors1 = baselines1;
%         for k = 1:numel(stat)
%             w = stat(k).lam(:)';
%             ww = w/sum(w);
%             ipix = stat(k).ipix(:)';
%             if ~isempty(ipix)
%                 scalefactors1(k) = ww *ops.sdmov(ipix)';
%                 baselines1(k) = w*mimg1(ipix)';
%             end
%         end        
        
        save(sprintf('%s/F_%s_%s_plane%d.mat', ops.ResultsSavePath, ...
            ops.mouse_name, ops.date, ops.iplane),  'ops',  'stat',...
            'Fcell', 'FcellNeu','scalefactors','baselines','-v7.3');
    end
    
    if ops.DeleteBin
        fclose('all');
        delete(ops.RegFile);        % delete temporary bin file
    end
    
end

% clean up
fclose all;