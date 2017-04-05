function create_diagnostic_figures(ops0, db,cfg)

if ischar(db) && ischar(ops0)
    load(db);
    load(ops0);
end

ops = build_ops3(db, ops0);

root = ops.ResultsSavePath;
fregops =  sprintf('regops_%s_%s.mat', ops.mouse_name, ops.date);
allops = load(fullfile(root, fregops));
allops = allops.ops1;

tic;

fprintf('\n\nCreating diagnostic images\n');
for plane = 1:ops.nplanes
    
    fprintf('... plane %i\n',plane);
    
    ops1  = allops{plane};
    limits = [0,cumsum(ops1.Nframes)];    
    fpath = sprintf('%s/F_%s_%s_plane%d_proc.mat', ops1.ResultsSavePath, ...
        ops1.mouse_name, ops1.date, plane);
    if exist(fpath, 'file')
        A = load(fpath,'Fcell');
    else
        fpath = sprintf('%s/F_%s_%s_plane%d.mat', ops1.ResultsSavePath, ...
            ops1.mouse_name, ops1.date, plane);
        A = load(fpath,'Fcell');
    end    
    
    mult_align = nan(size(ops1.mimg_end,3),2);
    
    fold_root = [ops.klab_info.outpath,filesep,'diagnostics'];
    if ~exist(fold_root,'dir')
        mkdir(fold_root);
    end          
    
    all_str=[];
    
    for file = 1:size(ops1.mimg_end,3)
        
        [~,str]=fileparts(ops1.klab_info.sbxfiles{file});
        all_str{file} = str;
        titlesrt = sprintf('%s (file %i), plane %i',str,file,plane);
        
        handle =  figure('position',[667,511,1275,987],'Visible','off');
        colormap('gray');
        imagesc(0.5*ops1.mimg_end(:,:,file)+0.5*ops1.mimg_beg(:,:,file));hold on;
        colorbar;
        hold on;
        title(titlesrt,'interpreter','none');
        xlabel('x (pixel)')
        ylabel('y (pixel)')
        rectangle('Position',[ops1.xrange(1),ops1.yrange(1),range(ops1.xrange),range(ops1.yrange)],'EdgeColor','g','FaceColor','none');    
                        
        if ops.nplanes>1
            fold = [fold_root,filesep,'mean_intensity',filesep,['plane_',num2str(plane)]];
        else
            fold = [fold_root,filesep,'mean_intensity'];
        end
        if ~exist(fold,'dir')
            mkdir(fold);
        end
        saveas(handle,[fold,filesep,sprintf('mean_intensity_file%i_plane%i_%s.png',file,plane,str)]);
        close(handle);

        %% motion correction        
        
        frame_ind = (limits(file)+1):limits(file+1);
        N = length(frame_ind);
        
        if ops.nplanes>1
            fold = [fold_root,filesep,'motion_correction',filesep,['plane_',num2str(plane)]];
        else
            fold = [fold_root,filesep,'motion_correction'];
        end        
        if ~exist(fold,'dir')
            mkdir(fold);
        end
        
        mot = ops1.DS(frame_ind,:);
        
        bad_frames = find(ops1.badframes(frame_ind));
        
        mult_align(file,1:2) = mean(mot);
                
        handle = figure('Visible','off','position',[1000,770,680,728],'PaperPositionMode','auto');
        
        subplot(3,1,1);                      
        plot(mot);hold on;
        plot(bad_frames,mot(bad_frames,:),'g*');        
        legend('Y','X','bad','location','best');
        xlabel('Frame');
        ylabel('Shift (pixels)');
        title(titlesrt,'interpreter','none');
        axis tight;
        box on;
        
        %% RMS and correlation
        subplot(3,1,2);
                
        corvals = ops1.CorrFrame(frame_ind,:);  
        ax = plotyy(1:N,sqrt(sum(detrend(mot,'constant').^2,2)),1:N,corvals);
        xlabel('Frame');
        ylabel(ax(1),'Shift norm (pixels)');
        ylabel(ax(2),'Correlation with mean');
        axis(ax,'tight');
        box on;
        
        % mean ROI intensity
        subplot(3,1,3);        
        m = mean(A.Fcell{file},1);
        ax = plotyy(1:N,m,1:N,100*(m-mean(m))/mean(m)); hold on;       
        xlabel('Frame');        
        ylabel(ax(2),'% change from mean');
        ylabel(ax(1),'Mean raw signal');
        axis(ax,'tight');
        box on;
        title(sprintf('Mean ROI (total %i) signal stats',size(A.Fcell{file},1)));        
        
        saveas(handle,[fold,filesep,sprintf('motion_correction_file%i_plane%i_%s.png',file,plane,str)]);        
        close(handle);                        
    end      
    
    handle =  figure('position',[259         642        1833         825],'Visible','off');
    plot(1:size(mult_align,1),mult_align,'o-');hold on;
    plot(1:size(mult_align,1),zeros(1,size(mult_align,1)),'k:');
    axis tight;
    ylabel('Mean relative shift (pixels)');
    legend('Y','X','location','best');
    title(sprintf('Experiment %s (%i files), plane %i',ops1.klab_info.experiment_ID,size(mult_align,1),plane),'interpreter','none');
    set(handle.CurrentAxes,'XTick',1:size(mult_align,1),'XTicklabel',all_str,'XTickLabelRotation',90,'ticklabelinterpreter','none','FontSize',12);    
    saveas(handle,[fold_root,filesep,sprintf('multialignment_plane%i_%s.png',plane,ops1.klab_info.experiment_ID)]);

end

fprintf('all done (took %is)\n\n',round(toc));

if nargin==3
    load(cfg);    
    cfg.processing_stage = 5;
    save(cfg.CONFIGFILE,'cfg');
end

end