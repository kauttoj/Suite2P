function create_diagnostic_figures(ops, db)

ops = build_ops3(db, ops);

root = ops.ResultsSavePath;
fregops =  sprintf('regops_%s_%s.mat', ops.mouse_name, ops.date);
allops = load(fullfile(root, fregops));
allops = allops.ops1;

tic;

fprintf('\n\nCreating diagnostic images\n');
for plane = 1:ops.planesToProcess
    
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
        imagesc(ops1.mimg_end(:,:,file));hold on;
        colorbar;
        %line()
        
        title(titlesrt,'interpreter','none');        
        fold = [fold_root,filesep,'mean_intensity_plane',num2str(plane)];
        if ~exist(fold,'dir')
            mkdir(fold);
        end
        saveas(handle,[fold,filesep,sprintf('mean_intensity_file%i_plane%i_%s.png',file,plane,str)]);
        close(handle);

        %% motion correction        
        
        frame_ind = (limits(file)+1):limits(file+1);
        N = length(frame_ind);
        
        fold = [fold_root,filesep,'motion_correction_plane',num2str(plane)];
        if ~exist(fold,'dir')
            mkdir(fold);
        end        
        
        mot = ops1.DS(frame_ind,:);
        
        mult_align(file,1:2) = mean(mot);
                
        handle = figure('Visible','off','position',[781         491         899        1007],'PaperPositionMode','auto');
        
        subplot(3,1,1);
        plot(detrend(mot));
        legend('Y','X','location','best');
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
        axes(ax(1)); ylabel('Shift norm (pixels)');
        axes(ax(2)); ylabel('Correlation with mean');
        axis(ax,'tight');
        box on;
        
        % mean ROI intensity
        subplot(3,1,3);        
        title(sprintf('Mean ROI (total %i) signal stats',size(A.Fcell{file},1)));        
        m = mean(A.Fcell{file},1);
        ax = plotyy(1:N,m,1:N,100*(m-mean(m))/mean(m)); hold on;       
        xlabel('Frame');        
        axes(ax(2)); ylabel('% change from mean');
        axes(ax(1)); ylabel('Mean raw signal');
        axis(ax,'tight');
        title('mean ROI signal')
        box on;
        
        saveas(handle,[fold,filesep,sprintf('motion_correction_file%i_plane%i_%s.png',file,plane,str)]);        
        close(handle);                        
    end      
    
    handle =  figure('position',[259         642        1833         825],'Visible','off');
    plot(1:size(mult_align,1),mult_align,'o-');hold on;
    plot(1:size(mult_align,1),zeros(1,size(mult_align,1)),'k:');
    axis tight;
    ylabel('Mean relative shift (pixels)');
    legend('Y','X','location','best');
    title(sprintf('Experiment %s (%i files), plane %i',ops1.klab_info.experiment_ID,size(mult_align,1),plane));
    set(handle.CurrentAxes,'XTick',1:size(mult_align,1),'XTicklabel',all_str,'XTickLabelRotation',90,'ticklabelinterpreter','none','FontSize',12);    
    saveas(handle,[fold_root,filesep,sprintf('multialignment_plane%i_%s.png',plane,ops1.klab_info.experiment_ID)]);

end

fprintf('all done (took %is)\n\n',round(toc));


end