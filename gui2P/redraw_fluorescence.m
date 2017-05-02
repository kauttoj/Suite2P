function redraw_fluorescence(h)

axes(h.axes4)

cla(h.axes4,'reset');

blue = [0,0.4470,0.7410];
green = [0,200,50]/255;
red = [0.8500,0.3250,0.0980];

ichosen = h.dat.F.ichosen;
F = [];
Fneu = [];
border = zeros(1,1+numel(h.dat.Fcell));
for j = 1:numel(h.dat.Fcell)      
    F    = cat(2, F, h.dat.Fcell{j}(ichosen, :));
    
    Fneu = cat(2, Fneu, h.dat.FcellNeu{j}(ichosen, :));    
    
    border(j+1) = length(F);
end

F_dff = F - h.dat.stat(ichosen).neuropilCoefficient*Fneu + median(F);

Fneu = Fneu*h.dat.stat(ichosen).neuropilCoefficient;

DS = 0*F_dff;

if h.is_ROI_making_mode==3
    
    if h.custom_ROI.sig_isred
        plot(h.custom_ROI.sig,'color',red,'linewidth',2);
        legend('RED channel signal');
    else
        plot(h.custom_ROI.sig,'color',green,'linewidth',2);
        legend('GREEN channel signal');
    end
    
    axis tight
    
    ylabel('Mask mean signal','Fontsize',12);
    
else
    
    if h.show_dff==-1
        
        for i=1:(length(border)-1)
            ind = ((border(i)+1):border(i+1));
            DS(ind) = sqrt(sum(detrend(h.dat.ops.DS(ind,:)).^2,2));
        end
        plot(DS,'color',blue);
        ylabel('motion correction (pixels)','Fontsize',12);
        
        axis tight
        
    elseif h.show_dff==1 && ~isempty(F_dff)
        
        F_dff_spikes = 0*F_dff;
        F_dff_spikes(h.dat.stat(ichosen).st) = h.dat.stat(ichosen).c;
        
        for i=1:(length(border)-1)
            ind = ((border(i)+1):border(i+1));
            F_temp = F_dff(ind);
            th = prctile(F_temp,80);
            F0 = mean(F_temp(F_temp<th));
            F_dff(ind) = (F_temp-F0)/F0;
        end
        
        x = 1:length(F_dff_spikes);
        
        if strcmp(h.zoom_handle.Enable,'off')        
           F_dff = my_conv_local(medfilt1(double(F_dff), 3), 3);
        end
            
        [AX,H1,H2] = plotyy(x,F_dff_spikes,x,F_dff);
        
        ylabel(AX(1),'Spike amplitude','Fontsize',12);
        ylabel(AX(2),'(F-F0)/F0','Fontsize',12);
        
        axis(AX(2),'tight');
        axis(AX(1),'tight');
        
        ylim_left = ylim(AX(1));
        ylim_right = ylim(AX(2));
        
        % Set same ratio for first axis
        [ylim_left,ylim_right] = fix_zero(ylim_left,ylim_right);
        set(AX(1),'Ylim',ylim_left);
        set(AX(2),'Ylim',ylim_right);
        
        %% CANNOT CHANGE Y-AXES AS UISTACK WILL DISTURB UI EVENTS
        %uistack(AX(2),'bottom')
        %set(AX(1), 'Color', 'none');
        %set(AX(2), 'Color', 'w');                
        
    else
        
        if strcmp(h.zoom_handle.Enable,'off')        
           F = my_conv_local(medfilt1(double(F), 3), 3);
        end
        
        plot(F,'color',blue);
        
        axis tight
        hold on
        
        if isfield(h.dat, 'FcellNeu')
            if strcmp(h.zoom_handle.Enable,'off')      
                Fneu = my_conv_local(medfilt1(double(Fneu), 3), 3);
            end
            plot(Fneu);
        end
        ylabel('Raw signal','Fontsize',12);
        
    end
end

hold on;
y = ylim();
r=range(y);
for i=2:(length(border))
   line(border(i)*[1,1],y,'color',[0.6,0.6,0.6]);
   text(mean(border(i+[-1,0])),y(2)-0.05*r,sprintf('%i',i-1),'horizontalalign','center','color',[0.6,0.6,0.6]);
end

box on

end

function [ylim_left,ylim_right] = fix_zero(ylim_left,ylim_right)
 
pos1 = -ylim_left(1)/range(ylim_left);
pos2 = -ylim_right(1)/range(ylim_right);

pos_new = (pos1+pos2)/2;

if pos1 < pos_new
    s1 = ylim_left(1) - (ylim_left(2)*pos_new)/(pos_new - 1);
    ylim_left(1) = ylim_left(1) - s1;
else
    s1 = -ylim_left(1)/pos_new - range(ylim_left);
    ylim_left(2) = ylim_left(2) + s1;
end
if pos2 < pos_new
    s2 = ylim_right(1) - (ylim_right(2)*pos_new)/(pos_new - 1);
    ylim_right(1) = ylim_right(1) - s2;
else
    s2 = -ylim_right(1)/pos_new - range(ylim_right);
    ylim_right(2) = ylim_right(2) + s2;
end

end