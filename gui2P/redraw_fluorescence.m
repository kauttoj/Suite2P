function redraw_fluorescence(h)
axes(h.axes4)
hold off

ichosen = h.dat.F.ichosen;
F = [];
Fneu = [];
border = 0;
for j = 1:numel(h.dat.Fcell)
    F    = cat(2, F, h.dat.Fcell{j}(ichosen, :));
    border = [border,length(F)];
    Fneu = cat(2, Fneu, h.dat.FcellNeu{j}(ichosen, :));
end

plot(my_conv_local(medfilt1(double(F), 3), 3))
axis tight
hold on

if isfield(h.dat, 'FcellNeu')
    plot(my_conv_local(medfilt1(double(Fneu), 3), 3))
end

y = ylim();
r=range(y);
for i=2:(length(border))
   line(border(i)*[1,1],y,'color',[0.6,0.6,0.6]);
   text(mean(border(i+[-1,0])),y(2)-0.05*r,sprintf('%i',i-1),'horizontalalign','center','color',[0.6,0.6,0.6]);
end

box off
set(gca, 'xcolor', 'w')
% plot([0 NT], [0 0], 'k', 'Linewidth', 2)
% axis off
