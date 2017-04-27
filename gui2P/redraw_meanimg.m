
function redraw_meanimg(h)

if h.dat.procmap
    I = h.dat.mimg_proc(:,:,h.dat.map);
else    
    I = h.dat.mimg(:,:,h.dat.map);
end

mu = median(I(:));
sd1 = mean(abs(I(I<mu+1e-7) - mu));
sd2 = 1e-7 + mean(abs(I(I>mu-1e-7) - mu));

axes(h.axes2); 

cm = gray(128);
lims = 5*[-sd1 sd2];
I = I - mu;
I(I>lims(2))=lims(2);
I(I<lims(1))=lims(1);
%im=imagesc(I, mu + 5*[-sd1 sd2]);
I = floor(1+ (I - lims(1)) / (lims(2) -lims(1)) * (size(cm,1)-1));

cm = [cm;[0,1,0]];
if h.add_segment_halo
    mask = zeros(h.dat.cl.Ly,  h.dat.cl.Lx);
    mask(h.dat.stat(h.dat.F.ichosen).ipix)=1;
    I(add_cell_halo(I,mask)>0)=size(cm,1);
end
image(I);
colormap(cm);
xlim([h.dat.xlim]); ylim([h.dat.ylim]);
axis off

%colormap('gray')
axes(h.axes3); 
%imagesc(I, mu + 5*[-sd1 sd2]);
image(I);
colormap(cm);
xlim([h.dat.xlim]); ylim([h.dat.ylim]);
axis off
colormap(cm);

drawnow