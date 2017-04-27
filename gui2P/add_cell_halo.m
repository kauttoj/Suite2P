function I = add_cell_halo(I,mask)

if nnz(mask)>0
    SE = strel('disk',2);
    mask = imdilate(mask,SE)-mask;
    if size(I,3)==1
        I = mask;
        return;
    end
    mask = mask(:);
    mask = [mask,2*mask,mask];   
    mask = reshape(mask,size(I));
    I(mask==1)=0;
    I(mask==2)=1;
end

end