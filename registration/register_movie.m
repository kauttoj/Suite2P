function dreg = register_movie(data, ops, ds)

orig_class = class(data);

[h, w, nFrames] = size(data);
if getOr(ops, 'useGPU', 0)
    nFramesPerBatch = getBatchSize(w*h);
else
    nFramesPerBatch = 1000; % if not on GPU, should have plenty of RAM available
end
nBatches = ceil(nFrames/nFramesPerBatch);
startFrame = 1:nFramesPerBatch:nFrames;
endFrame = min(startFrame+nFramesPerBatch-1, nFrames);
dreg = zeros(size(data), orig_class);

globalmin = max(min(data(:))*0.95,0);
globalmax = max(data(:))*1.05;

for iBatch = 1:nBatches
    idx = startFrame(iBatch):endFrame(iBatch);
    if ops.useGPU
        dataBatch = gpuArray(single(data(:,:,idx)));
    else
        dataBatch = data(:,:,idx);
    end
    [Ly, Lx, NT] = size(dataBatch);
    
    Ny = ifftshift([-fix(Ly/2):ceil(Ly/2)-1]);
    Nx = ifftshift([-fix(Lx/2):ceil(Lx/2)-1]);
    [Nx,Ny] = meshgrid(Nx,Ny);
    Nx = Nx / Lx;
    Ny = Ny / Ly;
    
    if ops.useGPU
        dregBatch = gpuArray.zeros(size(dataBatch), orig_class);
    else
        dregBatch = zeros(size(dataBatch), orig_class);
    end
    
    if ops.useGPU
        dsBatch = gpuArray(permute(ds(idx, :), [3, 2, 1]));
        Nx = gpuArray(single(Nx));
        Ny = gpuArray(single(Ny));
    else
        dsBatch = ds(idx, :);
    end
    
    if ops.useGPU % do it batch-by-batch
        dph         = 2*pi*(bsxfun(@times, dsBatch(1,1,:), Ny) + ...
            bsxfun(@times, dsBatch(:,2,:), Nx));
        fdata       = fft2(dataBatch);
        dregBatch = real(ifft2(fdata .* exp(1i * dph)));
    else % do it frame-by-frame
        for i = 1:NT
            dph         = 2*pi*(dsBatch(i,1)*Ny + dsBatch(i,2)*Nx);
            fdata       = fft2(single(dataBatch(:,:,i)));
            dregBatch(:,:,i) = real(ifft2(fdata .* exp(1i * dph)));
        end
    end
    
    % limit values into sane ones (sanity check)
    ind = dregBatch < globalmin;
    dregBatch(ind) =  globalmin;
    
    ind = dregBatch > globalmax;
    dregBatch(ind) = globalmax;
    
    if ops.useGPU
        dregBatch = gather_try(dregBatch);
    end
    dreg(:,:,idx) = dregBatch;    
    
end




