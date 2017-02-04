function [data,means,trend] = klab_detrend(data,framerate,T)
% [data,means,trend] = klab_detrend(data,framerate,T)
% detrending and baseline estimation
% Very slow adn dumb but safer and more predictable than using high-pass filter
% since we only consider lowest 5-55% of signal values and avoid peaks. If all signal is used, high peaks can
% create distorted trend and even make signal worse than it was...

% Data is mirrored at boundaries.

% Note: Window should be at least 50s and far outside real neural signal variation.

% ver 1, lots of room for optimization...
% 2/3/2017, Janne Kauttonen

CHUNKSIZE = 10000;

if nargin == 0
    % testdata
   framerate = 15.48;   
   data = cumsum(randn(ceil(3*60/(1/framerate)),512*512));
end

if nargin<3
    T = 100; % seconds
end

ss = size(data);

% compute baseline using these percentiles of data
MI_PRC = 5; % lower percentile (skip possible outlier at the bottom)
MA_PRC = 55; % upper percentile (all peaks should be above this)

dt = 1/framerate;
L = round(T/dt);

L = L + (mod(L,2)==0);

NEED_RESHAPE = 0;
if length(size(data))==3   
   data = reshape(data,[],ss(3))';
   NEED_RESHAPE = 1;
end

N_pixels = size(data,2);
N_chunks = ceil(N_pixels/CHUNKSIZE);
chunklims = round(linspace(1,N_pixels+1,N_chunks+1));
N_chunks = length(chunklims)-1;

N_timepoints = size(data,1);

lim1 = floor(L/2);
lim2 = N_timepoints - lim1;

L_ind = (1:L)-lim1-1;

%data1 = data;

% figure;hold on;
% plot(data(:,10));

stepping = floor(L/10);

timepoints = 1:stepping:N_timepoints;

means = zeros(1,N_pixels);

trend = [];
GET_TREND = 0;
if nargout>2  
    GET_TREND = 1;
    trend = zeros(size(data),'single');    
end

if timepoints(end)<N_timepoints
    timepoints(end) = N_timepoints;
end

if N_timepoints < round(L*1.5)
    
    warning('Too few frames (%i points with window %i), skipping detrend!',N_timepoints,L)
    
    L = N_timepoints;
    
    tempdata = data;
    mi = prctile(tempdata,MI_PRC);
    ma = prctile(tempdata,MA_PRC);
    
    mi = repmat(mi,L,1);
    ma = repmat(ma,L,1);
    
    tempdata(tempdata<mi | tempdata>ma)=nan;
    data=repmat(nanmedian(tempdata,1),N_timepoints,1);
    if GET_TREND
        trend = data;
    end
    
    means = mean(data);
    
else
    
    tic;
    fprintf('Starting detrend (window length %is, stepping %is)\n',round(T),round(T/10));
    for chunk = 1:N_chunks
        lims = chunklims(chunk):(chunklims(chunk+1)-1);
        data_orig = data(:,lims);
        for tt=1:length(timepoints)
            
            t = timepoints(tt);
            
            ind = L_ind+t;
                       
            if t<=lim1 % left edge
                ind1 = -(ind(ind<1))+2;
                ind2 = ind(ind>0);
                tempdata = [(data_orig(ind1,:));data_orig(ind2,:)];
            elseif t>lim2 % right edge
                ind1 = 2*N_timepoints -(ind(ind>N_timepoints)) + 1;
                ind2 = ind(ind<N_timepoints+1);
                tempdata = [data_orig(ind2,:);(data_orig(ind1,:))];
            else % middle
                tempdata = data_orig(ind,:);
            end
            
            mi = prctile(tempdata,MI_PRC);
            ma = prctile(tempdata,MA_PRC);
            
            mi = repmat(mi,L,1);
            ma = repmat(ma,L,1);
            
            tempdata(tempdata<mi | tempdata>ma)=nan;
            
            data(t,lims)=nanmedian(tempdata,1);
            
        end
        
        data(:,lims) = interp1(timepoints,data(timepoints,lims),1:N_timepoints,'PCHIP');
        
        if GET_TREND
            trend(:,lims) = data(:,lims);
        end
        
        m = mean(data(:,lims));
        
        means(lims) = m;
        
        data(:,lims) = data_orig - data(:,lims);
        
        data(:,lims) = bsxfun(@plus, data(:,lims),m);
                
        if mod(chunk,floor(N_chunks/10))==0
            fprintf('...chunk %i/%i done\n',chunk,N_chunks)
        end
    end
end

if NEED_RESHAPE
    data = reshape(data',ss);
    means = reshape(means,ss(1:2));
    if GET_TREND
        trend = reshape(trend',ss);
    end
end

fprintf('All done! (took %is)\n',round(toc));

end

