function [data,means,trend] = klab_detrend(data,framerate,T,TYPE)
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
%TYPE = 'MEAN'; % 'MEAN'

if nargin == 0
    % testdata
    framerate = 15.48;
    %data = max(0,1000+cumsum(randn(ceil(2*60/(1/framerate)),20*34)));
    %data = reshape(data,20,34,[]);
    %data(randsample(numel(data),round(numel(data)*0.2)))=0;
    A = load('example_data');
    data = single(A.data);
    cols = [31918,10735,11334,8647];    
end

if strcmp(TYPE,'VARIANCE')
    TYPE = 2;
    data = max(data,0.1);
elseif strcmp(TYPE,'MEAN')
    TYPE = 1;
else
    error('Unknown Type!')
end

if nargin<3
    T = 50; % seconds
end

ss = size(data);

% compute baseline using these percentiles of data
MI_PRC = 10; % lower percentile (skip possible outliers at the bottom)
MA_PRC = 60; % upper percentile (all peaks should be above this)

dt = 1/framerate;
L = round(T/dt);

L = L + (mod(L,2)==0);

NEED_RESHAPE = 0;
if length(size(data))==3
    data = reshape(data,[],ss(3))';
    NEED_RESHAPE = 1;
end

ZERO_CUTPOINT = 0;
if max(max(data))>=0
    ZERO_CUTPOINT = 1;
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
if nargout>2 || nargin == 0
    GET_TREND = 1;
    trend = zeros(size(data),'single');
end

if timepoints(end)<N_timepoints
    timepoints(end) = N_timepoints;
end

LO_ind = round(L*MI_PRC/100);
HI_ind = round(L*MA_PRC/100);
IND_range = LO_ind:HI_ind;

% 		if(T*cfg.TR>=240)
% 			SGlen=round(240/cfg.TR);
% 		else
% 			SGlen=T; % if we have less than 4 minutes, let's use all the data
% 		end
% 		disp(['Performing Savitzky-Golay detrending over ' num2str(SGlen) ' timepoints']);
%
% 		if(mod(SGlen,2)==0) SGlen=SGlen-1; end % it needs to be odd
% 		trend=sgolayfilt(tempdata,3,SGlen);
% 		for v=1:size(tempdata,2) % foreach voxel
% 			if(var(trend(:,v))==0) continue; end
% 			[aa bb res]=regress(tempdata(:,v),[trend(:,v)  ones(T,1)]);
% 			tempdata(:,v)=res;
% 		end
% 		for row=1:T

%yy2 = smooth(1:N_timepoints,data,L,'rloess');



if N_timepoints < round(L*1.25)
    
    warning('Too few frames (%i points with window %i), skipping detrend!',N_timepoints);
    
    data_orig = data;
    
    tempdata = sort(data,1,'ascend');
    
    LO_ind = round(N_timepoints*MI_PRC/100);
    HI_ind = round(N_timepoints*MA_PRC/100);
    IND_range = LO_ind:HI_ind;
    
    tempdata = tempdata(IND_range,:);
    
    data=repmat(nanmedian(tempdata,1),N_timepoints,1);
    if GET_TREND
        trend = data;
    end
    
    m = mean(data);
    
    data = data_orig - data;
    if TYPE == 1 % MEAN
        data = bsxfun(@plus,data,m);
    elseif TYPE == 2 % VARIANCE
        data = max(0.1,bsxfun(@plus,data,m));
        data = bsxfun(@times,data,1./m);
    else
        error('unknown type')
    end
    
else
    
    tic;
    fprintf('Starting detrend (window length %is, step %is)\n',round(T),stepping*dt);
    for chunk = 1:N_chunks
        
        lims = chunklims(chunk):(chunklims(chunk+1)-1);
        data_orig = data(:,lims);
        medpoints = zeros(length(timepoints),length(lims));
        
        parfor tt=1:length(timepoints)
            
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
            
            tempdata = sort(tempdata,1,'ascend');
            tempdata = tempdata(IND_range,:);
            
            medpoints(tt,:)=nanmedian(tempdata,1);
            
        end
        
        data(timepoints,lims)=medpoints;
        
        data(:,lims) = interp1(timepoints,data(timepoints,lims),1:N_timepoints,'PCHIP');
        
        if GET_TREND
            trend(:,lims) = data(:,lims);
        end
        
        m = mean(data(:,lims));
        
        data(:,lims) = data_orig - data(:,lims);
        if TYPE == 1 % MEAN
            data(:,lims) = bsxfun(@plus,data(:,lims),m);
        elseif TYPE == 2 % VARIANCE
            data(:,lims) = max(0.1,bsxfun(@plus,data(:,lims),m));
            data(:,lims) = bsxfun(@times,data(:,lims),1./m);
        else
            error('unknown type')
        end               
        
%         if mod(chunk,floor(N_chunks/10))==0
%             fprintf('...chunk %i/%i done\n',chunk,N_chunks)
%         end
        
    end
end

if NEED_RESHAPE
    data = reshape(data',ss);
    means = reshape(means,ss(1:2));
    if GET_TREND
        trend = reshape(trend',ss);
    end
end

if ZERO_CUTPOINT==1
    data = max(0,data);  
end

fprintf('All done! (took %is)\n',round(toc));

end



function y = prctile(x,p,dim)
%PRCTILE Percentiles of a sample.
%   Y = PRCTILE(X,P) returns percentiles of the values in X.  P is a scalar
%   or a vector of percent values.  When X is a vector, Y is the same size
%   as P, and Y(i) contains the P(i)-th percentile.  When X is a matrix,
%   the i-th row of Y contains the P(i)-th percentiles of each column of X.
%   For N-D arrays, PRCTILE operates along the first non-singleton
%   dimension.
%
%   Y = PRCTILE(X,P,DIM) calculates percentiles along dimension DIM.  The
%   DIM'th dimension of Y has length LENGTH(P).
%
%   Percentiles are specified using percentages, from 0 to 100.  For an N
%   element vector X, PRCTILE computes percentiles as follows:
%      1) The sorted values in X are taken as the 100*(0.5/N), 100*(1.5/N),
%         ..., 100*((N-0.5)/N) percentiles.
%      2) Linear interpolation is used to compute percentiles for percent
%         values between 100*(0.5/N) and 100*((N-0.5)/N)
%      3) The minimum or maximum values in X are assigned to percentiles
%         for percent values outside that range.
%
%   PRCTILE treats NaNs as missing values, and removes them.
%
%   Examples:
%      y = prctile(x,50); % the median of x
%      y = prctile(x,[2.5 25 50 75 97.5]); % a useful summary of x
%
%   See also IQR, MEDIAN, NANMEDIAN, QUANTILE.

%   Copyright 1993-2015 The MathWorks, Inc.


if ~isvector(p) || numel(p) == 0 || any(p < 0 | p > 100) || ~isreal(p)
    error(message('stats:prctile:BadPercents'));
end

% Figure out which dimension prctile will work along.
sz = size(x);
if nargin < 3
    dim = find(sz ~= 1,1);
    if isempty(dim)
        dim = 1;
    end
    dimArgGiven = false;
else
    % Permute the array so that the requested dimension is the first dim.
    nDimsX = ndims(x);
    perm = [dim:max(nDimsX,dim) 1:dim-1];
    x = permute(x,perm);
    % Pad with ones if dim > ndims.
    if dim > nDimsX
        sz = [sz ones(1,dim-nDimsX)];
    end
    sz = sz(perm);
    dim = 1;
    dimArgGiven = true;
end


% Drop X's leading singleton dims, and combine its trailing dims.  This
% leaves a matrix, and we can work along columns.
nrows = sz(dim);
ncols = numel(x) ./ nrows;
x = reshape(x, nrows, ncols);

x = sort(x,1);
n = sum(~isnan(x), 1); % Number of non-NaN values in each column

% For columns with no valid data, set n=1 to get nan in the result
n(n==0) = 1;

IND = round(nrows*p/100);

y=x(IND,:);

% Reshape Y to conform to X's original shape and size.
szout = sz; szout(dim) = numel(p);
y = reshape(y,szout);

% undo the DIM permutation
if dimArgGiven
    y = ipermute(y,perm);
end

% If X is a vector, the shape of Y should follow that of P, unless an
% explicit DIM arg was given.
if ~dimArgGiven && isvector(x)
    y = reshape(y,size(p));
end

end

function y = interpColsSame(x, p, n)
%INTERPCOLSSAME An aternative approach of 1-D linear interpolation which is
%   faster than using INTERP1Q and can deal with invalid data so long as
%   all columns have the same number of valid entries (scalar n).

% Make p a column vector. Note that n is assumed to be scalar.
if isrow(p)
    p = p';
end

% Form the vector of index values (numel(p) x 1)
r = (p/100)*n;
k = floor(r+0.5); % K gives the index for the row just before r
kp1 = k + 1;      % K+1 gives the index for the row just after r
r = r - k;        % R is the ratio between the K and K+1 rows

% Find indices that are out of the range 1 to n and cap them
k(k<1 | isnan(k)) = 1;
kp1 = bsxfun( @min, kp1, n );

% Use simple linear interpolation for the valid precentages
y = bsxfun(@times, 0.5-r, x(k,:)) + bsxfun(@times, 0.5+r, x(kp1,:));

% Make sure that values we hit exactly are copied rather than interpolated
exact = (r==-0.5);
if any(exact)
    y(exact,:) = x(k(exact),:);
end

end

function y = interpColsDiffer(x, p, n)
%INTERPCOLSDIFFER A simple 1-D linear interpolation of columns that can
%deal with columns with differing numbers of valid entries (vector n).

[nrows, ncols] = size(x);

% Make p a column vector. n is already a row vector with ncols columns.
if isrow(p)
    p = p';
end

% Form the grid of index values (numel(p) x numel(n))
r = (p/100)*n;
k = floor(r+0.5); % K gives the index for the row just before r
kp1 = k + 1;      % K+1 gives the index for the row just after r
r = r - k;        % R is the ratio between the K and K+1 rows

% Find indices that are out of the range 1 to n and cap them
k(k<1 | isnan(k)) = 1;
kp1 = bsxfun( @min, kp1, n );

% Convert K and Kp1 into linear indices
offset = nrows*(0:ncols-1);
k = bsxfun( @plus, k, offset );
kp1 = bsxfun( @plus, kp1, offset );

% Use simple linear interpolation for the valid precentages.
% Note that NaNs in r produce NaN rows.
y = (0.5-r).*x(k) + (0.5+r).*x(kp1);

% Make sure that values we hit exactly are copied rather than interpolated
exact = (r==-0.5);
if any(exact(:))
    y(exact) = x(k(exact));
end

end