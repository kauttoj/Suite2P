function varargout = new_main(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @new_main_OpeningFcn, ...
                   'gui_OutputFcn',  @new_main_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

warning('on','all');

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before new_main is made visible.
function new_main_OpeningFcn(hObject, eventdata, h, varargin)
% This function has no output args, see OutputFcn.
% varargin   command line arguments to new_main (see VARARGIN)

%% TOGGLE TO FREEZE CLASSIFIER
FREEZE_CLASSIFIER = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n ---- GUI version 5/5/2017 (klab mod)----\n\n');

axes(h.axes2);
set(gca, 'xtick', [], 'ytick', [])
set(gca, 'xcolor', 'k', 'ycolor', 'k')
axes(h.axes3);
set(gca, 'xtick', [], 'ytick', [])
set(gca, 'xcolor', 'k', 'ycolor', 'k')
axes(h.axes4);
set(gca, 'xtick', [], 'ytick', [])
set(gca, 'xcolor', 'k', 'ycolor', 'k')

zoom_handle = zoom();
h.zoom_handle = zoom_handle;

h.output = hObject;
h.FREEZE_CLASSIFIER = FREEZE_CLASSIFIER;

h.show_dff = 0;
h.is_ROI_making_mode = 0;
h.add_segment_halo = 0;

guidata(hObject, h);

set(gcf, 'units','normalized','outerposition',[0.01 0.05 0.98 0.93]);
% UIWAIT makes new_main wait for user response (see UIRESUME)
% uiwait(h.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = new_main_OutputFcn(hObject, eventdata, h) 
% varargout  cell array for returning output args (see VARARGOUT);
varargout{1} = h.output;

function pushbutton87_KeyPressFcn(hObject, eventdata, h)

% main LOAD call
function pushbutton17_Callback(hObject, eventdata, h)
% [filename1,filepath1]=uigetfile('\\zserver\Lab\Share\Marius\', 'Select Data File');

flag = 0;
try
    if isfield(h, 'dat') && isfield(h.dat, 'filename')
        root = fileparts(h.dat.filename);
    else
        root = 'D:\DATA\F\';
    end
    [filename1,filepath1]=uigetfile(fullfile(root, 'F*.mat'), 'Select Data File');
    set(h.figure1, 'Name', filename1);
    
    % construct dat. Everything is loaded except F and Fneu.
    refresh_stats(h,2);
    fprintf('\nLoading file %s ...',fullfile(filepath1, filename1));
    h.dat = load(fullfile(filepath1, filename1));
    fprintf(' done\n');
    
    flag = 1;
catch err
    fprintf('\nFile load failed!! Reason: %s\n',err.message);
end

if flag
    % if the user selected a file, do all the initializations
    rng('default')
    
    % keyboard;
    if isfield(h.dat, 'dat')
        
        h.dat = h.dat.dat;        
        if isfield(h.dat,'classifier_backup')
            h.statLabels = h.dat.classifier_backup.statLabels;
            h.is_shared_classifier = 0;
            rootS2p_PATH=fileparts(which('run_pipeline'));
            run([rootS2p_PATH,filesep,'SHARED_CLASSIFIER_PATHS.m']);
            if exist('CLASSIFIER_DATAFILE','var')
                for i=1:length(CLASSIFIER_DATAFILE)
                    if strcmp(h.dat.cl.fpath,CLASSIFIER_DATAFILE(i).file)
                        h.is_shared_classifier = 1;
                        break;
                    end
                end
            end
            h = identify_classifier(h);
            old_stat = h.dat.stat;
            h = classROI(h);            
            h.dat.stat = old_stat;  % do not use classified prediction, use existing          
        end
        h.quadvalue = zeros(3);
        for j = 1:3
            for i = 1:3
                set(h.(sprintf('Q%d%d', j,i)), 'BackgroundColor',[.92 .92 .92]);
            end
        end
        h.dat.ylim = [0 h.dat.cl.Ly];
        h.dat.xlim = [0 h.dat.cl.Lx];        
        
    else
        h.dat.filename = fullfile(filepath1, filename1);
        h.dat.cl.Ly       = numel(h.dat.ops.yrange);
        h.dat.cl.Lx       = numel(h.dat.ops.xrange);
        
        % make up iclut here
        try
            [h.dat.res.iclust, h.dat.res.lambda, h.dat.res.lambda0] =...
                getIclust(h.dat.stat, h.dat.cl);
        catch
        end
        h.dat.res.iclust = reshape(h.dat.res.iclust, h.dat.cl.Ly, h.dat.cl.Lx);
        %     h.dat.res.lambda = reshape(h.dat.res.lambda, h.dat.cl.Ly, h.dat.cl.Lx);
        
        h.dat.ops.Nk = numel(h.dat.stat);
        h.dat.cl.rands_orig   = .1 + .8 * rand(1, h.dat.ops.Nk);
        h.dat.cl.rands        = h.dat.cl.rands_orig;
        
        if isfield(h.dat.ops, 'clustrules')
            h.dat.clustrules = h.dat.ops.clustrules;
        end
        
        % set up classifier
        h.dat.cl.threshold  = 0.5;
        h.dat.cl.red_threshold  = nan;
        h                   = identify_classifier(h);
        h                   = classROI(h);
        
        % set all quadrants as not visited
        h.quadvalue = zeros(3);
        for j = 1:3
            for i = 1:3
                set(h.(sprintf('Q%d%d', j,i)), 'BackgroundColor',[.92 .92 .92]);
            end
        end
        
        h.dat.ylim = [0 h.dat.cl.Ly];
        h.dat.xlim = [0 h.dat.cl.Lx];
        
        if ~isfield(h.dat.stat,'redcell')
            for j = 1:numel(h.dat.stat)
                h.dat.stat(j).redcell = 0;
                h.dat.stat(j).redprob = 0;
            end
        end
        
        h.dat.F.ichosen = 1;
        
        % loop through redcells and set h.dat.cl.rands(h.dat.F.ichosen) = 0
        for j = find([h.dat.stat.redcell])
            h.dat.cl.rands(j) = 0;
        end
        
        % x and y limits on subquadrants
        h.dat.figure.x0all = round(linspace(0, 19/20*h.dat.cl.Lx, 4));
        h.dat.figure.y0all = round(linspace(0, 19/20*h.dat.cl.Ly, 4));
        h.dat.figure.x1all = round(linspace(1/20 * h.dat.cl.Lx, h.dat.cl.Lx, 4));
        h.dat.figure.y1all = round(linspace(1/20 * h.dat.cl.Ly, h.dat.cl.Ly, 4));
        
    end
    
    % activate all pushbuttons
    pb = [84 93 101 86 87 89 90 92 103 98 95 96 102 99 100 1 2 104 113 114 115];
    for j = 1:numel(pb)
        set(eval(sprintf('h.pushbutton%d', pb(j))),'Enable','on')
    end
    pb = [11 12 13 21 22 23 31 32 33];
    for j = 1:numel(pb)
        set(eval(sprintf('h.Q%d', pb(j))),'Enable','on')
    end
    set(h.full,'Enable', 'on');
    set(h.edit50,'Enable', 'on');
    set(h.edit51,'Enable', 'on');
    set(h.edit50,'String', num2str(h.dat.cl.threshold));
    if ~isfield(h.dat.cl,'red_threshold')
       h.dat.cl.red_threshold = 0.5;
    end
    set(h.edit51,'String', num2str(h.dat.cl.red_threshold));
    set_Bcolor(h, 1);
    set_maskCcolor(h, 1);
    % select unit normalized ROI brightness
    h.dat.cl.vmap = 'unit';
    set_maskBcolor(h, 1);
    set(h.full, 'BackgroundColor', [1 0 0])
    
    % setup different views of GUI
    h.dat.maxmap = 2;
    ops = h.dat.ops;
    if isfield(ops, 'mimg1') && ~isempty(ops.mimg1)
        h.dat.mimg(:,:,h.dat.maxmap) = ops.mimg1(ops.yrange, ops.xrange);
        h.dat.mimg_proc(:,:,h.dat.maxmap) = normalize_image(h.dat.mimg(:,:,h.dat.maxmap));
    end
    h.dat.mimg(:,:,5) = 0;
    
    h.dat.maxmap = h.dat.maxmap + 1;
    if isfield(ops, 'mimgRED') && ~isempty(ops.mimgRED)
        h.dat.mimg(:,:,h.dat.maxmap) = ops.mimgRED(ops.yrange, ops.xrange);
        h.dat.mimg_proc(:,:,h.dat.maxmap) = normalize_image(h.dat.mimg(:,:,h.dat.maxmap));
    end
    h.dat.maxmap = h.dat.maxmap + 1;
    if isfield(ops, 'mimgREDcorrected') && ~isempty(ops.mimgREDcorrected)
        h.dat.mimg(:,:,h.dat.maxmap) = ops.mimgREDcorrected;
        h.dat.mimg_proc(:,:,h.dat.maxmap) = normalize_image(h.dat.mimg(:,:,h.dat.maxmap));
    end
    h.dat.maxmap = h.dat.maxmap + 1;
    if isfield(ops, 'Vcorr') && ~isempty(ops.Vcorr)
        h.dat.mimg(:,:,h.dat.maxmap) = ops.Vcorr;
        h.dat.mimg_proc(:,:,h.dat.maxmap) = normalize_image(h.dat.mimg(:,:,h.dat.maxmap));
    end
    if isfield(ops, 'common_corrmap') && ~isempty(ops.common_corrmap)
        h.dat.mimg(:,:,h.dat.maxmap) = ops.common_corrmap(ops.yrange,ops.xrange);
        h.dat.mimg_proc(:,:,h.dat.maxmap) = normalize_image(h.dat.mimg(:,:,h.dat.maxmap));
    end
    
    h.dat.procmap = 0;
    h.dat.map = 1;
    
    N_cells = 0;
    for i=1:length(h.dat.stat)
        N_cells = N_cells + double(h.dat.stat(i).iscell);
    end
    h.N_cells = N_cells;
    
    redraw_fluorescence(h);
    redraw_figure(h);
    
    hManager = uigetmodemanager(h.figure1);
    try
        set(hManager.WindowListenerHandles, 'Enable', 'off');  % HG1
    catch
        [hManager.WindowListenerHandles.Enabled] = deal(false);  % HG2
    end
%     set(h.figure1,'WindowButtonDownFcn', @(x,y) myCallbackFcn(x,y,h));
%     set(h.figure1,'KeyPressFcn',@(x,y) myCallbackFcn(x,y,h));
    
    guidata(hObject,h)
    
    refresh_stats(h,1);
end

function pushbutton2_Callback(hObject, eventdata, h)
% variance explained mask
h.dat.cl.vmap = 'unit';
set_maskBcolor(h, 1)

redraw_figure(h);
guidata(hObject,h);

function set_maskBcolor(h, ih)
set(h.pushbutton1, 'BackgroundColor', .94 * [1 1 1]); 
set(h.pushbutton2, 'BackgroundColor', .94 * [1 1 1]); 

switch ih
    case 1
        set(h.pushbutton2, 'BackgroundColor', [1 0 0]); 
    case 2
        set(h.pushbutton1, 'BackgroundColor', [1 0 0]); 
end

function pushbutton1_Callback(hObject, eventdata, h)
% unit vector mask
h.dat.cl.vmap = 'var';
set_maskBcolor(h, 2)

redraw_figure(h);
guidata(hObject,h);

function pushbutton84_Callback(hObject, eventdata, h)
% save proc file and rules file
DATALIMIT = 200000;

h.dat.F.trace = [];
dat = h.dat;

h.st0(:,1) = double([h.dat.stat.iscell]);
st = cat(1, h.st, h.st0);
S = size(st,1);
if S>DATALIMIT
    warning('Classification datafile has %i samples, saving only latest %i!',S,DATALIMIT)
    st = st((end-DATALIMIT+1):end,:);
end    

dat.classifier_backup.st = st;
dat.classifier_backup.prior = h.prior;
dat.classifier_backup.statLabels = h.statLabels;
dat.classifier_backup.timestamp = char(datetime('now'));
dat.classifier_backup.classifier_path = h.dat.cl.fpath;
dat.classifier_backup.DATALIMIT = DATALIMIT;

refresh_stats(h,2);
try
    fprintf('\nSaving PROC file %s ...',[h.dat.filename(1:end-4) '_proc.mat']);
    save([h.dat.filename(1:end-4) '_proc.mat'],'dat','-v7.3')
    fprintf(' done\n'); 
catch err    
    warning('\nFailed to write PROC file %s: %s',[h.dat.filename(1:end-4) '_proc.mat'],err.message);
end
refresh_stats(h,1);
%
%
statLabels  = h.statLabels;
prior       = h.prior;
st          = cat(1, h.st, h.st0);
refresh_stats(h,2);

if h.FREEZE_CLASSIFIER==0
    
    try
        fprintf('\nUpdating classifier (total %i samples) file %s ...',size(st,1),h.dat.cl.fpath);
        save(h.dat.cl.fpath, 'st', 'statLabels', 'prior')
        fprintf(' done\n');
    catch err
        warning('\nFailed to write classifier file %s: %s',h.dat.cl.fpath,err.message);
    end
    
else
    
    fprintf('\n NOTE: Classifier updating is turned off (flag FREEZE_CLASSIFIER==1)\n\n')
    
end

refresh_stats(h,1);


function figure1_ResizeFcn(hObject, eventdata, h)

function Q11_Callback(hObject, eventdata, h)
iy = 1; ix = 1;
quadrant(hObject, h, iy, ix)
paint_quadbutton(h, iy, ix);

function Q12_Callback(hObject, eventdata, h)
iy = 1; ix = 2;
quadrant(hObject, h, iy, ix)
paint_quadbutton(h, iy, ix);

function Q13_Callback(hObject, eventdata, h)
iy = 1; ix = 3;
quadrant(hObject, h, iy, ix)
paint_quadbutton(h, iy, ix);

function Q21_Callback(hObject, eventdata, h)
iy = 2; ix = 1;
quadrant(hObject, h, iy, ix)
paint_quadbutton(h, iy, ix);

function Q22_Callback(hObject, eventdata, h)
iy = 2; ix = 2;
quadrant(hObject, h, iy, ix)
paint_quadbutton(h, iy, ix);

function Q23_Callback(hObject, eventdata, h)
iy = 2; ix = 3;
quadrant(hObject, h, iy, ix)
paint_quadbutton(h, iy, ix);

function Q31_Callback(hObject, eventdata, h)
iy = 3; ix = 1;
quadrant(hObject, h, iy, ix)
paint_quadbutton(h, iy, ix);

function Q32_Callback(hObject, eventdata, h)
iy = 3; ix = 2;
quadrant(hObject, h, iy, ix)
paint_quadbutton(h, iy, ix);

function Q33_Callback(hObject, eventdata, h)
iy = 3; ix = 3;
quadrant(hObject, h, iy, ix);
paint_quadbutton(h, iy, ix);

function paint_quadbutton(h, iy, ix)
set(h.full, 'BackgroundColor', .92 * [1 1 1])

for j = 1:3
    for i = 1:3
        if h.quadvalue(j,i)==1
            set(h.(sprintf('Q%d%d', j,i)), 'BackgroundColor', [.4 .4 .4]); 
        end
    end
end
set(h.(sprintf('Q%d%d', iy,ix)), 'BackgroundColor', [1 0 0]); 

function full_Callback(hObject, eventdata, h)
h.dat.ylim = [0 h.dat.cl.Ly];
h.dat.xlim = [0 h.dat.cl.Lx];
if h.dat.map==1
    redraw_figure(h);
else
    redraw_meanimg(h);
end
set(h.full, 'BackgroundColor', [1 0 0]);
for i = 1:3
    for j = 1:3
        if h.(sprintf('Q%d%d', j,i)).BackgroundColor(1) >.99 
            set(h.(sprintf('Q%d%d', j,i)), 'BackgroundColor', [.4 .4 .4])
        end
    end
end
guidata(hObject,h);

function quadrant(hObject, h, iy, ix)
h.dat.ylim = [h.dat.figure.y0all(iy) h.dat.figure.y1all(iy+1)];
h.dat.xlim = [h.dat.figure.x0all(ix) h.dat.figure.x1all(ix+1)];
h.quadvalue(iy, ix) = 1;

guidata(hObject,h);

refresh_figure(h);

function figure1_WindowKeyPressFcn(hObject, eventdata, h)
switch eventdata.Key
    case 'f'
        % flip currently selected unit
        h.dat.stat(h.dat.F.ichosen).iscell = 1 - ...
            h.dat.stat(h.dat.F.ichosen).iscell;
        if h.dat.maxmap==1
            redraw_figure(h);
        end
        guidata(hObject,h);
    case 'q'
        pushbutton87_Callback(hObject, eventdata, h);
    case 'e'
        pushbutton89_Callback(hObject, eventdata, h);
    case 'r'
        if size(h.dat.mimg, 3)>2
            pushbutton90_Callback(hObject, eventdata, h);
        end
    case 't'
        pushbutton92_Callback(hObject, eventdata, h);
    case 'w'
        pushbutton103_Callback(hObject, eventdata, h);
    case 'p'
        pushbutton86_Callback(hObject, eventdata, h);
    case 'a'
        pushbutton98_Callback(hObject, eventdata, h);
    case 's'
        pushbutton95_Callback(hObject, eventdata, h);
    case 'd'
        pushbutton96_Callback(hObject, eventdata, h);
    case 'z'
        pushbutton102_Callback(hObject, eventdata, h);
    case 'x'
        pushbutton99_Callback(hObject, eventdata, h);
    case 'c'
        pushbutton100_Callback(hObject, eventdata, h);
    case 'v'
        pushbutton104_Callback(hObject, eventdata, h);
end

function figure1_WindowButtonDownFcn(hObject, eventdata, h)
z = round(eventdata.Source.CurrentAxes.CurrentPoint(1,:));
x = round(z(1));
y  = round(z(2));

if x>=1 && y>=1 && x<=h.dat.cl.Lx && y<=h.dat.cl.Ly 
    
    if h.is_ROI_making_mode>0        
        ratio = h.dat.cl.Lx/h.dat.cl.Ly;
        if h.is_ROI_making_mode==1 || ...
                ( h.is_ROI_making_mode>1 && norm(h.custom_ROI.Center-[x,y])>h.custom_ROI.Radius )
            h.custom_ROI.Center = [x,y];
            h.custom_ROI.Radius = 10;%h.circleRadius;            
            r = h.custom_ROI.Radius;
            try
                delete(h.custom_ROI.plothandle);
            catch
                
            end                
            h.custom_ROI.plothandle = rectangle('Position',[h.custom_ROI.Center-[r,r/ratio],2*r,2*r/ratio],'Curvature',[1,1],'EdgeColor','g','LineWidth',3);
            h.is_ROI_making_mode = 2;
        else
            r = h.custom_ROI.Radius;   
            try
                delete(h.custom_ROI.plothandle);
            catch
                
            end            
            h.custom_ROI.plothandle = rectangle('Position',[h.custom_ROI.Center-[r,r/ratio],2*r,2*r/ratio],'Curvature',[1,1],'EdgeColor','r','LineWidth',4);            
            h = getROISignal(h);
            h.is_ROI_making_mode = 3;
            redraw_fluorescence(h);
            refresh_figure(h);
        end
        guidata(hObject,h);
        return;
    end
    
    if ~(h.dat.res.iclust(y,x)>0)
        return;
    end        
    
    h.dat.F.ichosen = h.dat.res.iclust(y, x);
    ichosen = h.dat.F.ichosen;
    
    switch eventdata.Source.SelectionType
        case 'alt'
            % flip currently selected unit
            h.dat.stat(ichosen).iscell = 1 - ...
                h.dat.stat(ichosen).iscell;
            
            if h.dat.stat(ichosen).iscell==1
                h.N_cells = h.N_cells + 1;
            else
                h.N_cells = h.N_cells - 1;
            end
            refresh_stats(h,1);
        case 'extend'
            h.dat.stat(ichosen).redcell = 1 -  h.dat.stat(ichosen).redcell;
            
            if h.dat.stat(ichosen).redcell ==1
                h.dat.cl.rands(ichosen) = 0;
            else
                h.dat.cl.rands(ichosen) = h.dat.cl.rands_orig(ichosen);
            end
%             h.dat.img1.H       = reshape(h.dat.cl.rands(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);
%             h.dat.img2.H       = reshape(h.dat.cl.rands(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);
            
            if h.dat.stat(ichosen).redcell
                display('red')
            else
                display('not red')
            end
    end
    
    redraw_fluorescence(h);
    
%     Sat = ones(h.dat.cl.Ly, h.dat.cl.Lx);
%     Sat(h.dat.res.iclust==h.dat.F.ichosen) = 0;
%     h.dat.img1.Sat     = Sat;
%     h.dat.img2.Sat     = Sat;
%     h = buildLambdaValue(h);
    
    if h.add_segment_halo==0
        h.dat.map = 1;
        set_Bcolor(h);
    end

    refresh_figure(h);        
    
    guidata(hObject,h);
       
    str = [];
    labels = [h.statLabels(2:end), {'iscell'}, {'redcell'}];
    for j =1:length(labels)
       if isfield(h.dat.stat, labels{j})
           sl = eval(sprintf('h.dat.stat(ichosen).%s', labels{j}));
           strnew = sprintf('%s = %2.2f \n', labels{j}, sl);
           str = cat(2, str, strnew);
       end
    end
    if h.dat.stat(ichosen).iscell
        str = cat(2, str, sprintf('Segment %i (cell ID %i)',ichosen,sum([h.dat.stat(1:ichosen).iscell])));
    else
        str = cat(2, str, sprintf('Segment %i',ichosen));
    end
    str = cat(2, str, sprintf('\n(x,y)=(%i,%i)',round(h.dat.stat(ichosen).med(2)),round(h.dat.stat(ichosen).med(1))));
    
    set(h.text54,'String', str);
    
end

function pushbutton86_Callback(hObject, eventdata, h)

h.dat.procmap = 1 -  h.dat.procmap;
if h.dat.map>1
    redraw_meanimg(h);
end

if h.dat.procmap>0
    set(h.pushbutton86, 'BackgroundColor', [1 0 0]); 
else
    set(h.pushbutton86, 'BackgroundColor', .94 * [1 1 1]); 
end
guidata(hObject,h);

function set_Bcolor(h, ih)
pb = [87 103 89 90 92];

if nargin==1
    ih = h.dat.map;
end
 
for j = 1:length(pb)
    if j==ih
        set(h.(sprintf('pushbutton%d', pb(ih))), 'BackgroundColor', [1 0 0]); 
    else
        if h.(sprintf('pushbutton%d', pb(j))).BackgroundColor(1)>.99
            set(h.(sprintf('pushbutton%d', pb(j))), 'BackgroundColor', .94 * [1 1 1]);
        end
    end
end

function set_maskCcolor(h, ih)
pb = [98 95 96 102 99 100 104]; 
 
%set_Bcolor(h, 1)

for j = 1:length(pb)
    if j==ih
        set(h.(sprintf('pushbutton%d', pb(ih))), 'BackgroundColor', [1 0 0]); 
    else
        if h.(sprintf('pushbutton%d', pb(j))).BackgroundColor(1)>.99
            set(h.(sprintf('pushbutton%d', pb(j))), 'BackgroundColor', .94 * [1 1 1]);
        end
    end
end


function pushbutton102_Callback(hObject, eventdata, h)
hval = [h.dat.stat.mimgProjAbs];
h.dat.cl.rands   = .1 + .8 * min(1, hval/mean(hval));
h.dat.cl.rands_orig = h.dat.cl.rands;
refresh_figure(h);
set_maskCcolor(h, 4);
guidata(hObject,h);

function pushbutton103_Callback(hObject, eventdata, h)
h.dat.map = 5;
redraw_meanimg(h);
set_Bcolor(h, 2);
guidata(hObject,h);

function pushbutton96_Callback(hObject, eventdata, h)
hval = [h.dat.stat.skew];
h.dat.cl.rands   = .1 + .8 * min(1, hval/mean(hval));
h.dat.cl.rands_orig = h.dat.cl.rands;
refresh_figure(h);
set_maskCcolor(h, 3);
guidata(hObject,h);

function pushbutton95_Callback(hObject, eventdata, h)
h.dat.cl.rands   = .1 + .8 * [h.dat.stat.cellProb];
h.dat.cl.rands_orig = h.dat.cl.rands;
refresh_figure(h);
set_maskCcolor(h, 2);
guidata(hObject,h);

function pushbutton98_Callback(hObject, eventdata, h)
h.dat.cl.rands   = .1 + .8 * rand(1, h.dat.ops.Nk);
h.dat.cl.rands(logical([h.dat.stat.redcell])) = 0;
h.dat.cl.rands_orig = h.dat.cl.rands;
refresh_figure(h);
set_maskCcolor(h, 1);
guidata(hObject,h);

function pushbutton99_Callback(hObject, eventdata, h)
hval = max(0, [h.dat.stat.cmpct]-1);
h.dat.cl.rands   = .1 + .8 * min(1, .5 * hval/mean(hval));
h.dat.cl.rands_orig = h.dat.cl.rands;
refresh_figure(h);
set_maskCcolor(h, 5);
guidata(hObject,h);

function pushbutton100_Callback(hObject, eventdata, h)
hval = [h.dat.stat.footprint];
h.dat.cl.rands   = .1 + .8 * min(1, .5 * hval/mean(hval));
h.dat.cl.rands_orig = h.dat.cl.rands;
refresh_figure(h);
set_maskCcolor(h, 6);
guidata(hObject,h);

% --- Executes on button press in pushbutton104.
function pushbutton104_Callback(hObject, eventdata, h)
hval = [h.dat.stat.redprob];
if isnan(h.dat.cl.red_threshold)
    if sum([h.dat.stat.redcell]) > 0
        Th = min(hval(logical([h.dat.stat.redcell])));
    else
        Th = 1;
    end
    for i=1:length(h.dat.stat)
        h.dat.stat(i).redcell = hval(i)>0.5;
    end        
else
    hval = hval>=h.dat.cl.red_threshold;
    Th = 1;    
    for i=1:length(h.dat.stat)
        h.dat.stat(i).redcell = hval(i);
    end    
end
h.dat.cl.rands   = max(0, min(1, .7 - .7*hval/Th));
h.dat.cl.rands_orig = h.dat.cl.rands;
redraw_figure(h);
set_maskCcolor(h, 7);
guidata(hObject,h);


function pushbutton87_Callback(hObject, eventdata, h)
 h.dat.map = 1;
redraw_figure(h);
set_Bcolor(h, 1);

guidata(hObject,h);

function pushbutton89_Callback(hObject, eventdata, h)
h.dat.map = 2;
redraw_meanimg(h);
set_Bcolor(h, 3);
guidata(hObject,h);

function pushbutton90_Callback(hObject, eventdata, h)
 h.dat.map = 3;
 redraw_meanimg(h);
set_Bcolor(h, 4);
 guidata(hObject,h);
 
% RED CORRECTED BUTTON
function pushbutton92_Callback(hObject, eventdata, h)
h.dat.map = 4;
redraw_meanimg(h);
set_Bcolor(h, 5);
guidata(hObject,h);

function pushbutton93_Callback(hObject, eventdata, h)
rootS2p = which('run_pipeline');
if isempty(rootS2p)
    error('could not identify Suite2p location! where is new_main.m?')
end
rootS2p_PATH=fileparts(rootS2p);
rootS2p = fileparts(rootS2p);
rootS2p = fullfile(rootS2p, 'configFiles');

[filename1,filepath1]   = uigetfile(fullfile(rootS2p, '*.mat'), 'Select classifier file');
if filename1
    
    run([rootS2p_PATH,filesep,'SHARED_CLASSIFIER_PATHS.m']);
    if exist('CLASSIFIER_DATAFILE','var')        
        h.is_shared_classifier = 0;
        for i=1:length(CLASSIFIER_DATAFILE)
            if strcmp(filename1,CLASSIFIER_DATAFILE(i).file)
                h.is_shared_classifier = 1;
                break;
            end
        end
    end    
    
    h.dat.cl.fpath          = fullfile(filepath1, filename1);
    h                       = classROI(h);
    
    redraw_figure(h);
    
%     hload = load(h.dat.cl.fpath);
%     h.st = hload.st;
%     h.save_cls = h.dat.cl.fpath;
%     h.prior = prior;

    guidata(hObject,h);
end

function edit50_Callback(hObject, eventdata, h)
h.dat.cl.threshold = str2double(get(h.edit50,'String'));
for j = 1:length(h.dat.stat)
    h.dat.stat(j).iscell = h.dat.stat(j).cellProb > h.dat.cl.threshold;
end
h.N_cells = sum([h.dat.stat.iscell]);
redraw_figure(h);
refresh_stats(h,1);
guidata(hObject,h);


function edit50_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbutton101_Callback(hObject, eventdata, h)
% new classifier button
rootS2p = which('run_pipeline');
if isempty(rootS2p)
    warndlg('Suite2p location not in PATH! where is run_pipeline.m?')
    [filename1,rootS2p]=uigetfile(root, 'Select run_pipeline.m');
    
end
rootS2p = fileparts(rootS2p);
rootS2p = fullfile(rootS2p, 'configFiles');

def_name = fullfile(rootS2p, 'cl_new.mat');

[filename1,filepath1]   = uigetfile(fullfile(rootS2p, 'priors', '*.mat'), ...
        'Select priors file, will be used to seed the new classifier');
if filename1
    prior_file = fullfile(filepath1, filename1);
else
    error('you must choose a priors file')
end

[FileName,PathName] = uiputfile('*.mat', 'Create new classifier', def_name); 

if FileName
    load(prior_file);
    st                      = [];
    save(fullfile(PathName, FileName), 'st', 'prior', 'statLabels')
    
    h.dat.cl.fpath          = fullfile(PathName, FileName);
    h                       = classROI(h);
    
    redraw_figure(h);

    hload = load(h.dat.cl.fpath);
    h.st = hload.st;
    h.save_cls = h.dat.cl.fpath;
    
    guidata(hObject,h);
else
    warndlg('you did not create a new classifier file');
    error('you did not create a new classifier file')
end

function pushbutton64_Callback(hObject, eventdata, h)

msg{1} = ['First, you need to load an F_*_*_*.mat produced by Suite2p. ']; 
msg{3} = ['Now start exploring the data with the various visualizations. Press the help buttons next to each category to learn the options.'];
msg{5} = ['You can label an ROI as cell/not cell, which determines the image axes it is shown in.'];
msg{7} = ['Right-click on an ROI to flip its label.'];
msg{9} = ['The more you do this, the better the automated classification becomes.'];
msg{11} = ['hint: you can keep any help box open while interacting with the GUI. '];

msgbox(msg, 'Instructions!');

function pushbutton108_Callback(hObject, eventdata, handles)
% msg{1} = ['Timecourse of selected ROI and neuropil fluorescence across experiment.'];
% msg{3} = ['The selected cell has 0 saturation (gray) in ROI image.'];
% msg{5} = ['Traces have been smoothed for display purposes.'];
% msg{7} = ['Statistics shown are same variables as used in "mask color" section.'];
% msgbox(msg, 'Fluorescence instructions');
handles.show_dff = 1-(handles.show_dff==1);
guidata(hObject,handles);
redraw_fluorescence(handles);
refresh_figure(handles);

function pushbutton106_Callback(hObject, eventdata, handles)
msg{1} = ['This applies only if you select "ROIs" under "background".'];
msg{3} = ['Selection determines the color/hue for each ROI, between 0.1 and 0.9.'];
msg{5} = ['Reds are reserved to indicate labels made on the red/secondary color channel (middle mouse-click on an ROI to switch ON/OFF).'];
msg{7} = ['RANDOM: color chosen randomly'];
msg{9} = ['CLASSIFIER: probability assigned by classifier'];
msg{11} = ['SKEW: skewness of activity, after neuropil correction and some smoothing'];
msg{13} = ['MEANIMG: weighting of activity mask onto mean image'];
msg{15} = ['CMPCT: compactness of ROI pixels. Smallest is 1, for disks.'];
msg{17} = ['FOOT: "footprint" of ROI; ~ number of correlated neighboring pixels'];
msg{19} = ['RED: probability of being a red-tagged cell, assigned by algorithm'];
msg{21} = ['hint: the letters in paranthesis are keyboard shortcuts.'];

msgbox(msg, 'Mask color instructions');

function pushbutton107_Callback(hObject, eventdata, handles)
msg{1} = ['This applies only if you select "ROIs" under "background".'];
msg{3} = ['Selection determines the brightness for all pixels inside ROIs.'];
msg{5} = ['UNIT NORM: values are normalized to unit norm per ROI'];
msg{7} = ['VARIANCE: fraction of variance explained in each pixel by parent ROI'];

msgbox(msg, 'Mask brightness instructions');

function pushbutton105_Callback(hObject, eventdata, handles)
msg{1} = ['ROI: shows the masks identified by Suite2p, colored according to the property selected under "mask color"'];
msg{3} = ['CORR: shows the correlation map between a pixel and its nearby pixels'];
msg{5} = ['MEAN: average registered image'];
msg{7} = ['RED: average registered image of "red"/secondary color channel'];
msg{9} = ['RED - GREEN: subtracts off the contamination from "green"/primary channel'];
msg{11} = ['PROC: switch to toggle image contrast normalization'];
msg{13} = ['hint: the letters in paranthesis are keyboard shortcuts'];

msgbox(msg, 'Background instructions');


function pushbutton110_Callback(hObject, eventdata, handles)
msg{1} = ['The classifier assigns the initial "iscell" labels to the ROIs.'];
msg{3} = ['Initially, it might not work very well, but will improve as you make choices in the GUI (and save the "proc" files).'];
msg{5} = ['To begin, press "new classifier?" and choose one of the provided priors.'];
msg{7} = ['The last used classifier will be automatically selected, every time you load a new plane.'];
msg{9} = ['As you process datasets, a database of your cells is built, and the classifier is re-trained.'];
msg{11} = ['The prior counts for about 300 cells, so your choices will start making a difference after a few hundred manually validated cells.'];
msg{13} = ['The system allows you to build personalized classifiers for different kinds of data (i.e. somas, dendrites, boutons, different brain areas, calcium indicators, or zoom levels).'];

msgbox(msg, 'Classifier instructions');

function pushbutton111_Callback(hObject, eventdata, handles)
msg{1} = ['Zoom in on portion of the image. Quadrants have 10% overlap.'];

msg{3} = ['Buttons become dark grey after visiting a quadrant.'];

msgbox(msg, 'ZOOM panel instructions');

function refresh_stats(handles,state)

nrois= nan;
ncells = nan;
if isfield(handles,'st0')
    nrois = size(handles.st0,1);
    ncells = handles.N_cells;
end
switch state
    case 1
        thestr = sprintf('%d Segments\n%d (%i%%) Cells\n\nReady.',nrois,ncells,round(100*ncells/nrois));
        set(handles.datastatsID,'BackGroundcolor','g')
    case 2
        thestr = sprintf('%d Segments\n%d (%i%%) Cells\n\nPROCESSING!!',nrois,ncells,round(100*ncells/nrois));
        set(handles.datastatsID,'BackGroundcolor','r')
end
set(handles.datastatsID,'String',thestr);
drawnow();

function edit51_Callback(hObject, eventdata, h)
% hObject    handle to edit51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit51 as text
%        str2double(get(hObject,'String')) returns contents of edit51 as a double
h.dat.cl.red_threshold = str2double(get(h.edit51,'String'));
for j = 1:length(h.dat.stat)
    h.dat.stat(j).isred = h.dat.stat(j).cellProb > h.dat.cl.red_threshold;
end
redraw_figure(h);

guidata(hObject,h);


% --- Executes during object creation, after setting all properties.
function edit51_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton112.
function pushbutton112_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton112 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.show_dff = -1;
guidata(hObject,handles);
redraw_fluorescence(handles);
refresh_figure(handles);

% --- Executes on button press in pushbutton113.
function pushbutton113_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton113 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(handles.zoom_handle.Enable,'on')
    handles.zoom_handle.Enable='off';
    set(handles.(sprintf('pushbutton%d', 113)), 'BackgroundColor',0.94*[1 1 1]);    
else
    handles.zoom_handle.Enable='on';
    set(handles.(sprintf('pushbutton%d', 113)), 'BackgroundColor',[1 0 0]);    
end       
redraw_fluorescence(handles);
guidata(hObject,handles);

% --- Executes on button press in pushbutton114.
function pushbutton114_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton114 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% hObject    handle to pushbutton114 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'is_ROI_making_mode')
    handles.is_ROI_making_mode = 0;
end

if strcmp(handles.zoom_handle.Enable,'on')
    handles.zoom_handle.Enable='off';
    set(handles.(sprintf('pushbutton%d', 114)), 'BackgroundColor',0.94*[1 1 1]);        
end

if handles.is_ROI_making_mode
    handles.is_ROI_making_mode = 0;
    set(handles.(sprintf('pushbutton%d', 114)), 'BackgroundColor',0.94*[1 1 1]);    
    try
       delete(handles.custom_ROI.plothandle);
    catch
        
    end
else
    handles.is_ROI_making_mode = 1;
    set(handles.(sprintf('pushbutton%d', 114)), 'BackgroundColor',[1 0 0]);    
end       

guidata(hObject, handles);

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton114.
function pushbutton114_ButtonDownFcn(hObject, eventdata, handles)




% --- Executes on scroll wheel click while the figure is in focus.
function figure1_WindowScrollWheelFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)
if handles.is_ROI_making_mode==2
    ratio = handles.dat.cl.Lx/handles.dat.cl.Ly;
    handles.custom_ROI.Radius = max(3,min(35,handles.custom_ROI.Radius - (eventdata.VerticalScrollAmount*eventdata.VerticalScrollCount/3)));
    r = handles.custom_ROI.Radius;
    try
       delete(handles.custom_ROI.plothandle);
    catch
        
    end    
    handles.custom_ROI.plothandle = rectangle('Position',[handles.custom_ROI.Center-[r,r/ratio],2*r,2*r/ratio],'Curvature',[1,1],'EdgeColor','g','LineWidth',3);   
    guidata(hObject, handles);
end

guidata(hObject, handles);


function h = getROISignal(h)

path = fileparts(h.dat.filename);
[~,file] = fileparts(h.dat.ops.RegFile);

if any(h.dat.map-[3,4]==0)
    
    file = [path,filesep,file,'_RED.bin'];
    if ~exist(file,'file')
        error('Registered RED binary nor found! Requesting file:\n  %s\n',file);
    end    
    isred = 1;
    
else
    
    file = [path,filesep,file,'.bin'];
    if ~exist(file,'file')
        error('Registered GREEN binary nor found! Requesting file:\n  %s\n',file);
    end
    isred = 0;
        
end

TEMP_ops = h.dat.ops;

TEMP_ops.RegFile = file;

a = -h.custom_ROI.Radius:h.custom_ROI.Radius;
[x,y] = meshgrid(h.custom_ROI.Center(1) + a,h.custom_ROI.Center(2) + a);

good = find(not(sqrt( (x-h.custom_ROI.Center(1)).^2 + (y-h.custom_ROI.Center(2)).^2 ) >  h.custom_ROI.Radius));

coord = [x(good),y(good)];

coord(:,1) = coord(:,1) + TEMP_ops.xrange(1);
coord(:,2) = coord(:,2) + TEMP_ops.yrange(1);
coord = coord - 1;

fprintf('\nExtracting mean signal of the selected ROI (this can take several minutes)...\n\n');

refresh_stats(h,2);
h.custom_ROI.sig = readBin(TEMP_ops,0,'all',h.dat.ops.iplane,isred,coord);

h.custom_ROI.sig_isred = isred;

h.custom_ROI.sig_sourcefile = h.dat.ops.RegFile;
refresh_stats(h,1);


% --- Executes on button press in pushbutton115.
function pushbutton115_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton115 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.add_segment_halo
    handles.add_segment_halo = 0;
    set(handles.(sprintf('pushbutton%d', 115)), 'BackgroundColor',0.94*[1 1 1]);    
else
    handles.add_segment_halo = 1;
    set(handles.(sprintf('pushbutton%d', 115)), 'BackgroundColor',[1 0 0]);    
end      
refresh_figure(handles);
guidata(hObject,handles);

function refresh_figure(handles)

if handles.dat.map==1
    redraw_figure(handles);
else
    redraw_meanimg(handles);
end
refresh;

