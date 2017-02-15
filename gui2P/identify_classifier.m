function [h]= identify_classifier(h)
% load classifier files
fprintf('\nLoading classifier data\n');

rootS2p = which('run_pipeline');
rootS2p_PATH=fileparts(rootS2p);
if isempty(rootS2p)
    warndlg('Suite2p location not in PATH! where is run_pipeline.m?')
    [filename1,rootS2p]=uigetfile(root, 'Select Data File');    
end
rootS2p = fileparts(rootS2p);
rootS2p = fullfile(rootS2p, 'configFiles');

flag = 0;
if isfield(h.dat.ops, 'classifier') && ~isempty(h.dat.ops.classifier)
    if ~exist(fullfile(rootS2p, h.dat.ops.classifier), 'file')
%         warndlg('specified classifier database not found, reverting  to last used');
    else
        def_file = dat.ops.classifier;
        flag = 1;
    end
else
%     warning('no specified classifier database, reverting  to last used');
end
h.is_shared_classifier = 0;
if (flag==0)
    try
        run([rootS2p_PATH,filesep,'SHARED_CLASSIFIER_PATHS.m']);
        if exist('CLASSIFIER_DATAFILE','var')
            
            isok = zeros(1,length(CLASSIFIER_DATAFILE));
            for i=1:length(CLASSIFIER_DATAFILE)
                if exist(CLASSIFIER_DATAFILE(i).file,'file')
                    isok(i)=1;
                end
            end
            fprintf('... found %i shared classification files:\n',sum(isok));
            if sum(isok)==0
                assert(0);
            end            
            sel=[];k=0;
            for i=1:length(CLASSIFIER_DATAFILE)
                if isok(i)
                    k=k+1;
                    fprintf('..... %i: %s\n',k,CLASSIFIER_DATAFILE(i).file);
                    if isempty(sel) && h.dat.ops.nplanes==CLASSIFIER_DATAFILE(i).planes
                        sel = [k,i];
                    end
                end
            end
            if ~isempty(sel)                
                def_file = CLASSIFIER_DATAFILE(sel(2)).file;
                h.dat.cl.fpath  = def_file;
                hload = load(h.dat.cl.fpath);
                if ~isfield(hload, 'st') || ~isfield(hload, 'statLabels') || ~isfield(hload, 'prior')
                    error('Classifier file %s is corrupt, should contain fields st, prior and statLabels!!', def_file);
                end
                fprintf('... classifier %i with %i planes is loaded (total %i samples)\n',sel(1),CLASSIFIER_DATAFILE(sel(2)).planes,size(hload.st,1));
                h.st        = hload.st;
                h.prior     = hload.prior;
                h.statLabels = hload.statLabels;
                h.is_shared_classifier = 1;
            else
                fprintf('... none of the classifiers is suitable (plane number mismatch)!\n');
                assert(0);
            end
        end
        
    catch
        
        warning(' No shared path classifier found, trying to open a local classifier!')
        
        fs = dir(fullfile(rootS2p, '*.mat'));
        if isempty(fs)
            warndlg('no classifier found, please make a new one!')
            h           = make_classifier(h);
        else
            [~, imax]    = max([fs.datenum]);
            def_file = fs(imax).name;
            h.dat.cl.fpath  = fullfile(rootS2p, def_file);
            hload = load(h.dat.cl.fpath);
            if ~isfield(hload, 'st') || ~isfield(hload, 'statLabels') || ~isfield(hload, 'prior')
                error('found a non-classifier file in configFiles, called %s. \n Please remove and try again!', def_file)
            end
            h.st        = hload.st;
            h.prior     = hload.prior;
            h.statLabels = hload.statLabels;
            fprintf('... found local classification file %s\n',h.dat.cl.fpath);
        end
    end
end

fprintf(' Classifier is ready to go!\n\n');

end

function h = make_classifier(h)
% new classifier button
rootS2p = which('run_pipeline');
if isempty(rootS2p)
    warndlg('Suite2p location not in PATH! where is run_pipeline.m?')
    [filename1,rootS2p]=uigetfile(root, 'Select run_pipeline.m');
    
end
rootS2p = fileparts(rootS2p);
rootS2p = fullfile(rootS2p, 'configFiles');

def_name = fullfile(rootS2p, 'cl_new.mat');
[FileName,PathName] = uiputfile('*.mat', 'Create new classifier', def_name); 

if FileName
    load(fullfile(rootS2p, 'priors', 'prior_default.mat'));
    st = [];
    save(fullfile(PathName, FileName), 'st', 'prior', 'statLabels')
    
    h.st = st;
    h.prior = prior;
    h.statLabels = statLabels;
    
    h.dat.cl.fpath          = fullfile(PathName, FileName);
    h                       = classROI(h);
    
%     h = splitROIleftright(h);
%     h = buildLambdaValue(h);
%     redraw_figure(h);
    
    hload = load(h.dat.cl.fpath);
    h.st        = hload.st;
    h.statLabels = hload.statLabels;

    h.prior     = hload.prior;
else
    error('you must make a new classifier first!')
end

end
