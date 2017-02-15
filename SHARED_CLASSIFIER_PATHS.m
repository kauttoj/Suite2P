% this is path for shared Suite2p classifier data files
% these files contain labeled training data to determine status of each segment (cell or not)

% this dataset is updated each time you save *proc.mat files
% copy of this dataset (before update) is also included in each *proc.mat file

% for now the selection of dataset is made by number of planes!!

CLASSIFIER_DATAFILE(1).file = 'Z:\commonKlabData\SUITE2P_classifier_files\standard_single_plane_WAVE.mat';
CLASSIFIER_DATAFILE(1).planes = 1;
CLASSIFIER_DATAFILE(1).bidirectional = 0;

CLASSIFIER_DATAFILE(2).file = 'Z:\commonKlabData\SUITE2P_classifier_files\six_plane_IARPA_WAVE.mat';
CLASSIFIER_DATAFILE(2).planes = 6;
CLASSIFIER_DATAFILE(2).bidirectional = 1;

