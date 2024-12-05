% ds000117 data
% Main script to execute different smoothing methods on 3D data
%
%--------------------------------------------------------------------------
% ess0001 F-contrast named Canonical HRF effects of interest, which tests
% if there is any significant activation across the defined regressors.
%
% represent the linear combinations of the parameter estimates (e.g., 
% contrasts like Faces>Scrambled Faces or conditions like Famous, Unfamiliar, etc.).
% to see the impact of smoothing on the statistical results derived from specific conditions or task contrasts. 
% The edges between activated and non-activated areas (often at the gray/white matter boundary) are highlighted in these images
% con_0002 represents the contrast for Faces > Scrambled Faces.
% con_0003 is the contrast for Famous.
% con_0004 corresponds to Unfamiliar.
% con_0005 represents Scrambled.
%
%--------------------------------------------------------------------------
% REFERENCES
% BIDS: Gorgolewski, K.J., et al. (2016). BIDS: The brain imaging data structure. A standard for organizing and describing outputs of neuroimaging experiments. Scientific Data, 3, 160044.
% ds000117: Wakeman, D.G. & Henson, R.N. (2015). A multi-subject, multi-modal human neuroimaging dataset. Sci. Data 2:150001 doi: 10.1038/sdata.2015.1
% ds000117: Henson, R.N., Wakeman, D.G., Litvak, V. & Friston, K.J. (2011). A Parametric Empirical Bayesian framework for the EEG/MEG inverse problem: generative models for multisubject and multimodal integration. Frontiers in Human Neuroscience, 5, 76, 1-16.
% ds000117: Chapter 42 of the SPM12 manual (http://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf)
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------
%% In Silico Data
% Cleaning environment & setting up SPM path
close all; clear; clc;
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');
[param, flag] = aj_smooth_default();

% Paths to access to the data
ds_dir = 'D:\Master_Thesis\Data\InSilicoData';
param.outDerivName = 'xxx';
BIDS_ph = fullfile(ds_dir, 'derivatives', param.outDerivName);

% Use BIDS to get data information
BIDS = spm_BIDS(ds_dir);
nsub = length(BIDS.subjects);

% Set up the MPM names list
MPMs_listname = {'MTsat', 'PDmap', 'R1map', 'R2starmap'};
nMPMnames = length(MPMs_listname);

% Looking for warped MPMs and TC segmentation maps for each subject
wMPM_paths = cell(nsub,1);
TCseg_paths = cell(nsub,1);


%% Public OpenNeuro dataset: ds000117 (fMRI)
% Cleaning environment & setting up SPM path
close all; clear; clc;
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');
[param, flag] = aj_smooth_default();

% Paths to access to the data
work_bids_root = 'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\work_copies\openneuro.org\ds000117';
param.outDerivName = '1stat';
stat_path = fullfile(work_bids_root, 'derivatives', param.outDerivName);
preproc_path = fullfile(work_bids_root, 'derivatives', 'preprocessing');

% Use BIDS to get data information
BIDS_stat = spm_BIDS(stat_path);
BIDS_preproc = spm_BIDS(preproc_path);
nsub = length(BIDS_stat.subjects); % stat is the most limiting for number of subjects

% Looking for constrat MPMs and modulated warped TC for each subject
conMPM_paths = cell(nsub,1);
mwTC_paths = cell(nsub,1);
for i = 1:nsub
    conMPM_paths{i} = spm_select('FPListRec', fullfile(BIDS_stat.subjects(i).path,'func'), '^.*con.*\.nii$');
    mwTC_paths{i} = spm_select('FPListRec', fullfile(BIDS_preproc.subjects(i).path, 'anat'), '^.*mwc.*\.nii$'); % well sorted
end

% Initialize cell arrays to store the smoothed results for all subjects
gs_imgaussfilt3_paths = cell(nsub,1);
gs_spm_paths = cell(nsub,1);
tws_paths = cell(nsub,1);
smwTC_paths = cell(nsub,1);
tspoon_paths = cell(nsub,1);

% Start parallel pool if not already open
if isempty(gcp('nocreate'))
    parpool; % Create a default parallel pool
end

% Call smoothing functions for each subject with loaded NIfTI data
parfor i = 1:nsub
    fprintf('Executing smoothing for subject %d...\n', i);
    
    [gs_imgaussfilt3_paths{i}, gs_spm_paths{i},...  % Gaussian results
        tws_paths{i}, smwTC_paths{i},...            % TWS results
        tspoon_paths{i}] = ...                      % TSPOON results
        aj_smoothing(conMPM_paths{i}, mwTC_paths{i}, param, flag, 3); % 3D so dim = 3
end

delete(gcp('nocreate'));

%% Published Callaghan dataset: AgingData (qMRI)
% REFERENCES
% M.F. Callaghan and al. (2014) http://dx.doi.org/10.1016/j.neurobiolaging.2014.02.008
% https://hackmd.io/u_vOEzA8TzS1yGj52V6Txg
% https://github.com/CyclotronResearchCentre/BIDS_AgingData

% Cleaning environment & setting up SPM path
close all; clear; clc;
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');
[param, flag] = aj_smooth_default();

% Paths to access to the data
ds_dir = 'D:\Master_Thesis\Data\BIDS_AgingData';
param.outDerivName = 'SPM12_dartel';
BIDS_warped_data = fullfile(ds_dir, 'derivatives', param.outDerivName);

% Use BIDS to get data information
BIDS_stat = spm_BIDS(BIDS_warped_data);
nsub = length(BIDS_stat.subjects);

% Set up the MPM names list
MPMs_listname = {'MTsat', 'PDmap', 'R1map', 'R2starmap'};
nMPMnames = length(MPMs_listname);

% Looking for warped MPMs and TC segmentation maps for each subject
wMPM_paths = cell(nsub,1);
TCseg_paths = cell(nsub,1);
parfor i = 1:nsub
    CSF_path_i = spm_select('FPListRec', fullfile(BIDS_stat.subjects(i).path,'anat'), '^*CSF_probseg*\.nii$');
    GM_path_i = spm_select('FPListRec', fullfile(BIDS_stat.subjects(i).path,'anat'), '^*GM_probseg*\.nii$');
    WM_path_i = spm_select('FPListRec', fullfile(BIDS_stat.subjects(i).path,'anat'), '^*WM_probseg*\.nii$');
    TCseg_paths{i} = char(GM_path_i, WM_path_i, CSF_path_i); % well sorted
    
    for ii = 1:nMPMnames
        % Get the file list for current subject and MPM
        file_list = spm_select('FPListRec', fullfile(BIDS_stat.subjects(i).path,'anat'), ['^*' MPMs_listname{ii} '*\.nii$']);
        if ~isempty(file_list)
            % Ensure wMPM_paths{i} is a character array
            if isempty(wMPM_paths{i})
                wMPM_paths{i} = file_list; % Initialize with first file list
            else
                % Concatenate character arrays vertically
                wMPM_paths{i} = char(wMPM_paths{i}, file_list); 
            end
        end
    end
end

% Initialize cell arrays to store the smoothed results for all subjects
gs_imgaussfilt3_paths = cell(nsub,1);
gs_spm_paths = cell(nsub,1);
tws_paths = cell(nsub,1);
smwTC_paths = cell(nsub,1);
tspoon_paths = cell(nsub,1);

% Start parallel pool if not already open
if isempty(gcp('nocreate'))
    parpool; % Create a default parallel pool
end

% Parallel loop for smoothing
parfor i = 51:138
    fprintf('Executing smoothing for subject %d...\n', i);
    
    [gs_imgaussfilt3_paths{i}, gs_spm_paths{i},...  % Gaussian results
        tws_paths{i}, smwTC_paths{i},...            % TWS results
        tspoon_paths{i}] = ...                      % TSPOON results
        aj_smoothing(wMPM_paths{i}, TCseg_paths{i}, param, flag, 3); % 3D so dim = 3
end

delete(gcp('nocreate'));
