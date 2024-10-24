%----------------------------------------------------------------------------------------------------
% This script automates the preprocessing of fMRI data and subsequent
% conversion of outputs to the BIDS format for an entire dataset (ds000117). It uses
% parallel computing to speed up the processing across multiple subjects.
% The preprocessing includes standard steps such as realignment, coregistration,
% segmentation, normalization, and smoothing. After preprocessing, the results
% are reformatted to comply with the BIDS standard.
%----------------------------------------------------------------------------------------------------
% References:
% SPM: Friston, K.J., et al. (1994). Statistical parametric maps in functional imaging: A general linear approach. Human Brain Mapping, 2, 189-210.
% matlabbatch{1}: Friston, K.J., et al. (1996). Movement-related effects in fMRI time-series. Magnetic Resonance in Medicine, 35, 346-355.
% matlabbatch{2}: Ashburner, J. & Friston, K.J. (2005). Unified segmentation. NeuroImage, 26(3), 839-851.
% matlabbatch{3}: Ashburner, J. & Friston, K.J. (1997). Multimodal image coregistration and partitioning—a unified framework. NeuroImage, 6(3), 209-217.
% matlabbatch{4}: Ashburner, J. & Friston, K.J. (1999). Nonlinear spatial normalization using basis functions. Human Brain Mapping, 7(4), 254-266.
% 
% BIDS: Gorgolewski, K.J., et al. (2016). BIDS: The brain imaging data structure. A standard for organizing and describing outputs of neuroimaging experiments. Scientific Data, 3, 160044.
% ds000117: Wakeman, D.G. & Henson, R.N. (2015). A multi-subject, multi-modal human neuroimaging dataset. Sci. Data 2:150001 doi: 10.1038/sdata.2015.1
% ds000117: Henson, R.N., Wakeman, D.G., Litvak, V. & Friston, K.J. (2011). A Parametric Empirical Bayesian framework for the EEG/MEG inverse problem: generative models for multisubject and multimodal integration. Frontiers in Human Neuroscience, 5, 76, 1-16.
% ds000117: Chapter 42 of the SPM12 manual (http://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf)
% Codes: Multimodal integration of M/EEG and f/MRI data in SPM12 of gllmflndn ; https://github.com/spm/MultimodalScripts/tree/master
%----------------------------------------------------------------------------------------------------

close all; clear; clc;
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12'); % Add SPM12 to MATLAB path

%% Step 1: Define dataset paths and extract subjects
disp('Initializing dataset paths and subject extraction...');

base_dir = 'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117';
work_bids_root = fullfile(base_dir, 'work_copies', 'openneuro.org', 'ds000117');
orig_bids_root = fullfile(base_dir, 'downloads', 'openneuro.org', 'ds000117');
session = 'ses-mri';

% Extract subjects from BIDS dataset
subs = spm_BIDS(work_bids_root, 'subjects'); 
nsub = numel(subs);

% Format subjects into appropriate string list
subjects = cell(nsub, 1);
for s = 3:6 %1:nsub
    subjects{s} = sprintf('sub-%s', subs{s});
end

% Restore work copies for all subjects
disp('Restoring work copies for all subjects...');
aj_restore_work_copies(subjects);

%% Step 2: Preprocess fMRI data for each subject using parallel computing
disp('Starting fMRI preprocessing with parallel processing...');

% Extract run information (number of runs for the task)
runs = spm_BIDS(work_bids_root,'runs', 'modality','func', 'type','bold', 'task','facerecognition'); 
nrun = numel(runs);

% Initialize SPM configuration
spm_jobman('initcfg'); 
spm('defaults', 'fmri');

% Start parallel pool if not already open
numworkers = gcp('nocreate');
if isempty(numworkers)
    parpool; % Create a default parallel pool
    numworkers = gcp;
end

% Preprocess each subject in parallel
% parfor (s = 1:nsub, numworkers)
parfor s = 3:6 %1:nsub
    % Display current subject being processed
    disp(['Processing subject: ' sprintf('sub-%02d', s)]);

    % Define the job file for SPM preprocessing
    jobfile = {fullfile(base_dir, 'aj_preproc_master_job.m')};
    inputs  = cell(nrun + 2, 1); 
    
    % Loop over each run and prepare the scan inputs
    for r = 1:nrun
        inputs{r} = cellstr(spm_select('ExtFPList', ...
            fullfile(work_bids_root, sprintf('sub-%02d', s), 'ses-mri', 'func'), ...
            sprintf('^sub-.*run-%02d_bold\\.nii', r), Inf));
        disp(['Run ' num2str(r) ' inputs: ' inputs{r}{1}]);
    end
    
    % Load the anatomical T1w images
    inputs{nrun+1} = cellstr(spm_select('ExtFPList', ...
        fullfile(work_bids_root, sprintf('sub-%02d', s), 'ses-mri', 'anat'), ...
        '^sub-.*_T1w\.nii$', Inf));
    inputs{nrun+2} = cellstr(spm_select('ExtFPList', ...
        fullfile(work_bids_root, sprintf('sub-%02d', s), 'ses-mri', 'anat'), ...
        '^sub-.*_T1w\.nii$', Inf));
    
    % Display anatomical inputs
    disp(['Anatomical inputs: ' inputs{nrun+1}{1}, ', ', inputs{nrun+2}{1}]);
    
    % Run the preprocessing job for the current subject
    spm_jobman('run', jobfile, inputs{:});
end

% Close the parallel pool to free up resources
disp('Closing the parallel pool...');
delete(gcp('nocreate'));

%% Step 3: Post-process and BIDSify the results for each subject
disp('Post-processing results and converting to BIDS format...');

for s = 3:6 %1:nsub
    subject = subjects{s};  % Get current subject
    disp(['BIDSifying results for subject: ' subject]);
    
    % Convert the results to BIDS format
    aj_bidsify_batchresults(orig_bids_root, work_bids_root, subject, session);
end

disp('Processing complete.');
