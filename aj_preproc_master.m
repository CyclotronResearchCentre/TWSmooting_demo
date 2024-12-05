%--------------------------------------------------------------------------
% Main script to preprocess data from ds000117
%
% PROCESS
% This script automates the preprocessing of fMRI data and subsequent
% conversion of outputs to the BIDS format for an entire dataset ds000117.
% It uses parallel computing to speed up the processing across multiple 
% subjects. The preprocessing includes standard steps such as realignment, 
% coregistration, segmentation, normalization, and smoothing. After 
% preprocessing, the results are reformatted to comply with the BIDS 
% standard.
%
% REFERENCES
% SPM: Friston, K.J., et al. (1994). Statistical parametric maps in functional imaging: A general linear approach. Human Brain Mapping, 2, 189-210.
% matlabbatch{1}: Friston, K.J., et al. (1996). https://doi.org/10.1002/mrm.1910350312
% matlabbatch{2}: Ashburner, J. & al. (2005). https://doi.org/10.1016/j.neuroimage.2005.02.018
% matlabbatch{3}: Ashburner, J. & al. (1997). https://doi.org/10.1006/nimg.1997.0290
% matlabbatch{4}: Ashburner, J. & Friston, K.J. (1999). https://doi.org/10.1002/(sici)1097-0193(1999)7:4&#x0003c;254::aid-hbm4&#x0003e;3.0.co;2-g
% 
% BIDS: Gorgolewski, K.J. & al. (2016). https://doi.org/10.1038/sdata.2016.44
% ds000117: Wakeman, D.G. & al. (2015). https://doi.org/10.1038/sdata.2015.1
% ds000117: Henson, R.N., Wakeman, D.G., Litvak, V. & Friston, K.J. (2011). 
%           A Parametric Empirical Bayesian framework for the EEG/MEG inverse problem:
%           generative models for multisubject and multimodal integration. 
%           Frontiers in Human Neuroscience, 5, 76, 1-16.
% ds000117: Chapter 42 of the SPM12 manual (http://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf)
% Codes: Multimodal integration of M/EEG and f/MRI data in SPM12 of gllmflndn ; https://github.com/spm/MultimodalScripts/tree/master
%--------------------------------------------------------------------------
% FUTURE DEV
% use bidsme instead of aj_bidsify_batchresults
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------
%% Step 0: Reset environment + ask user 
close all; clear; clc;
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12'); % Add SPM12 to MATLAB path

% Defining subjects to process
proc_sub = 1:4;

%% Step 1: Define dataset paths and extract subjects
disp('Initializing dataset paths and subject extraction...');

base_dir = 'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117';
work_bids_root = fullfile(base_dir, 'work_copies', 'openneuro.org', 'ds000117');
orig_bids_root = fullfile(base_dir, 'downloads', 'openneuro.org', 'ds000117');

% Extract subjects from BIDS dataset
subs = spm_BIDS(work_bids_root, 'subjects'); 
nsub = numel(subs);

% Format subjects into appropriate string list
subjects = cell(length(proc_sub), 1);
for s = proc_sub
    subjects{s} = sprintf('sub-%s', subs{s});
end

% Restore work copies for all subjects
disp('Restoring work copies for all subjects...');
aj_restore_work_copies(subjects, orig_bids_root, work_bids_root);

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
parfor s = proc_sub
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
    end
    
    % Load the anatomical T1w images
    inputs{nrun+1} = cellstr(spm_select('ExtFPList', ...
        fullfile(work_bids_root, sprintf('sub-%02d', s), 'ses-mri', 'anat'), ...
        '^sub-.*_T1w\.nii$', Inf));
    inputs{nrun+2} = cellstr(spm_select('ExtFPList', ...
        fullfile(work_bids_root, sprintf('sub-%02d', s), 'ses-mri', 'anat'), ...
        '^sub-.*_T1w\.nii$', Inf));
    
    % Run the preprocessing job for the current subject
    spm_jobman('run', jobfile, inputs{:});
end

disp('PREPROCESSING: DONE');
disp('_____________________________________________________________________');

%% Step 3: Post-process and BIDSify the results for each subject
disp('Post-processing results and converting to BIDS format...');

for s = proc_sub
    subject = subjects{s};  % Get current subject
    disp(['BIDSifying results for subject: ' subject]);
    
    % Convert the results to BIDS format
    aj_bidsify_batchresults(orig_bids_root, work_bids_root, subject);
end

%% Step 4: Create SPM's multiple conditions files
% This part creates a text file that contains the 6 movement parameters for
% each scan, which was created during the realignment (preprocessing part).
% These parameters will be added to the GLM to capture residual 
% motion-related artifacts in the data

disp('Create SPMs multiple conditions files...');

derivpath = 'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\work_copies\openneuro.org\ds000117\derivatives';
trialtypes = {'Famous','Unfamiliar','Scrambled'}; % impose order

for s = proc_sub
    for r = 1:nrun  
            d = spm_load(char(spm_BIDS(work_bids_root,'data',...
                'modality','func','type','events','sub',subs{s},'run',runs{r}))); 
            clear conds 
        for t = 1:numel(trialtypes) 
                    conds.names{t} = trialtypes{t}; 
                    conds.durations{t} = 0; 
                    conds.onsets{t} = d.onset(strcmpi(d.stim_type,trialtypes{t}));  
        end 
            save(fullfile(derivpath,'preprocessing',...
                sprintf('sub-%s',subs{s}),'ses-mri','func',...
                sprintf('sub-%s_run-%s_spmdef.mat',subs{s},runs{r})),'-struct','conds'); 
    end 
end

disp('_____________________________________________________________________');

%% Step 5: First-level statistics
disp('Starting first-level statistics processing with parallel processing...');

% Initialize SPM configuration
spm_jobman('initcfg'); 
spm('defaults', 'fmri');

% Start parallel pool if not already open
numworkers = gcp('nocreate');
if isempty(numworkers)
    parpool; % Create a default parallel pool
    numworkers = gcp;
end

parfor s = proc_sub
    jobfile = {fullfile(base_dir,'aj_1stat_master_job.m')}; 
    inputs  = cell(nrun*3+1,1);
    inputs{1} = {fullfile(derivpath,'1stat',sprintf('sub-%02d',s),'ses-mri','func')};
    for r = 2:3:nrun*3+1
        inputs{r} = cellstr(spm_select('ExtFPList',...
            fullfile(derivpath,'preprocessing',sprintf('sub-%02d',s),'ses-mri','func'),...
            sprintf('^wsub-.*run-%02d_bold\\.nii', (r+1)/3), Inf));
        inputs{r+1} = cellstr(spm_select('FPList',...
            fullfile(derivpath,'preprocessing', sprintf('sub-%02d',s),'ses-mri','func'),...
            sprintf('^sub-%02d_run-%02d_spmdef\\.mat', s, (r+1)/3)));
        inputs{r+2} = cellstr(spm_select('FPList',...
            fullfile(derivpath,'preprocessing',sprintf('sub-%02d',s),'ses-mri','func'),...
            sprintf('^rp.*run-%02d.*\\.txt', (r+1)/3)));
    end
    spm_jobman('run', jobfile, inputs{:});
end

disp('_____________________________________________________________________');

%% Step 6: Close the parallel pool to free up resources
disp('Closing the parallel pool...');
delete(gcp('nocreate'));

disp('Processing complete. READY FOR SMOOTHING!!!');
