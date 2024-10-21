%----------------------------------------------------------------------------------------------------
% Main script to preprocess fMRI data for a selected subject using SPM12
% This script automates the preprocessing steps including realignment, 
% segmentation, coregistration and normalization of functional and anatomical data 
% from a BIDS dataset.
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
% clear; close all; clc;

addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');

base_dir = 'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117';

% List of open inputs
nrun = 1; % enter the number of runs here
jobfile = {fullfile(base_dir, 'aj_batch_preproc_F_job.m')};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});

work_bids_root = fullfile(base_dir, 'work_copies', 'openneuro.org', 'ds000117');
orig_bids_root = fullfile(base_dir, 'downloads', 'openneuro.org', 'ds000117');
subject = 'sub-02 ';
session = 'ses-mri';

aj_bidsify_batchresults(orig_bids_root, work_bids_root, subject, session);
