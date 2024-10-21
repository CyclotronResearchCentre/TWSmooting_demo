%% Description
% This script performs preprocessing steps on functional and anatomical 
% neuroimaging data using SPM. It allows the user to select a functional 
% image, validates the selection, retrieves corresponding anatomical 
% images, and performs various preprocessing steps such as realignment, 
% slice timing correction, coregistration, segmentation, normalization, 
% and smoothing.

%% Cleaning & SPM path
close all;
clear;
clc;

addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');

%% Get the paths to the NIfTI files of N subjects
N_func = 16;
N_anat = 16;
func_files = aj_get_bids_files(N_func, 'func');
anat_files = aj_get_bids_files(N_anat, 'anat');

%% File selection by the user
user_choice = input('Enter the number of the functional file you want to process: ');

% Validate the user input for functional files
if user_choice < 1 || user_choice > length(func_files)
    fprintf('Invalid selection. Please run the script again and select a valid number.\n');
    return;
end

% Select the corresponding functional file
selected_func_file = func_files{user_choice};
V_func = spm_vol(selected_func_file); % Load the header information for the selected functional volume

if length(V_func) > 1
    V_func = V_func(1); % Select the first volume
end

nslices = V_func.dim(3); % Get the number of slices from the functional volume
fprintf('Number of slices: %d\n', nslices);

%% Matching with anatomical files => pay attention to not take a FLASH anat nii
% Retrieve the anatomical file corresponding to the same subject
% It is assumed that anat and func files are ordered by subjects : FALSE
selected_anat_file = anat_files{user_choice};

% Check if the anatomical file exists
if isempty(selected_anat_file)
    fprintf('No anatomical image found for the selected functional image.\n');
    return;
end

V_anat = spm_vol(selected_anat_file); % Load the header information for the selected anatomical volume

fprintf('Selected anatomical file: %s\n', selected_anat_file);

%% Get the voxel size in mm
func_floor_voxel_size_mm = aj_voxel_size(V_func.mat);
anat_floor_voxel_size_mm = aj_voxel_size(V_anat.mat);

%% Get TR and slice timing parameters based on image metadata
tr_value = regexp(V_func.descrip, 'TR=(\d+)ms', 'tokens');

if ~isempty(tr_value)
    tr = str2double(tr_value{1}{1});  % Convert to double
    ta = tr - (tr / nslices);
    fprintf('TR value: %d ms\n', tr);
    fprintf('TA value: %d ms\n', ta);
else
    fprintf('TR not found in the description.\n');
end

%% Slice acquisition order from JSON documentation file (ds000117)
url = 'https://s3.amazonaws.com/openneuro.org/ds000117/task-facerecognition_bold.json?versionId=dyXAWn7z8ubaAkwlNbzAZ0EQu422yI70';
filename = 'task-facerecognition_bold.json';
websave(filename, url);

fid = fopen(filename);
raw = fread(fid, inf);
str = char(raw');  % Convert to string
fclose(fid);

jsonData = jsondecode(str); % Decode JSON content
sliceTiming = jsonData.SliceTiming; % Extract SliceTiming information
[sortedSliceTiming, originalIndices] = sort(sliceTiming);

disp('Original indices:');
disp(originalIndices);

%% Display
% spm_image('Display', selected_func_file);
% spm_image('Display', selected_anat_file);
 
%% Initialise SPM
spm('Defaults','fMRI');
spm_jobman('initcfg');

%% Preprocessing steps
clear matlabbatch

% Realign
matlabbatch{1}.spm.spatial.realign.estwrite.data{1} = cellstr(selected_func_file);

% Slice Timing Correction
matlabbatch{2}.spm.temporal.st.scans{1} = cellstr(spm_file(selected_func_file,'prefix','r'));
matlabbatch{2}.spm.temporal.st.nslices = nslices;
matlabbatch{2}.spm.temporal.st.tr = tr / 1000;  % Convert to seconds
matlabbatch{2}.spm.temporal.st.ta = ta / 1000;  % Convert to seconds
matlabbatch{2}.spm.temporal.st.so = originalIndices;
matlabbatch{2}.spm.temporal.st.refslice = round(nslices / 2);

% Coregister: Align functional to anatomical images
matlabbatch{3}.spm.spatial.coreg.estimate.ref    = cellstr(spm_file(selected_anat_file));
matlabbatch{3}.spm.spatial.coreg.estimate.source = cellstr(spm_file(selected_func_file,'prefix','mean'));

% Segment the anatomical image
matlabbatch{4}.spm.spatial.preproc.channel.vols  = cellstr(selected_anat_file);
matlabbatch{4}.spm.spatial.preproc.channel.write = [0 1];  % Save bias-corrected image
matlabbatch{4}.spm.spatial.preproc.warp.write    = [1 1];  % Save forward deformation fields

% Normalize the functional images using deformation fields from anatomical segmentation
matlabbatch{5}.spm.spatial.normalise.write.subj.def = cellstr(spm_file(selected_anat_file,'prefix','y_','ext','nii'));
matlabbatch{5}.spm.spatial.normalise.write.subj.resample = cellstr(spm_file(selected_func_file,'prefix','ar')); % Resample the slice-timed functional image
% matlabbatch{5}.spm.spatial.normalise.write.subj.resample = cellstr(char(spm_file(selected_func_file,'prefix','ar'), spm_file(selected_func_file,'prefix','mean')));
matlabbatch{5}.spm.spatial.normalise.write.woptions.vox = func_floor_voxel_size_mm;

% Normalize the anatomical image as well
matlabbatch{6}.spm.spatial.normalise.write.subj.def = cellstr(spm_file(selected_anat_file,'prefix','y_','ext','nii'));
matlabbatch{6}.spm.spatial.normalise.write.subj.resample = cellstr(spm_file(selected_anat_file,'prefix','m','ext','nii'));
matlabbatch{6}.spm.spatial.normalise.write.woptions.vox = anat_floor_voxel_size_mm;

% Smooth the functional images
% matlabbatch{7}.spm.spatial.smooth.data = cellstr(spm_file(selected_func_file,'prefix','war'));
% matlabbatch{7}.spm.spatial.smooth.fwhm = [8 8 8];  % Smoothing with a full-width half-maximum (FWHM) kernel of 8 mm

% Run the batch
spm_jobman('run', matlabbatch);

%% HELP FUNCTIONS
function floor_voxel_size_mm = aj_voxel_size(voxel_matrix)    
    % Extract voxel sizes
    voxel_mm = voxel_matrix(1:3, 1:3); % Get voxel sizes
    voxel_size_mm = sqrt(sum(voxel_mm.^2)); % Calculate voxel size
    
    % Display voxel sizes
    fprintf('Voxel size (in mm): [%.2f, %.2f, %.2f]\n', voxel_size_mm(1), voxel_size_mm(2), voxel_size_mm(3));
    
    % Floor the rounded voxel sizes
    floor_voxel_size_mm = floor(round(voxel_size_mm, 2)); % Round to avoid floating-point errors
    fprintf('Floor voxel size (in mm): [%.2f, %.2f, %.2f]\n', floor_voxel_size_mm(1), floor_voxel_size_mm(2), floor_voxel_size_mm(3));
end
