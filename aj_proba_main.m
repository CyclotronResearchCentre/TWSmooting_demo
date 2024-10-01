% Main script to run 1D, 2D, and 3D tissue probability calculations

clear all;
close all;
clc;

addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');

[param, flag] = aj_proba_default();

current_path = pwd;  % Get the current working directory

% Step 1: Get all available 1D, 2D, and 3D files
files_1D = dir(fullfile(current_path, 'phantom_1D_*.mat'));
files_2D = dir(fullfile(current_path, 'phantom_2D_*.mat'));
files_3D = dir(fullfile(current_path, 'phantom_3D_*.nii'));

% Check if there are any files
if isempty(files_1D)
    fprintf('No 1D files found.\n');
else
    fprintf('Found %d 1D files.\n', length(files_1D));
end

if isempty(files_2D)
    fprintf('No 2D files found.\n');
else
    fprintf('Found %d 2D files.\n', length(files_2D));
end

if isempty(files_3D)
    fprintf('No 3D files found.\n');
else
    fprintf('Found %d 3D files.\n', length(files_3D));
end

% Step 2: Process each 1D, 2D, and 3D file
for i = 1:min([length(files_1D), length(files_2D), length(files_3D)])
    % Load 1D data
    data_1D_file = fullfile(current_path, files_1D(i).name);
    fprintf('Loading 1D data from %s...\n', data_1D_file);
    load(data_1D_file, 'data_1D');
    fprintf('1D data loaded successfully.\n');

    % Load 2D data
    data_2D_file = fullfile(current_path, files_2D(i).name);
    fprintf('Loading 2D data from %s...\n', data_2D_file);
    load(data_2D_file, 'data_2D');
    fprintf('2D data loaded successfully.\n');

    % Load 3D data using SPM
    data_3D_file = fullfile(current_path, files_3D(i).name);
    fprintf('Loading 3D data from %s...\n', data_3D_file);
    V = spm_vol(data_3D_file);  % Get volume information from the NIFTI file
    data_3D = spm_read_vols(V);  % Read the volume data
    fprintf('3D data loaded successfully.\n');

    % Call 1D, 2D & 3D probability calculation
    fprintf('Running 1D tissue probability calculation for phantom %d...\n', i);
    data_GmWmCsfSculpt_1D = aj_create_proba_1d(data_1D, param, flag.plot_fig);
    fprintf('Running 2D tissue probability calculation for phantom %d...\n', i);
    data_GmWmCsfSculpt_2D = aj_create_proba_2d(data_2D, param, flag.plot_fig);
    fprintf('Running 3D tissue probability calculation for phantom %d...\n', i);
    data_GmWmCsfSculpt_3D = aj_create_proba_3d(data_3D, param, flag.plot_fig);
    
    % Save 1D data as .mat file
    data_1D_filename = sprintf('proba_1D_%d.mat', i);
    save(data_1D_filename, 'data_GmWmCsfSculpt_1D');
    fprintf('Saved: %s\n', data_1D_filename);
    
    % Save 2D data as .mat file
    data_2D_filename = sprintf('proba_2D_%d.mat', i);
    save(data_2D_filename, 'data_GmWmCsfSculpt_2D');
    fprintf('Saved: %s\n', data_2D_filename);
    
    % Save 2D data as .mat file
    data_3D_filename = sprintf('proba_3D_%d.mat', i);
    save(data_3D_filename, 'data_GmWmCsfSculpt_3D');
    fprintf('Saved: %s\n', data_3D_filename);
end

fprintf('Processing complete for all available 1D, 2D, and 3D files.\n');


% Main script to run 1D, 2D, and 3D tissue probability calculations

% clear all;
% close all;
% clc;
% 
% addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');
% 
% [param, flag] = aj_proba_default();
% 
% current_path = pwd;  % Get the current working directory
% 
% % Step 1: Ask the user whether to process all files or specific ones
% choice = input('Do you want to process all files (type "all") or specific files (type "specific")? ', 's');
% 
% % Initialize file arrays
% files_1D = dir(fullfile(current_path, 'phantom_1D_*.mat'));
% files_2D = dir(fullfile(current_path, 'phantom_2D_*.mat'));
% files_3D = dir(fullfile(current_path, 'phantom_3D_*.nii'));
% 
% if strcmp(choice, 'all')
%     % Process all files
%     fprintf('Processing all found files...\n');
% 
%     % Check if there are any files
%     if isempty(files_1D)
%         fprintf('No 1D files found.\n');
%     else
%         fprintf('Found %d 1D files.\n', length(files_1D));
%     end
% 
%     if isempty(files_2D)
%         fprintf('No 2D files found.\n');
%     else
%         fprintf('Found %d 2D files.\n', length(files_2D));
%     end
% 
%     if isempty(files_3D)
%         fprintf('No 3D files found.\n');
%     else
%         fprintf('Found %d 3D files.\n', length(files_3D));
%     end
% 
%     % Process each file
%     for i = 1:min([length(files_1D), length(files_2D), length(files_3D)])
%         % Load 1D data
%         data_1D_file = fullfile(current_path, files_1D(i).name);
%         fprintf('Loading 1D data from %s...\n', data_1D_file);
%         load(data_1D_file, 'data_1D');  % Load the 'data_1D' variable
%         fprintf('1D data loaded successfully.\n');
% 
%         % Load 2D data
%         data_2D_file = fullfile(current_path, files_2D(i).name);
%         fprintf('Loading 2D data from %s...\n', data_2D_file);
%         load(data_2D_file, 'data_2D');  % Load the 'data_2D' variable
%         fprintf('2D data loaded successfully.\n');
% 
%         % Load 3D data using SPM
%         data_3D_file = fullfile(current_path, files_3D(i).name);
%         fprintf('Loading 3D data from %s...\n', data_3D_file);
%         V = spm_vol(data_3D_file);  % Get volume information from the NIFTI file
%         data_3D = spm_read_vols(V);  % Read the volume data
%         fprintf('3D data loaded successfully.\n');
% 
%         % Call 1D probability calculation
%         fprintf('Running 1D tissue probability calculation for phantom %d...\n', i);
%         data_GmWmCsfSculpt_1D = aj_create_proba_1d(data_1D, param, flag.plot_fig);
% 
%         % Call 2D probability calculation
%         fprintf('Running 2D tissue probability calculation for phantom %d...\n', i);
%         data_GmWmCsfSculpt_2D = aj_create_proba_2d(data_2D, param, flag.plot_fig);
% 
%         % Call 3D probability calculation
%         fprintf('Running 3D tissue probability calculation for phantom %d...\n', i);
%         data_GmWmCsfSculpt_3D = aj_create_proba_3d(data_3D, param, flag.plot_fig);
%         
%     end
% 
% elseif strcmp(choice, 'specific')
%     % Process specific files
%     fprintf('You chose to process specific files.\n');
% 
%     % Ask the user to input the indices of the files to process
%     indices = input('Please enter the indices of the files you want to process (e.g., [1, 2, 3]): ');
% 
%     % Validate indices
%     if max(indices) > min([length(files_1D), length(files_2D), length(files_3D)])
%         error('One or more indices are out of range.');
%     end
% 
%     % Process the specified files
%     for i = indices
%         % Load 1D data
%         data_1D_file = fullfile(current_path, files_1D(i).name);
%         fprintf('Loading 1D data from %s...\n', data_1D_file);
%         load(data_1D_file, 'data_1D');  % Load the 'data_1D' variable
%         fprintf('1D data loaded successfully.\n');
% 
%         % Load 2D data
%         data_2D_file = fullfile(current_path, files_2D(i).name);
%         fprintf('Loading 2D data from %s...\n', data_2D_file);
%         load(data_2D_file, 'data_2D');  % Load the 'data_2D' variable
%         fprintf('2D data loaded successfully.\n');
% 
%         % Load 3D data using SPM
%         data_3D_file = fullfile(current_path, files_3D(i).name);
%         fprintf('Loading 3D data from %s...\n', data_3D_file);
%         V = spm_vol(data_3D_file);  % Get volume information from the NIFTI file
%         data_3D = spm_read_vols(V);  % Read the volume data
%         fprintf('3D data loaded successfully.\n');
% 
%         % Call 1D probability calculation
%         fprintf('Running 1D tissue probability calculation for phantom %d...\n', i);
%         data_GmWmCsfSculpt_1D = aj_create_proba_1d(data_1D, param, flag.plot_fig);
% 
%         % Call 2D probability calculation
%         fprintf('Running 2D tissue probability calculation for phantom %d...\n', i);
%         data_GmWmCsfSculpt_2D = aj_create_proba_2d(data_2D, param, flag.plot_fig);
% 
%         % Call 3D probability calculation
%         fprintf('Running 3D tissue probability calculation for phantom %d...\n', i);
%         data_GmWmCsfSculpt_3D = aj_create_proba_3d(data_3D, param, flag.plot_fig);
%     end
% else
%     error('Invalid choice. Please type "all" or "specific".');
% end
% 
% fprintf('Processing complete.\n');
