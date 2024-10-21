% ds000117 data
% Main script to execute different smoothing methods on 3D data

%% Cleaning environment & setting up SPM path
close all;
clear;
clc;

addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');

%% Find nii data from BIDS derivatives folder
% data_path = fullfile(pwd, 'data', 'ds000117');
% wfunc_paths = spm_select('FPList', data_path, '^wsub-01_ses-mri_task-facerecognition_run-.*\.nii$');
% wmean_paths = spm_select('FPList', data_path, '^wmeansub-01_ses-mri_task-facerecognition_run-.*\.nii$');
% mpm_paths = spm_select('FPList', data_path, '^mwc.*sub-01_ses-mri_acq-mprage_T1w\.nii$');

% Define the root directory of the BIDS dataset and derivatives
bids_root = 'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\work_copies\openneuro.org\ds000117';
derivatives_path = fullfile(bids_root, 'derivatives', 'preprocessing');

% Load the BIDS structure using spm_BIDS
bids_info = spm_BIDS(bids_root);

% Define the subject and session you are working on
subject = 'sub-03';
session = 'ses-mri';

% Fetch the files from the derivatives folder using aj_BIDS_select
wfunc_paths = aj_BIDS_select(bids_info, 'sub', subject, 'ses', session, ...
    'modality', 'func', 'pattern', 'wsub-*.nii', 'derivatives', derivatives_path);
wmean_paths = aj_BIDS_select(bids_info, 'sub', subject, 'ses', session, ...
    'modality', 'func', 'pattern', 'wmeansub-*.nii', 'derivatives', derivatives_path);
mpm_paths = aj_BIDS_select(bids_info, 'sub', subject, 'ses', session, ...
    'modality', 'anat', 'pattern', 'mwc*.nii', 'derivatives', derivatives_path);

%% Load the 3D data for functional images and anatomical masks
[param, flag] = aj_smooth_default();

for i = 1:numel(wfunc_paths)
    % Load the functional data
    wfunc_vol = spm_vol(wfunc_paths{i}); % Get volume information
    wfunc_3D = spm_read_vols(wfunc_vol(1)); % Read volume data

    % Load the mean functional data (if needed)
    wmean_vol = spm_vol(wmean_paths{i});
    wmean_3D = spm_read_vols(wmean_vol);
    
    % Load the anatomical probability map (MPM)
    mpm_vol = spm_vol(mpm_paths{i}); 
    mpm_3D = spm_read_vols(mpm_vol);

    % Perform the Gaussian smoothing
    dims = size(mpm_3D);
    [nb_tissue, idx_min] = min(dims);
    permute_order = [idx_min, setdiff(1:4, idx_min)];
    if idx_min ~= 1
        mpm_3D = permute(mpm_3D, permute_order);
    end

    dim = 3; % Data are in 3D

    % Call smoothing functions with loaded NIfTI data
    [gsP_signal, twsP_signal, f_twsP_signal, tosP_signal, f_tosP_signal] = ...
        aj_smoothing(wfunc_3D, mpm_3D, param, flag, dim);

    % Display results if the flag is activated
    if flag.plot_fig
        aj_smoothing_plot(wfunc_paths{i}, gsP_signal, twsP_signal, f_twsP_signal, tosP_signal, f_tosP_signal, [], dim);
    end

    % Sauvegarder les résultats dans un fichier .mat
    if flag.save_data
        output_filename = sprintf('results_3D_%s.mat', wfunc_paths(i,end-4:end));  % Générer un nom de fichier
        save(output_filename, 'gsP_signal', 'twsP_signal', 'tosP_signal');
        fprintf('Results saved in %s\n', output_filename);
    end

    % Réorganiser les signaux par type de tissu
    [ggsP_GmWmCsf, ttwsP_signal, ttosP_signal] = aj_reorganize_signals_by_tissue(gsP_signal, twsP_signal, tosP_signal, dim, nb_tissue);
end

disp('Processing completed for all 3D files.');



%% BACK UP 3D Smoothing
for i = 1:1%length(ph_files)
    % Load 3D data and tissue probabilities
    data_file = fullfile(current_path, ph_files(i).name); 
    fprintf('Loading data from %s...\n', data_file);
    load(data_file, 'ph_data');

    % Ensure the required fields exist
    if isfield(ph_data, 'ph_3D') && isfield(ph_data, 'noisy_proba_map_3D')
        wfunc_3D = ph_data.ph_3D; % [N x N x N] matrix
        mpm_3D = ph_data.noisy_proba_map_3D; % [N x N x N x nb_tissue] matrix
        dim = 3; % since 3D data was collected from the ph_data structure
    else
        error('3D data missing in file %s\n', ph_files(i).name);
    end
    
    % Ensure proba_3D has dimensions [nb_tissue x N x N x N]
    dims = size(mpm_3D);
    [nb_tissue, idx_min] = min(dims);
    permute_order = [idx_min, setdiff(1:4, idx_min)];
    if idx_min ~= 1
        mpm_3D = permute(mpm_3D, permute_order);
    end
    
    % Call the smoothing functions for the 3 types
    [gsP_signal, twsP_signal, f_twsP_signal, tosP_signal, f_tosP_signal] = ...
        aj_smoothing(wfunc_3D, mpm_3D, param, flag, dim);

    % Display results if the flag is set
    if flag.plot_fig
        aj_smoothing_plot(ph_files(i).name, gsP_signal, twsP_signal, f_twsP_signal, tosP_signal, f_tosP_signal, ph_1D, dim);
    end

    % Save the results into a specific .mat file
    if flag.save_data
        output_filename = sprintf('results_3D_%s.mat', ph_files(i).name(1:end-4));
        save(output_filename, 'gsP_signal', 'twsP_signal', 'tosP_signal');
        fprintf('Results saved in %s\n', output_filename);
    end
    
    % Reorganize signals per type rather than per subject
    [ggsP_GmWmCsf, ttwsP_signal, ttosP_signal] = aj_reorganize_signals_by_tissue(gsP_signal, twsP_signal, tosP_signal, dim, nb_tissue);
end

disp('Processing completed for all 3D files.');