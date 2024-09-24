% Script to generate 'n' phantoms and extract 1D, 2D, and 3D data

clear all;
close all;
clc;

% Ajouter le chemin vers SPM
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');

[model, param, flag] = aj_phantom_default();

% Calculate resolution (mm per voxel)
resolution = param.FOV / param.grid_size;
fprintf('Resolution: %.2f mm per voxel\n', resolution);

ellipsoids = model.modified_shepp_logan();
ellipsoids(:, 2:7) = ellipsoids(:, 2:7) * resolution;

% Initialize containers for the results
phantoms_1D = cell(param.n, 1);
phantoms_2D = cell(param.n, 1);
phantoms_3D = cell(param.n, 1);

% Generate 'n' phantoms
for i = 1:param.n
    fprintf('Generating phantom %d/%d...\n', i, param.n);

    % Step 1: Create the original 3D phantom
    [phantom_3D, ellipse] = aj_create_phantom_3d(ellipsoids, resolution * param.grid_size);
    
    % Step 2: Add anatomical variability and noise
    [phantom_anatomical, phantom_noise, phantom_var, ellipse_anatomical] = ...
        aj_add_anatomical_variability(phantom_3D, ellipse, param.jitter_range, param.jitter_factor, param.prenoise_level, param.sm_kern, param.noise_range, param.noise_level, flag.noise_before_smoothing, flag.plot_fig);
    
    % Step 3: Extract 1D and 2D data from the noisy phantom
    [data_1D, data_2D] = aj_extract_data_from_phantom(phantom_var, flag.plot_fig);

    % Store the results
    phantoms_3D{i} = phantom_var;  % 3D phantom with anatomical variability and noise
    phantoms_1D{i} = data_1D;  % 1D profile
    phantoms_2D{i} = data_2D;  % 2D slice
    
    % Save the 3D phantom as a NIFTI file
    nifti_filename = sprintf('phantom_3D_%d.nii', i);
    
    % Create the NIFTI volume structure
    V = struct();
    V.fname = nifti_filename;  % Filename
    V.dim = size(phantom_var);  % Dimensions of the volume
    V.dt = [spm_type('float32') 0];  % Data type
    V.mat = eye(4);  % Identity transformation matrix

    % Write the volume to file
    spm_write_vol(V, phantom_var);

    % Optionnel : Vérifier que le fichier a bien été créé
    fprintf('Saved: %s\n', nifti_filename);
end

fprintf('Phantom generation completed.\n');
