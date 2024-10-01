% Main script to generate 'n' phantoms and extract 1D, 2D, and 3D data

clear all;
close all;
clc;

addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');

% Step 0: Clean up previous run files
delete('phantom_3D_*.nii');
delete('phantom_1D_*.mat');
delete('phantom_2D_*.mat');
fprintf('Previous run files deleted.\n');

[model, param, flag] = aj_phantom_default();

% Get the default model type (e.g., 'modified_shepp_logan')
model_func = model.(param.model_type);

% Calculate resolution (mm per voxel)
resolution = param.FOV / param.grid_size;
fprintf('Resolution: %.2f mm per voxel\n', resolution);

ellipsoids = model_func();
ellipsoids(:, 2:7) = ellipsoids(:, 2:7) * resolution;
model_func = ellipsoids;

% Initialize containers for the results
phantoms_1D = cell(param.n, 1);
phantoms_2D = cell(param.n, 1);
phantoms_3D = cell(param.n, 1);
ellipses_all = cell(param.n, 1); % Initialize a cell array to store ellipses for each phantom

% Generate 'n' phantoms
for i = 1:param.n
    fprintf('Generating phantom %d/%d...\n', i, param.n);

    % Step 1: Create the original 3D phantom
    [phantom_3D, ellipse] = aj_create_phantom_3d(model_func, param.FOV); % param.FOV to get [mm] unit
    
    % Step 2: Add anatomical variability and noise
    [phantom_anatomical, phantom_noise, phantom_var, ellipse_anatomical] = ...
        aj_add_anatomical_variability(phantom_3D, ellipse, param.jitter_range, param.jitter_factor, param.prenoise_level, param.sm_kern, param.noise_range, param.noise_level, flag.noise_before_smoothing, flag.plot_fig);
    
    % Store the ellipses with anatomical variability
    ellipses_all{i} = ellipse_anatomical; % Save ellipses for the current phantom
    
    % Step 3: Extract 1D and 2D data from the noisy phantom
    [data_1D, data_2D] = aj_extract_data_from_phantom(phantom_var, flag.plot_fig);

    % Store the results
    phantoms_3D{i} = phantom_var;  % 3D phantom with anatomical variability and noise
    phantoms_1D{i} = data_1D;  % 1D profile
    phantoms_2D{i} = data_2D;  % 2D slice
    
    % Save 1D data as .mat file
    data_1D_filename = sprintf('phantom_1D_%d.mat', i);
    save(data_1D_filename, 'data_1D');
    fprintf('Saved: %s\n', data_1D_filename);
    
    % Save 2D data as .mat file
    data_2D_filename = sprintf('phantom_2D_%d.mat', i);
    save(data_2D_filename, 'data_2D');
    fprintf('Saved: %s\n', data_2D_filename);
    
    % Save the 3D phantom as a NIFTI file
    nifti_filename = sprintf('phantom_3D_%d.nii', i);
    V = struct(); % Create the NIFTI volume structure
    V.fname = nifti_filename;  % Filename
    V.dim = size(phantom_var);  % Dimensions of the volume
    V.dt = [spm_type('float32') 0];  % Data type
    V.mat = eye(4);  % Identity transformation matrix
    spm_write_vol(V, phantom_var); % Write the volume to file
    fprintf('Saved: %s\n', nifti_filename);
end

% Save ellipses_all as a .mat file
save('ellipses_all.mat', 'ellipses_all');
fprintf('Saved: ellipses_all.mat\n');

fprintf('Phantom generation completed.\n');

% Call statistical analysis if more than 2 phantoms are generated
if param.n >= 2
    aj_phantom_stat_analysis(phantoms_3D, param);
    aj_quantify_phantom_variability(phantoms_3D);
    aj_quantify_ellipse_variability(ellipses_all);
end
