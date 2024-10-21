% Main script to:   generate n 3D phantoms
%                   extract 1D and 2D phantoms from the 3D ones
%                   compute the probability map in 1D, 2D and 3D
%                   compute the noisy probability map in 1D, 2D and 3D
%                   save each individual phantom structure and nifti files

%% Step 0: Cleaning stuff & SPM path
close all;
clear;
clc;

delete('ph_data_*.mat');
delete('ph_3D_*.nii');
delete('ph_proba_map_3D_*.nii');
delete('ph_noisy_proba_map_3D_*.nii');

addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');

%% Step 1: Preprocessing
[model, param, flag] = aj_phantom_default();

% Initialize storage for RMSE values
if param.n > 1 && flag.ph_stat
    rmse_values = zeros(param.n * (param.n - 1) / 2, 1);  % Initialize array to store RMSE values
    index = 1;  % Initialize index for storing RMSE values
end

%% Step 2: Generating 'n' phantoms
for ph_id = 1:param.n
    fprintf('Processing phantom %d...\r', ph_id);

    % Create the original 3D phantom
    fprintf('Generating the original phantom...');
    [phantom_3D, ellipse] = aj_phantom_create3D(model.(param.model_type)(), param);
    fprintf(repmat('\b', 1, length('Generating the original phantom...')));
    
    % Add anatomical variability and noise
    fprintf('Adding anatomical variability and noise...');
    [anat_ph_3D, noisy_ph_3D, final_ph_3D, ellipse_anatomical] = ...
        aj_phantom_variability(phantom_3D, ellipse, param, flag);
    fprintf(repmat('\b', 1, length('Adding anatomical variability and noise...')));
    
    % Extract 1D and 2D data from the 3D phantom
    fprintf('Extracting 1D and 2D phantom from the 3D...');
    [final_ph_1D, final_ph_2D] = aj_phantom_extract1D2D(final_ph_3D, param, flag.plot_fig);
    [anat_ph_1D, anat_ph_2D] = aj_phantom_extract1D2D(anat_ph_3D, param, false);
    fprintf(repmat('\b', 1, length('Extracting 1D and 2D phantom from the 3D...')));
    
    % Compute probability maps for each dimension
    fprintf('Generating probability maps...');
    proba_map_1D = aj_proba_maps(anat_ph_1D, param, flag);
    proba_map_2D = aj_proba_maps(anat_ph_2D, param, flag);
    proba_map_3D = aj_proba_maps(anat_ph_3D, param, flag);
    fprintf(repmat('\b', 1, length('Generating probability maps...')));

    % Step 2.5: Add noise to the probability maps
    fprintf('Adding noise to the probability maps...');
    noisy_proba_map_1D = aj_proba_noise(proba_map_1D, param.proba_noise);
    noisy_proba_map_2D = aj_proba_noise(proba_map_2D, param.proba_noise);
    noisy_proba_map_3D = aj_proba_noise(proba_map_3D, param.proba_noise);
    fprintf(repmat('\b', 1, length('Adding noise to the probability maps...')));
    
    % Step 2.6: Save results in a struct for each subject
    fprintf('Saving results...');
    ph_data.ph_1D = final_ph_1D;
    ph_data.ph_2D = final_ph_2D;
    ph_data.ph_3D = final_ph_3D;
    
    ph_data.anat_ph_1D = anat_ph_1D;
    ph_data.anat_ph_2D = anat_ph_2D;
    ph_data.anat_ph_3D = anat_ph_3D;
    
    ph_data.noisy_proba_map_1D = noisy_proba_map_1D;
    ph_data.noisy_proba_map_2D = noisy_proba_map_2D;
    ph_data.noisy_proba_map_3D = noisy_proba_map_3D;
    fprintf(repmat('\b', 1, length('Saving results...')));
    
    fprintf(['Phantom ', num2str(ph_id), ' processing complete.\n']);
    
    %% Step 2.7: Save the subject data to a separate .mat file
    filename = sprintf('ph_data_%d.mat', ph_id);
    save(filename, 'ph_data');
    disp(['Data saved to ', filename]);
    
    %% Step 2.8: Save the 3D phantom and (noisy) probability maps as a NIFTI file
    nifti_filename = sprintf('ph_3D_%d.nii', ph_id);
    V = struct();
    V.fname = nifti_filename;
    V.dim = size(final_ph_3D);
    V.dt = [spm_type('float32') 0];  % Data type
    V.mat = eye(4);  % Identity transformation matrix
    spm_write_vol(V, final_ph_3D); % Write the volume to file
    fprintf('Saved %s\n', nifti_filename);
    
    nifti_filename = sprintf('ph_proba_map_3D_%d.nii', ph_id);
    niftiwrite(proba_map_3D, nifti_filename);
    disp(['Saved ', nifti_filename]);
    
    noisy_nifti_filename = sprintf('ph_noisy_proba_map_3D_%d.nii', ph_id);
    niftiwrite(noisy_proba_map_3D, noisy_nifti_filename);
    disp(['Saved ', noisy_nifti_filename]);
    
    %% Step 2.9: Compute RMSE values if more than one phantom exists
    if param.n > 1 && ph_id > 1  && flag.ph_stat
        % Compute RMSE between the current phantom and all previous phantoms
        for j = 1:ph_id-1
            previous_data = load(sprintf('ph_data_%d.mat', j));
            rmse_values(index) = aj_phantom_RMSE(anat_ph_3D, previous_data.ph_data.anat_ph_3D);
            fprintf('RMSE between phantom %d and phantom %d: %.4f\n', ph_id, j, rmse_values(index));
            index = index + 1;
        end
    end
end

% After all phantoms have been processed
if param.n > 1 && flag.plot_ph_stat && flag.ph_stat
    % Optionally, calculate average RMSE
    avg_rmse = mean(rmse_values);
    fprintf('Average RMSE across all pairs: %.4f\n', avg_rmse);
    
    figure;
    histogram(rmse_values, 'Normalization', 'pdf');
    title('Histogram of RMSE Values');
    xlabel('RMSE');
    ylabel('Probability Density');

    figure;
    qqplot(rmse_values);
    title('Q-Q Plot of RMSE Values');
    
    mean_rmse = mean(rmse_values);
    median_rmse = median(rmse_values);
    std_rmse = std(rmse_values);
    fprintf('Mean RMSE: %.4f, Median RMSE: %.4f, Std Dev RMSE: %.4f\n', mean_rmse, median_rmse, std_rmse);
end

fprintf('Phantom generation completed.\n');
