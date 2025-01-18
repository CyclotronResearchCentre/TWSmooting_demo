%--------------------------------------------------------------------------
% Main script to:   generate n 3D phantoms
%                   extract 1D and 2D phantoms from the 3D ones
%                   compute the probability map in 1D, 2D and 3D
%                   compute the noisy probability map in 1D, 2D and 3D
%                   save each individual phantom structure and nifti files
%--------------------------------------------------------------------------
% FUTURE DEV
% extract 1d and 2d data also on proba maps -> only 3D ph runs
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------

% Cleaning environment & setting up SPM path
close all; clear; clc;
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');

ds_dir = 'D:\Master_Thesis\Data\InSilicoData';
[model, param, flag] = aj_phantom_default();

% Initialize storage for RMSE values
if param.n > 1 && flag.ph_stat
    rmse_values = zeros(param.n * (param.n - 1) / 2, 1);  % Initialize array to store RMSE values
    index = 1;  % Initialize index for storing RMSE values
end

% Custom phantom generation
custom_model = model.(param.model_type)();
ph_paths = cell(param.n,1);
NoisyProbaMap_paths = cell(param.n,1);
for i = 1:param.n
    fprintf('Processing phantom %d...\r', i);

    % Creating the original 3D phantom
    phantom_3D = aj_phantom_create3D(custom_model,param);
    
    % Adding anatomical variability and noise
    [anat_ph_3D, final_ph_3D] = aj_phantom_variability(phantom_3D, custom_model, param, flag);
    
    % Compute probability maps for each dimension
    proba_map_3D = aj_proba_maps(anat_ph_3D, param, flag);

    % Add noise to the probability maps
    noisy_proba_map_3D = aj_proba_noise(proba_map_3D, param.proba_noise);
    
    fprintf(['Phantom ', num2str(i), ' processing complete.\n']);
    
    % Save the subject data to nifti file
    pth_out = fullfile(ds_dir, sprintf('sub-%02d',i), 'anat'); % try to BIDSify
    if ~exist(pth_out,'dir')
        mkdir(pth_out)
    elseif flag.del_prevResults
        delete(fullfile(pth_out,'sub-*.mat'));
        delete(fullfile(pth_out,'sub-*.nii'));
    end
    
    fprintf('Results are saved in %s\n', pth_out);
    
    V = struct();
    V.fname = fullfile(pth_out,sprintf('sub-%02d_3Dph.nii', i));
    V.dim = size(final_ph_3D);
    V.dt = [spm_type('float32') 0];
    V.mat = eye(4);  % Identity transformation matrix
    V.descrip = 'Custom Phantom Generator: phantom';
    spm_write_vol(V, final_ph_3D);
    ph_paths{i} = V.fname;
    
    TC_names = {'GM','WM','CSF','BoneSculpt'};
    for ii = 1:size(noisy_proba_map_3D,4)
        NoisyProbaMap_paths{i} = fullfile(pth_out,sprintf('sub-%02d_%s_ProbaMap.nii', i, TC_names{ii}));
        niftiwrite(noisy_proba_map_3D(:,:,:,ii), NoisyProbaMap_paths{i}); 
    end
    
    % Extract 1D and 2D data from the 3D phantom
    if flag.data_1D2D
        [final_ph_1D, final_ph_2D] = aj_phantom_extract1D2D(final_ph_3D, param, flag.plot_fig);
        [proba_map_1D, proba_map_2D] = aj_phantom_extract1D2D(proba_map_3D, param, false);
        [noisy_proba_map_1D, noisy_proba_map_2D] = aj_phantom_extract1D2D(noisy_proba_map_3D, param, false);
    
        % Save results in a struct for each subject
        ph_data.ph_1D = final_ph_1D;
        ph_data.ph_2D = final_ph_2D;
        ph_data.anat_ph_1D = anat_ph_1D;
        ph_data.anat_ph_2D = anat_ph_2D;
        ph_data.noisy_proba_map_1D = noisy_proba_map_1D;
        ph_data.noisy_proba_map_2D = noisy_proba_map_2D;
        save(fullfile(pth_out,sprintf('sub-%02d_1D2Dph.mat', i)), 'ph_data');
    end
    
    % Compute RMSE values if more than one phantom exists
    if param.n > 1 && i > 1 && flag.ph_stat
        % Compute RMSE between the current phantom and all previous phantoms
        for j = 1:i-1
            previous_data = load(sprintf('ph_data_%d.mat', j));
            rmse_values(index) = aj_phantom_RMSE(anat_ph_3D, previous_data.ph_data.anat_ph_3D);
            fprintf('RMSE between phantom %d and phantom %d: %.4f\n', i, j, rmse_values(index));
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
