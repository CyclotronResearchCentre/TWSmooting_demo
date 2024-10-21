% In silico data
% Main script to execute different smoothing methods on 1D, 2D and 3D data

%% Step 0: Cleaning environment & setting up SPM path
close all;
clear;
clc;

delete('results_phantom_1D_*.mat');
delete('results_phantom_2D_*.mat');
delete('results_phantom_3D_*.mat');

addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');

%% Step 1: Default parameter setup
[param, flag] = aj_smooth_default();

%% Step 2: Load data from ph_data_XD.mat files
current_path = pwd;
ph_files = dir(fullfile(current_path, 'ph_data_*.mat'));

if isempty(ph_files)
    error('No phantom file (ph_data_*.mat) found.');
else
    fprintf('Found %d phantom files.\n', length(ph_files));
end

%% 1D Smoothing
for i = 1:1%length(ph_files)
    % Load 1D data and tissue probabilities
    data_file = fullfile(current_path, ph_files(i).name); 
    fprintf('Loading data from %s...\n', data_file);
    load(data_file, 'ph_data');

    % Ensure the required fields exist
    if isfield(ph_data, 'ph_1D') && isfield(ph_data, 'noisy_proba_map_1D')
        ph_1D = ph_data.ph_1D;
        proba_1D = ph_data.noisy_proba_map_1D;
        dim = 1; % since 1D data was collected from the ph_data structure
    else
        error('1D data missing in file %s\n', ph_files(i).name);
    end
    
    % Ensure proba_1D has dimensions [nb_tissue x N]
    proba_1D = squeeze(proba_1D);
    nb_tissue = min(size(proba_1D));
    if size(proba_1D, 1) ~= nb_tissue
        proba_1D = proba_1D';
    end
    
    % Ensure ph_1D has dimensions [1 x N]
    if isrow(ph_1D)
        ph_1D = ph_1D';
    end
    
    % Call the smoothing functions for the 3 types
    [gsP_signal, twsP_signal, f_twsP_signal, tosP_signal, f_tosP_signal] = aj_smoothing(ph_1D, proba_1D, param, flag, dim);

    % Display results if the flag is set
    if flag.plot_fig
        aj_smoothing_plot(ph_files(i).name, gsP_signal, twsP_signal, f_twsP_signal, tosP_signal, f_tosP_signal, ph_1D, dim);
    end

    % Save the results into a specific .mat file
    if flag.save_data
        output_filename = sprintf('results_%s_1D.mat', ph_files(i).name(1:end-4));
        save(output_filename, 'gsP_signal', 'twsP_signal', 'tosP_signal');
        fprintf('Results saved in %s\n', output_filename);
    end
    
    % Reorganize signals per type rather than per subject
    [ggsP_GmWmCsf, ttwsP_signal, ttosP_signal] = aj_reorganize_signals_by_tissue(gsP_signal, twsP_signal, tosP_signal, dim, nb_tissue);
end

disp('Processing completed for all 1D files.');

%% 2D Processing
for i = 1:1%length(ph_files)
    % Load 2D data and tissue probabilities
    data_file = fullfile(current_path, ph_files(i).name); 
    fprintf('Loading data from %s...\n', data_file);
    load(data_file, 'ph_data');

    % Ensure the required fields exist
    if isfield(ph_data, 'ph_2D') && isfield(ph_data, 'noisy_proba_map_2D')
        ph_2D = ph_data.ph_2D; % [N x N] matrix
        proba_2D = ph_data.noisy_proba_map_2D; % [N x N x nb_tissue] matrix
        dim = 2; % since 2D data was collected from the ph_data structure
    else
        error('2D data missing in file %s\n', ph_files(i).name);
    end
    
    % Ensure proba_2D has dimensions [nb_tissue x N x N]
    dims = size(proba_2D);
    [nb_tissue, idx_min] = min(dims);
    permute_order = [idx_min, setdiff(1:3, idx_min)];
    if idx_min ~= 1
        proba_2D = permute(proba_2D, permute_order);
    end
    
    % Call the smoothing functions for the 3 types
    [gsP_signal, twsP_signal, f_twsP_signal, tosP_signal, f_tosP_signal] = ...
        aj_smoothing(ph_2D, proba_2D, param, flag, dim);
    
    % Display results if the flag is set
    if flag.plot_fig
        aj_smoothing_plot(ph_files(i).name, gsP_signal, twsP_signal, f_twsP_signal, tosP_signal, f_tosP_signal, ph_1D, dim);
    end

    % Save the results into a specific .mat file
    if flag.save_data
        output_filename = sprintf('results_2D_%s.mat', ph_files(i).name(1:end-4));
        save(output_filename, 'gsP_signal', 'twsP_signal', 'tosP_signal');
        fprintf('Results saved in %s\n', output_filename);
    end
    % Reorganize signals per type rather than per subject
    [ggsP_GmWmCsf, ttwsP_signal, ttosP_signal] = aj_reorganize_signals_by_tissue(gsP_signal, twsP_signal, tosP_signal, dim, nb_tissue);
end

disp('Processing completed for all 2D files.');

%% 3D Smoothing
clc; close all;
for i = 1:1%length(ph_files)
    % Load 3D data and tissue probabilities
    data_file = fullfile(current_path, ph_files(i).name); 
    fprintf('Loading data from %s...\n', data_file);
    load(data_file, 'ph_data');

    % Ensure the required fields exist
    if isfield(ph_data, 'ph_3D') && isfield(ph_data, 'noisy_proba_map_3D')
        ph_3D = ph_data.ph_3D; % [N x N x N] matrix
        proba_3D = ph_data.noisy_proba_map_3D; % [N x N x N x nb_tissue] matrix
        dim = 3; % since 3D data was collected from the ph_data structure
    else
        error('3D data missing in file %s\n', ph_files(i).name);
    end
    
    % Ensure proba_3D has dimensions [nb_tissue x N x N x N]
    dims = size(proba_3D);
    [nb_tissue, idx_min] = min(dims);
    permute_order = [idx_min, setdiff(1:4, idx_min)];
    if idx_min ~= 1
        proba_3D = permute(proba_3D, permute_order);
    end
    
    % Call the smoothing functions for the 3 types
    [gsP_signal, twsP_signal, f_twsP_signal, tosP_signal, f_tosP_signal] = ...
        aj_smoothing(ph_3D, proba_3D, param, flag, dim);

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
