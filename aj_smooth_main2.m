%% Step 0: Cleaning
clear all;
close all;
clc;

%% Step 1: Clean up previous run files
delete('results_phantom_1D_*.mat');
delete('results_phantom_2D_*.mat');
delete('results_phantom_3D_*.mat');
fprintf('Previous run files deleted.\n');

%% Step 2: Parameters Setting
[param, flag] = aj_smooth_default();

%% Step 3: Get all 1D, 2D and 3D available files
current_path = pwd;
ph_files_1D = dir(fullfile(current_path, 'phantom_1D_*.mat'));
proba_files_1D = dir(fullfile(current_path, 'proba_1D_*.mat'));
ph_files_2D = dir(fullfile(current_path, 'phantom_2D_*.mat'));
proba_files_2D = dir(fullfile(current_path, 'proba_2D_*.mat'));
ph_files_3D = dir(fullfile(current_path, 'phantom_3D_*.nii'));
proba_files_3D = dir(fullfile(current_path, 'proba_3D_*.mat'));

if isempty(ph_files_1D)
    error('No 1D file found.\n');
elseif isempty(ph_files_2D)
    error('No 2D file found.\n');
elseif isempty(ph_files_3D)
    error('No 3D file found.\n');
else
    fprintf('Found %d 1D files.\n', length(ph_files_1D));
    fprintf('Found %d 2D files.\n', length(ph_files_2D));
    fprintf('Found %d 3D files.\n', length(ph_files_3D));
end

if length(ph_files_1D) ~= length(proba_files_1D)
    error('Le nombre de fichiers de données 1D et de fichiers de probabilités 1D ne correspond pas.');
elseif length(ph_files_2D) ~= length(proba_files_2D)
    error('Le nombre de fichiers de données 2D et de fichiers de probabilités 2D ne correspond pas.');
elseif length(ph_files_3D) ~= length(proba_files_3D)
    error('Le nombre de fichiers de données 3D et de fichiers de probabilités 3D ne correspond pas.');
end

%% 1D Smoothing
disp('Début du traitement des fichiers 1D...');
if ~isempty(ph_files_1D)
    h = waitbar(0, 'Traitement des fichiers 1D...');
    total_files = length(ph_files_1D);
    
    for i = 1:total_files
        data_1D_file = fullfile(current_path, ph_files_1D(i).name); 
        proba_1D_file = fullfile(current_path, proba_files_1D(i).name);
        
        fprintf('Chargement des données 1D depuis %s...\n', data_1D_file);
        load(data_1D_file, 'data_1D');

        fprintf('Chargement des probabilités tissulaires 1D depuis %s...\n', proba_1D_file);
        load(proba_1D_file, 'data_GmWmCsfSculpt_1D');

        % Exemple de traitement (à adapter selon tes besoins de lissage)
        % Appliquer des opérations sur les données ici...
        % Exemple : gsP_signal = smooth(data_1D);
        % Exemple : twsP_signal = smooth(data_1D, 'lowess');
        % Exemple : tosP_signal = some_other_smoothing_function(data_1D);

        if flag.save_data
            output_filename = sprintf('results_1D_%s.mat', ph_files_1D(i).name(1:end-4));
            save(output_filename, 'gsP_signal', 'twsP_signal', 'tosP_signal');
        end
        
        waitbar(i / total_files, h, sprintf('Traitement des fichiers 1D (%d/%d)', i, total_files));
    end
    
    close(h);
else
    disp('Aucun fichier 1D trouvé.');
end

disp('Traitement terminé pour tous les fichiers 1D.');


%% 2D Smoothing
disp('Début du traitement des fichiers 2D...');
if ~isempty(ph_files_2D)
    h = waitbar(0, 'Traitement des fichiers 2D...');
    total_files = length(ph_files_2D);
    
    for i = 1:total_files
        data_2D_file = fullfile(current_path, ph_files_2D(i).name); 
        proba_2D_file = fullfile(current_path, proba_files_2D(i).name);
        
        fprintf('Chargement des données 2D depuis %s...\n', data_2D_file);
        load(data_2D_file, 'data_2D');

        fprintf('Chargement des probabilités tissulaires 2D depuis %s...\n', proba_2D_file);
        load(proba_2D_file, 'data_GmWmCsfSculpt_2D');

        % Exemple de traitement (à adapter selon tes besoins de lissage)
        % Appliquer des opérations sur les données ici...
        % Exemple : gsP_signal = smooth(data_2D);
        % Exemple : twsP_signal = smooth(data_2D, 'lowess');
        % Exemple : tosP_signal = some_other_smoothing_function(data_2D);

        if flag.save_data
            output_filename = sprintf('results_2D_%s.mat', ph_files_2D(i).name(1:end-4));
            save(output_filename, 'gsP_signal', 'twsP_signal', 'tosP_signal');
        end
        
        waitbar(i / total_files, h, sprintf('Traitement des fichiers 2D (%d/%d)', i, total_files));
    end
    
    close(h);
else
    disp('Aucun fichier 2D trouvé.');
end

disp('Traitement terminé pour tous les fichiers 2D.');


%% 3D Smoothing
disp('Début du traitement des fichiers 3D...');
if ~isempty(ph_files_3D)
    h = waitbar(0, 'Traitement des fichiers 3D...');
    total_files = length(ph_files_3D);
    
    for i = 1:total_files
        data_3D_file = fullfile(current_path, ph_files_3D(i).name); 
        proba_3D_file = fullfile(current_path, proba_files_3D(i).name);
        
        fprintf('Chargement des données 3D depuis %s...\n', data_3D_file);
        load(data_3D_file, 'data_3D');

        fprintf('Chargement des probabilités tissulaires 3D depuis %s...\n', proba_3D_file);
        load(proba_3D_file, 'data_GmWmCsfSculpt_3D');

        % Exemple de traitement (à adapter selon tes besoins de lissage)
        % Appliquer des opérations sur les données ici...
        % Exemple : gsP_signal = smooth(data_3D);
        % Exemple : twsP_signal = smooth(data_3D, 'lowess');
        % Exemple : tosP_signal = some_other_smoothing_function(data_3D);

        if flag.save_data
            output_filename = sprintf('results_3D_%s.mat', ph_files_3D(i).name(1:end-4));
            save(output_filename, 'gsP_signal');
        end
        
        waitbar(i / total_files, h, sprintf('Traitement des fichiers 3D (%d/%d)', i, total_files));
    end
    
    close(h);
else
    disp('Aucun fichier 3D trouvé.');
end

disp('Traitement terminé pour tous les fichiers 3D.');
