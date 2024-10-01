% Script principal pour exécuter les différentes méthodes de lissage sur
% des données 1D, 2D et 3D

%% Step 0 : cleaning
clear;
clc;
close all;

%% Step 1: Récupérer tous les fichiers 1D, 2D et 3D disponibles
current_path = pwd;
ph_files_1D = dir(fullfile(current_path, 'phantom_1D_*.mat'));
proba_files_1D = dir(fullfile(current_path, 'proba_1D_*.mat'));
ph_files_2D = dir(fullfile(current_path, 'phantom_2D_*.mat'));
proba_files_2D = dir(fullfile(current_path, 'proba_2D_*.mat'));
ph_files_3D = dir(fullfile(current_path, 'phantom_3D_*.nii'));
proba_files_3D = dir(fullfile(current_path, 'proba_3D_*.mat'));

if isempty(ph_files_1D)
    error('Aucun fichier 1D trouvé.');
elseif isempty(ph_files_2D)
    error('Aucun fichier 2D trouvé.');
elseif isempty(ph_files_3D)
    error('Aucun fichier 3D trouvé.');
end

if length(ph_files_1D) ~= length(proba_files_1D)
    error('Le nombre de fichiers de données 1D et de fichiers de probabilités 1D ne correspond pas.');
elseif length(ph_files_2D) ~= length(proba_files_2D)
    error('Le nombre de fichiers de données 2D et de fichiers de probabilités 2D ne correspond pas.');
elseif length(ph_files_3D) ~= length(proba_files_3D)
    error('Le nombre de fichiers de données 3D et de fichiers de probabilités 3D ne correspond pas.');
end

fprintf('Nombre de fichiers 1D trouvés : %d\n', length(ph_files_1D));
fprintf('Nombre de fichiers 2D trouvés : %d\n', length(ph_files_2D));
fprintf('Nombre de fichiers 3D trouvés : %d\n', length(ph_files_3D));

%% Default parameters
% Paramètres et flags communs de lissage (paramètres ajustables)
[param, flag] = aj_smooth_default();

%% 1D Smoothing
% Boucle à travers tous les fichiers 1D
for i = 1:length(ph_files_1D)
    % Charger les données 1D
    data_1D_file = fullfile(current_path, ph_files_1D(i).name);
    proba_1D_file = fullfile(current_path, proba_files_1D(i).name);
    
    fprintf('Chargement des données 1D depuis %s...\n', data_1D_file);
    load(data_1D_file, 'data_1D');
    
    fprintf('Chargement des probabilités tissulaires 1D depuis %s...\n', proba_1D_file);
    load(proba_1D_file, 'data_GmWmCsfSculpt_1D');
    
    % Vérifier la taille des données
%     disp(['Taille des données 1D : ', num2str(size(data_1D))]);
%     disp(['Taille des probabilités tissulaires 1D : ', num2str(size(data_GmWmCsfSculpt_1D))]);

    % --------------------------------------------------------------------
    % 1. Appliquer le lissage Gaussien standard
    % --------------------------------------------------------------------
    disp('Exécution du lissage Gaussien standard...');
    gsP_signal = aj_smooth_gaussian(data_1D, param);

    % --------------------------------------------------------------------
    % 2. Appliquer le lissage pondéré par les tissus (TWS)
    % --------------------------------------------------------------------
    disp('Exécution du lissage pondéré par les tissus...');
    twsP_signal = aj_smooth_TWS(data_1D, data_GmWmCsfSculpt_1D, param);

    % --------------------------------------------------------------------
    % 3. Appliquer le lissage TSPOON
    % --------------------------------------------------------------------
    disp('Exécution du lissage TSPOON...');
    tosP_signal = aj_smooth_TSPOON(data_1D, data_GmWmCsfSculpt_1D, param);

    % --------------------------------------------------------------------
    % Afficher les résultats pour ce fichier
    % --------------------------------------------------------------------
    if flag.plot_fig
        figure('Name', sprintf('Résultats du fichier : %s', ph_files_1D(i).name));

        subplot(3, 1, 1);
        plot(data_1D, 'b-', 'DisplayName', 'Données brutes');
        hold on;
        plot(gsP_signal, 'r-', 'DisplayName', 'Lissage Gaussien');
        legend;
        xlabel('Position');
        ylabel('Valeur du signal');
        title('Lissage Gaussien');

        subplot(3, 1, 2);
        plot(data_1D, 'b-', 'DisplayName', 'Données brutes');
        hold on;
        plot(twsP_signal(1,:), 'r-', 'DisplayName', 'Lissage Pondéré (GM)');
        plot(twsP_signal(2,:), 'g-', 'DisplayName', 'Lissage Pondéré (WM)');
        plot(twsP_signal(3,:), 'm-', 'DisplayName', 'Lissage Pondéré (CSF)');
        plot(twsP_signal(4,:), 'y-', 'DisplayName', 'Lissage Pondéré (Sculpt)');
        legend;
        xlabel('Position');
        ylabel('Valeur du signal');
        title('Lissage Pondéré par les Tissus');

        subplot(3, 1, 3);
        plot(data_1D, 'b-', 'DisplayName', 'Données brutes');
        hold on;
        plot(tosP_signal(1,:), 'r-', 'DisplayName', 'Lissage TSPOON (GM)');
        plot(tosP_signal(2,:), 'g-', 'DisplayName', 'Lissage TSPOON (WM)');
        plot(tosP_signal(3,:), 'm-', 'DisplayName', 'Lissage TSPOON (CSF)');
        legend;
        xlabel('Position');
        ylabel('Valeur du signal');
        title('Lissage TSPOON');

        hold off;
    end

    % Sauvegarder les résultats dans un fichier .mat spécifique à ce fichier
    if flag.save_data
        output_filename = sprintf('results_%s.mat', ph_files_1D(i).name(1:end-4));
        save(output_filename, 'gsP_signal', 'twsP_signal', 'tosP_signal');
    end
    
    fprintf('Résultats sauvegardés dans %s\n', output_filename);
end

disp('Traitement terminé pour tous les fichiers 1D.');


%% 2D Processing
% Boucle pour traiter les fichiers 2D
for i = 1:length(ph_files_2D)
    % Charger les données 2D
    data_2D_file = fullfile(current_path, ph_files_2D(i).name);
    proba_2D_file = fullfile(current_path, proba_files_2D(i).name);
    
    fprintf('Chargement des données 2D depuis %s...\n', data_2D_file);
    load(data_2D_file, 'data_2D');
    
    fprintf('Chargement des probabilités tissulaires 2D depuis %s...\n', proba_2D_file);
    load(proba_2D_file, 'data_GmWmCsfSculpt_2D');
    
    % Vérifier la taille des données
%     disp(['Taille des données 2D : ', num2str(size(data_2D))]);
%     disp(['Taille des probabilités tissulaires 2D : ', num2str(size(data_GmWmCsfSculpt_2D))]);

    % --------------------------------------------------------------------
    % 1. Appliquer le lissage Gaussien standard
    % --------------------------------------------------------------------
    disp('Exécution du lissage Gaussien standard...');
    gsP_signal = aj_smooth_gaussian(data_2D, param);

    % --------------------------------------------------------------------
    % 2. Appliquer le lissage pondéré par les tissus (TWS)
    % --------------------------------------------------------------------
    disp('Exécution du lissage pondéré par les tissus...');
    twsP_signal = aj_smooth_TWS(data_2D, data_GmWmCsfSculpt_2D, param);

    % --------------------------------------------------------------------
    % 3. Appliquer le lissage TSPOON
    % --------------------------------------------------------------------
    disp('Exécution du lissage TSPOON...');
    tosP_signal = aj_smooth_TSPOON(data_2D, data_GmWmCsfSculpt_2D, param);

    % --------------------------------------------------------------------
    % Afficher les résultats pour ce fichier
    % --------------------------------------------------------------------
    if flag.plot_fig
        figure('Name', sprintf('Résultats du fichier : %s', ph_files_2D(i).name));

        % Lissage Gaussien
        subplot(3, 1, 1);
        imagesc(gsP_signal);  % Affichage des données 2D
        axis image;  % Pour maintenir le rapport d'aspect
        colorbar;  % Ajouter une barre de couleur pour la légende
        title('Lissage Gaussien');
        xlabel('Position X');
        ylabel('Position Y');

        % Lissage pondéré par les tissus (TWS)
        subplot(3, 1, 2);
        imagesc(twsP_signal(:,:,1));  % Affichage pour le premier tissu (GM)
        axis image;
        colorbar;
        title('Lissage Pondéré par les Tissus (GM)');
        xlabel('Position X');
        ylabel('Position Y');

        % Vous pouvez ajouter des sous-graphes pour les autres tissus (WM, CSF, Sculpt)
        % Par exemple, pour le tissu WM :
        % subplot(3, 2, 3);
        % imagesc(twsP_signal(:,:,2));  % WM
        % axis image;
        % colorbar;
        % title('Lissage Pondéré par les Tissus (WM)');
        % xlabel('Position X');
        % ylabel('Position Y');

        % Lissage TSPOON
        subplot(3, 1, 3);
        imagesc(tosP_signal(:,:,1));  % Affichage pour le premier tissu (GM)
        axis image;
        colorbar;
        title('Lissage TSPOON (GM)');
        xlabel('Position X');
        ylabel('Position Y');

        % Ajouter des sous-graphes pour les autres tissus (WM, CSF)
        % subplot(3, 2, 4);
        % imagesc(tosP_signal(:,:,2));  % WM
        % axis image;
        % colorbar;
        % title('Lissage TSPOON (WM)');
    end

    % Sauvegarder les résultats dans un fichier .mat spécifique à ce fichier
    if flag.save_data
        output_filename = sprintf('results_%s.mat', ph_files_2D(i).name(1:end-4));
        save(output_filename, 'gsP_signal', 'twsP_signal', 'tosP_signal');
    end
    
    fprintf('Résultats sauvegardés dans %s\n', output_filename);
end

disp('Traitement terminé pour tous les fichiers 2D.');


%% 3D Smoothing
% Boucle pour traiter les fichiers 3D
% for i = 1:length(files_3D)
%     data_3D_file = fullfile(current_path, files_3D(i).name);
%     V = spm_vol(data_3D_file);  % Charger les données NIfTI
%     data_3D = spm_read_vols(V);  % Lire les volumes
% 
%     proba_3D_file = fullfile(current_path, strrep(files_3D(i).name, 'phantom_3D_', 'proba_3D_'));
%     V_proba = spm_vol(proba_3D_file);
%     data_GmWmCsfSculpt_3D = spm_read_vols(V_proba);  % Probabilités tissulaires 3D
% 
%     % Lissage Pondéré par les Tissus pour données 3D
%     twsP_signal = aj_smooth_TWS(data_3D, data_GmWmCsfSculpt_3D, param, flag);
% end

for i = 1:length(ph_files_3D)
    % Charger les données 3D
    data_3D_file = fullfile(current_path, ph_files_3D(i).name);
    proba_3D_file = fullfile(current_path, proba_files_3D(i).name);
    
    fprintf('Chargement des données 3D depuis %s...\n', data_3D_file);
    load(data_3D_file, 'data_3D');
    
    fprintf('Chargement des probabilités tissulaires 3D depuis %s...\n', proba_3D_file);
    load(proba_3D_file, 'data_GmWmCsfSculpt_3D');
    
    % Vérifier la taille des données
%     disp(['Taille des données 3D : ', num2str(size(data_3D))]);
%     disp(['Taille des probabilités tissulaires 3D : ', num2str(size(data_GmWmCsfSculpt_3D))]);

    % --------------------------------------------------------------------
    % 1. Appliquer le lissage Gaussien standard
    % --------------------------------------------------------------------
    disp('Exécution du lissage Gaussien standard...');
    gsP_signal = aj_smooth_gaussian(data_3D, param);

    % --------------------------------------------------------------------
    % 2. Appliquer le lissage pondéré par les tissus (TWS)
    % --------------------------------------------------------------------
    disp('Exécution du lissage pondéré par les tissus...');
    % twsP_signal = aj_smooth_TWS(data_3D, data_GmWmCsfSculpt_3D, param);

    % --------------------------------------------------------------------
    % 3. Appliquer le lissage TSPOON
    % --------------------------------------------------------------------
    disp('Exécution du lissage TSPOON...');
    % tosP_signal = aj_smooth_TSPOON(data_3D, data_GmWmCsfSculpt_3D, param);

    % --------------------------------------------------------------------
    % Afficher les résultats pour ce fichier
    % --------------------------------------------------------------------
    if flag.plot_fig
        figure('Name', sprintf('Résultats du fichier : %s', ph_files_3D(i).name));

        % Lissage Gaussien
        subplot(3, 1, 1);
        imagesc(gsP_signal);  % Affichage des données 2D
        axis image;  % Pour maintenir le rapport d'aspect
        colorbar;  % Ajouter une barre de couleur pour la légende
        title('Lissage Gaussien');
        xlabel('Position X');
        ylabel('Position Y');

        % Lissage pondéré par les tissus (TWS)
        subplot(3, 1, 2);
        imagesc(twsP_signal(:,:,1));  % Affichage pour le premier tissu (GM)
        axis image;
        colorbar;
        title('Lissage Pondéré par les Tissus (GM)');
        xlabel('Position X');
        ylabel('Position Y');

        % Vous pouvez ajouter des sous-graphes pour les autres tissus (WM, CSF, Sculpt)
        % Par exemple, pour le tissu WM :
        % subplot(3, 2, 3);
        % imagesc(twsP_signal(:,:,2));  % WM
        % axis image;
        % colorbar;
        % title('Lissage Pondéré par les Tissus (WM)');
        % xlabel('Position X');
        % ylabel('Position Y');

        % Lissage TSPOON
        subplot(3, 1, 3);
        imagesc(tosP_signal(:,:,1));  % Affichage pour le premier tissu (GM)
        axis image;
        colorbar;
        title('Lissage TSPOON (GM)');
        xlabel('Position X');
        ylabel('Position Y');

        % Ajouter des sous-graphes pour les autres tissus (WM, CSF)
        % subplot(3, 2, 4);
        % imagesc(tosP_signal(:,:,2));  % WM
        % axis image;
        % colorbar;
        % title('Lissage TSPOON (WM)');
    end

    % Sauvegarder les résultats dans un fichier .mat spécifique à ce fichier
    if flag.save_data
        output_filename = sprintf('results_%s.mat', ph_files_3D(i).name(1:end-4));
        save(output_filename, 'gsP_signal', 'twsP_signal', 'tosP_signal');
    end
    
    fprintf('Résultats sauvegardés dans %s\n', output_filename);
end

disp('Traitement terminé pour tous les fichiers 3D.');