% Main script to download and decompress BIDS files in parallel

addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');

%% Step 0: get URL and target_dir
[urls, target_directory] = aj_download_default();

% Vérifier si le répertoire cible existe, sinon le créer
if ~exist(target_directory, 'dir')
    mkdir(target_directory);
end

% Activer le pool parallèle
if isempty(gcp('nocreate'))
    parpool;  % Lancer un pool parallèle si aucun n'est actif
end

%% Téléchargement des fichiers en parallèle
fprintf('Starting parallel download...\n');
parfor i = 1:length(urls)
    url = urls{i};
    % Extraire le chemin relatif depuis l'URL
    file_path = urlsplit(url);
    % Créer le chemin complet de destination
    dest_file = fullfile(target_directory, file_path);
    % Créer les dossiers nécessaires pour le fichier
    dest_folder = fileparts(dest_file);
    if ~exist(dest_folder, 'dir')
        mkdir(dest_folder);
    end

    % Télécharger le fichier
    try
        fprintf('Downloading %s...\n', url);
        websave(dest_file, url);  % Téléchargement du fichier
        fprintf('Successfully downloaded to %s\n', dest_file);
    catch ME
        fprintf('Error downloading %s: %s\n', url, ME.message);
    end
end

%% Décompression des fichiers .gz en parallèle
fprintf('Starting parallel decompression...\n');
file_list = unzip_gz_files_parallel(target_directory);

%% HELP FUNCTIONS

% Fonction pour extraire le chemin relatif depuis une URL
function file_path = urlsplit(url)
    % Supprimer les paramètres de l'URL (tout ce qui suit '?')
    base_url = strtok(url, '?');
    
    % Extraire la partie après l'hôte (path relatif)
    split_url = split(base_url, '/');
    file_path = fullfile(split_url{4:end}); % Ignorer les trois premiers éléments (protocole + hôte)
end

% Fonction pour décompresser les fichiers .gz en parallèle
function file_list = unzip_gz_files_parallel(target_directory)
    % Cette fonction décompresse tous les fichiers GZ dans le répertoire cible 
    % et retourne une liste des fichiers décompressés.
    
    % Obtenir une liste de tous les fichiers .gz dans l'arborescence
    gz_files = dir(fullfile(target_directory, '**', '*.gz'));
    
    % Initialiser une cellule pour stocker les chemins des fichiers décompressés
    file_list = {};
    
    % Décompression en parallèle
    parfor i = 1:length(gz_files)
        gz_file_path = fullfile(gz_files(i).folder, gz_files(i).name);
        output_file = gz_file_path(1:end-3);  % Supprimer l'extension .gz pour le fichier de sortie
        
        try
            fprintf('Unzipping %s...\n', gz_file_path);
            gunzip(gz_file_path);  % Décompresser le fichier
            
            % Supprimer le fichier .gz après extraction (optionnel)
            delete(gz_file_path);
            
            % Ajouter le chemin du fichier décompressé à la liste
            file_list{i} = output_file;
        catch ME
            fprintf('Error unzipping %s: %s\n', gz_file_path, ME.message);
        end
    end
end
