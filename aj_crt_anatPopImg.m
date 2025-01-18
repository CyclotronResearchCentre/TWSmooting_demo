% Script to create anatomical mean image over the population (16 subjects)
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------
% Initialisation
clear; clc;
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12'); % Ajouter SPM au chemin

% Dossier contenant les images anatomiques
anat_dir = 'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\work_copies\openneuro.org\ds000117';
output_dir = fullfile(anat_dir, 'derivatives'); % Dossier de sortie pour l'image moyenne
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Liste des numéros de sujets (1 à 16)
num_subjects = 16;
subjects = arrayfun(@(x) sprintf('sub-%02d', x), 1:num_subjects, 'UniformOutput', false);

% Initialisation des variables pour le calcul moyen
sum_image = [];
count = 0;

% Options pour SPM
spm('defaults', 'FMRI');
spm_jobman('initcfg');

% Prétraitement et calcul de l'image moyenne
for i = 1:num_subjects
    % Chemin de l'image anatomique pour le sujet courant
    anat_path = fullfile(anat_dir, subjects{i}, 'ses-mri', 'anat', sprintf('%s_ses-mri_acq-mprage_T1w.nii', subjects{i}));
    
    % Vérification de l'existence du fichier
    if exist(anat_path, 'file')
        % Normalisation de l'image dans l'espace MNI avec SPM
        fprintf('Normalisation de %s...\n', anat_path);
        
        matlabbatch = {};
        matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = {anat_path}; % Image anatomique
        matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = {anat_path};
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm = {fullfile(spm('Dir'), 'tpm', 'TPM.nii')};
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
        matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70; 78 76 85];
        matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = [1 1 1];
        matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.interp = 4;
        
        spm_jobman('run', matlabbatch);
        
        % Image normalisée
        [pathstr, name, ext] = fileparts(anat_path);
        norm_anat_path = fullfile(pathstr, ['w' name ext]); % "w" préfixe ajouté par SPM
        
        % Lecture de l'image normalisée avec SPM
        V = spm_vol(norm_anat_path); % Métadonnées
        img = spm_read_vols(V); % Chargement des données en mémoire
        
        % Initialisation de la somme si c'est la première image
        if isempty(sum_image)
            sum_image = zeros(size(img)); % Initialisation avec les dimensions de l'image
        end
        
        % Ajout de l'image à la somme
        sum_image = sum_image + img;
        count = count + 1;
    else
        fprintf('Fichier non trouvé pour %s\n', subjects{i});
    end
end

% Calcul de l'image moyenne
if count > 0
    avg_image = sum_image / count;
    
    % Création d'une structure pour enregistrer l'image moyenne
    V_avg = V; % Copier les métadonnées du dernier fichier lu
    V_avg.fname = fullfile(output_dir, 'avg_T1w_MNI.nii'); % Chemin de l'image moyenne
    spm_write_vol(V_avg, avg_image); % Enregistrement de l'image moyenne
    fprintf('Image moyenne enregistrée dans %s\n', V_avg.fname);
else
    error('Aucune image n’a été trouvée ou normalisée. Vérifiez les chemins.');
end
