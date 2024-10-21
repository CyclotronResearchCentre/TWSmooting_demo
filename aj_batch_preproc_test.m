% Initialisation
nsub = 1; % Nombre de sujets
nrun = 9; % Nombre de runs (ici 9 basés sur vos chemins)
outpth = 'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\work_copies\openneuro.org\ds000117'; % Chemin de sortie
subdir = {'sub-02'}; % Liste des sujets

% Chemin du fichier de job
jobfile = {'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\aj_batch_preproc_test_job.m'};
jobs = repmat(jobfile, 1, nsub);

% Préparation des inputs
inputs = cell(nrun + 2, 1);

% Boucle sur les sujets
for s = 1:nsub
    % Changer le répertoire au répertoire du sujet
%     cd(fullfile(outpth, subdir{s}, 'ses-mri', 'func'));

    for r = 1:nrun
        inputs{r} = cellstr(spm_select('FPList', ...
            fullfile(outpth, subdir{s}, 'ses-mri', 'func'), ...
            sprintf('^%s_ses-mri_task-facerecognition_run-%02d_bold\\.nii$', subdir{s}, r)));
    end

   inputs{nrun + 1} = cellstr(spm_select('FPList', ...
        fullfile(outpth, subdir{s}, 'ses-mri', 'anat'), ...
        sprintf('^%s_ses-mri_acq-mprage_T1w\\.nii$', subdir{s})));

    % Vérifiez que nous avons récupéré les bons fichiers
    disp('Functional files:');
    disp(inputs(1:nrun));
    disp('Anatomical file:');
    disp(inputs{nrun + 1});
    
    % Vérification des fichiers fonctionnels
    for r = 1:nrun
        if isempty(inputs{r})
            error('Functional file for run %d is missing.', r);
        end
    end

    % Vérification du fichier anatomique
    if isempty(inputs{nrun + 1})
        error('No anatomical file');
    end
    
    % Exécution du job SPM
    spm_jobman('run', jobs{s}, inputs{:});
end

% Réinitialiser SPM pour chaque sujet
spm('defaults', 'FMRI');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Initialisation
% nrun = 9; % Nombre de runs (ici 9 basés sur vos chemins)
% outpth = 'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\work_copies\openneuro.org\ds000117'; % Chemin de sortie
% 
% % Récupération des sujets à partir de spm_BIDS
% subdir = spm_BIDS(outpth); % Assurez-vous que spm_BIDS retourne la liste des sous-répertoires
% nsub = numel(subdir); % Nombre de sujets
% 
% % Chemin du fichier de job
% jobfile = {'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\aj_batch_preproc_test_job.m'};
% jobs = repmat(jobfile, 1, nsub);
% 
% % Préparation des inputs
% inputs = cell(nrun + 2, 1);
% 
% % Parallélisation
% parfor s = 1:nsub
%     % Changer le répertoire au répertoire du sujet
%     cd(fullfile(outpth, subdir{s}, 'ses-mri', 'func'));
%     
%     % Récupération des fichiers fonctionnels
%     for r = 1:nrun
%         inputs{r} = cellstr(spm_select('FPList', ...
%             fullfile(outpth, subdir{s}, 'ses-mri', 'func'), ...
%             sprintf('^%s_ses-mri_task-facerecognition_run-%02d_bold\\.nii$', subdir{s}, r)));
%     end
%     
%     % Récupération des fichiers anatomiques
%     inputs{nrun + 1} = cellstr(spm_select('FPList', ...
%         fullfile(outpth, subdir{s}, 'ses-mri', 'anat'), ...
%         sprintf('^%s_ses-mri_acq-mprage_T1w\\.nii$', subdir{s})));
%     
%     % Vérifiez que nous avons récupéré les bons fichiers
%     disp(['Functional files for ', subdir{s}, ':']);
%     disp(inputs(1:nrun));
%     disp(['Anatomical file for ', subdir{s}, ':']);
%     disp(inputs{nrun + 1});
%     
%     % Exécution du job SPM
%     spm_jobman('run', jobs{s}, inputs{:});
% end
% 
% % Réinitialiser SPM pour chaque sujet
% spm('defaults', 'FMRI');
% 
% % Fermer le pool parallèle à la fin
% delete(gcp('nocreate'));