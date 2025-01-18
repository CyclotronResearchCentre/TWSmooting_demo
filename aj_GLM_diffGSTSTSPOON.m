% This script compute a GLM (Factorial Design) to spot the regions where
% the smoothed contrast maps are affected by the smoothing approach. The
% used contrast map is con0002 (Faces > Scrambled Faces) on which the three
% smoothing approaches have been applied (GS, TWS and TSPOON). This
% contrast map is obtained from the first-level GLM computed on the
% preprocessed (without smoothing) fMRI data.
%
% First the script computes GM and WM masks over the population using the
% "Majority and Greater than 20%" criterion. Then the Factorial Design is
% applied followed by the F and T contrasts.
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------
%% Create GM and WM masks for the population
% Take all mwc for each TC
% Mean them -> one mean mwc per TC
% Apply a mask by the "majority and above .2" criteria

clear;clc;

TCs = {'GM', 'WM', 'CSF'};
nTC = numel(TCs);

preproc_dir = 'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\work_copies\openneuro.org\ds000117\derivatives\preprocessing';

mwc_paths = cell(nTC,1);
for i = 1:nTC
    pattern = sprintf('^mwc%d.*ses-mri_acq-mprage_T1w\\.nii$',i);
    mwc_paths{i} = spm_select('FPListRec', preproc_dir, pattern);
end

mean_images = cell(1,nTC);
for i = 1:nTC
    % Charger les images mwc pour ce tissu
    files = mwc_paths{i};
    nFiles = size(files, 1);
    all_images = [];
    
    for j = 1:nFiles
        img = spm_read_vols(spm_vol(files(j,:))); % Lecture des données
        if isempty(all_images)
            all_images = zeros(size(img,1), size(img,2), size(img,3), nFiles);
        end
        all_images(:,:,:,j) = img;
    end
    
    % Calculer la moyenne des images
    mean_image = mean(all_images, 4);
    
    % Sauvegarder le masque comme fichier NIfTI
    out_file = fullfile(preproc_dir, sprintf('mean_mwc_%s.nii', TCs{i}));
    vol = spm_vol(files(1,:));
    vol.fname = out_file;
    spm_write_vol(vol, mean_image);
    mean_images{i} = out_file;
    
    fprintf('%s mean saved in: %s\n', TCs{i}, out_file);
end

addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\smoothing');
exMask_paths = aj_exMask(mean_images, preproc_dir);


%% Compute GLM diff: GS - TWS - TSPOON
clear;clc;
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');

% Initialize SPM configuration
spm_jobman('initcfg'); 
spm('defaults', 'fmri');

out_dir = 'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\work_copies\openneuro.org\ds000117\derivatives\AJ-diff_GSTWSTSPOON';

GS_dir = 'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\work_copies\openneuro.org\ds000117\derivatives\AJ-GS';
TWS_dir = 'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\work_copies\openneuro.org\ds000117\derivatives\AJ-TWS';
TSPOON_dir = 'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\work_copies\openneuro.org\ds000117\derivatives\AJ-TSPOON';
mask_dir = 'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\work_copies\openneuro.org\ds000117\derivatives\preprocessing';

smoothing_methods = {'GS', 'TWS', 'TSPOON'};
fcontrast = {'con_0002'};
TCs = {'GM', 'WM'};

% Combining parameters using combvec (Cartesian product)
[TCs_idx, fMRI_idx] = ndgrid(1:numel(TCs), 1:numel(fcontrast));

% Creating GS, TWS and TSPOON patterns
GS_patterns = cellstr(sprintf('^.*%s\\.nii$', fcontrast{1}));
TWS_patterns = arrayfun(@(i, j) sprintf('^%s_%s.*%s\\.nii$', ...
                        smoothing_methods{2}, TCs{i}, fcontrast{j}), ...
                        TCs_idx(:), fMRI_idx(:), 'UniformOutput', false);
TSPOON_patterns = arrayfun(@(i, j) sprintf('^%s_%s.*%s\\.nii$', ...
                        smoothing_methods{3}, TCs{i}, fcontrast{j}), ...
                        TCs_idx(:), fMRI_idx(:), 'UniformOutput', false);
Mask_patterns = arrayfun(@(i) sprintf('^Mask_mean_mwc_%s\\.nii$',...
                        TCs{i}),...
                        TCs_idx(:), 'UniformOutput', false);

for i = 1:numel(TCs)
    jobfile = {'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\aj_GLM_diffGSTSTSPOON_job.m'};
    
    inputs  = cell(5,1);
    % Output filename
    inputs{1} = cellstr(fullfile(out_dir, TCs{i})); % Factorial design specification: Directory - cfg_files;
    % GS files
    inputs{2} = cellstr(spm_select('FPListRec', GS_dir, GS_patterns)); % Factorial design specification: Scans - cfg_files
    % TWS files
    inputs{3} = cellstr(spm_select('FPListRec', TWS_dir, TWS_patterns{i})); % Factorial design specification: Scans - cfg_files
    % TSPOON files
    inputs{4} = cellstr(spm_select('FPListRec', TSPOON_dir, TSPOON_patterns{i})); % Factorial design specification: Scans - cfg_files
    % Mask file
    inputs{5} = cellstr(spm_select('ExtFPList', mask_dir, Mask_patterns{i})); % Factorial design specification: Explicit Mask - cfg_files
    
    spm_jobman('run', jobfile, inputs{:});
end
