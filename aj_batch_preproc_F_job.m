%-----------------------------------------------------------------------
% Job saved on 15-Oct-2024 10:47:06 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
%%
matlabbatch{1}.spm.spatial.realign.estwrite.data = {
                                                    {'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\work_copies\openneuro.org\ds000117\sub-02\ses-mri\func\sub-02_ses-mri_task-facerecognition_run-01_bold.nii,1'}
                                                    {'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\work_copies\openneuro.org\ds000117\sub-02\ses-mri\func\sub-02_ses-mri_task-facerecognition_run-02_bold.nii,1'}
                                                    {'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\work_copies\openneuro.org\ds000117\sub-02\ses-mri\func\sub-02_ses-mri_task-facerecognition_run-03_bold.nii,1'}
                                                    {'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\work_copies\openneuro.org\ds000117\sub-02\ses-mri\func\sub-02_ses-mri_task-facerecognition_run-04_bold.nii,1'}
                                                    {'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\work_copies\openneuro.org\ds000117\sub-02\ses-mri\func\sub-02_ses-mri_task-facerecognition_run-05_bold.nii,1'}
                                                    {'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\work_copies\openneuro.org\ds000117\sub-02\ses-mri\func\sub-02_ses-mri_task-facerecognition_run-06_bold.nii,1'}
                                                    {'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\work_copies\openneuro.org\ds000117\sub-02\ses-mri\func\sub-02_ses-mri_task-facerecognition_run-07_bold.nii,1'}
                                                    {'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\work_copies\openneuro.org\ds000117\sub-02\ses-mri\func\sub-02_ses-mri_task-facerecognition_run-08_bold.nii,1'}
                                                    {'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\work_copies\openneuro.org\ds000117\sub-02\ses-mri\func\sub-02_ses-mri_task-facerecognition_run-09_bold.nii,1'}
                                                    }';
%%
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [0 1];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
matlabbatch{2}.spm.spatial.preproc.channel.vols = {'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\work_copies\openneuro.org\ds000117\sub-02\ses-mri\anat\sub-02_ses-mri_acq-mprage_T1w.nii,1'};
matlabbatch{2}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{2}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{2}.spm.spatial.preproc.channel.write = [0 0];
matlabbatch{2}.spm.spatial.preproc.tissue(1).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,1')};
matlabbatch{2}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{2}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{2}.spm.spatial.preproc.tissue(1).warped = [0 1];
matlabbatch{2}.spm.spatial.preproc.tissue(2).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,2')};
matlabbatch{2}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{2}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{2}.spm.spatial.preproc.tissue(2).warped = [0 1];
matlabbatch{2}.spm.spatial.preproc.tissue(3).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,3')};
matlabbatch{2}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{2}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{2}.spm.spatial.preproc.tissue(3).warped = [0 1];
matlabbatch{2}.spm.spatial.preproc.tissue(4).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,4')};
matlabbatch{2}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{2}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{2}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{2}.spm.spatial.preproc.tissue(5).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,5')};
matlabbatch{2}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{2}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{2}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{2}.spm.spatial.preproc.tissue(6).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,6')};
matlabbatch{2}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{2}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{2}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{2}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{2}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{2}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{2}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{2}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{2}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{2}.spm.spatial.preproc.warp.write = [0 1];
matlabbatch{2}.spm.spatial.preproc.warp.vox = NaN;
matlabbatch{2}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
                                              NaN NaN NaN];
matlabbatch{3}.spm.spatial.coreg.estimate.ref = {'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\work_copies\openneuro.org\ds000117\sub-02\ses-mri\anat\sub-02_ses-mri_acq-mprage_T1w.nii,1'};
matlabbatch{3}.spm.spatial.coreg.estimate.source(1) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
matlabbatch{3}.spm.spatial.coreg.estimate.other(1) = cfg_dep('Realign: Estimate & Reslice: Realigned Images (Sess 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','cfiles'));
matlabbatch{3}.spm.spatial.coreg.estimate.other(2) = cfg_dep('Realign: Estimate & Reslice: Realigned Images (Sess 2)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{2}, '.','cfiles'));
matlabbatch{3}.spm.spatial.coreg.estimate.other(3) = cfg_dep('Realign: Estimate & Reslice: Realigned Images (Sess 3)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{3}, '.','cfiles'));
matlabbatch{3}.spm.spatial.coreg.estimate.other(4) = cfg_dep('Realign: Estimate & Reslice: Realigned Images (Sess 4)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{4}, '.','cfiles'));
matlabbatch{3}.spm.spatial.coreg.estimate.other(5) = cfg_dep('Realign: Estimate & Reslice: Realigned Images (Sess 5)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{5}, '.','cfiles'));
matlabbatch{3}.spm.spatial.coreg.estimate.other(6) = cfg_dep('Realign: Estimate & Reslice: Realigned Images (Sess 6)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{6}, '.','cfiles'));
matlabbatch{3}.spm.spatial.coreg.estimate.other(7) = cfg_dep('Realign: Estimate & Reslice: Realigned Images (Sess 7)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{7}, '.','cfiles'));
matlabbatch{3}.spm.spatial.coreg.estimate.other(8) = cfg_dep('Realign: Estimate & Reslice: Realigned Images (Sess 8)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{8}, '.','cfiles'));
matlabbatch{3}.spm.spatial.coreg.estimate.other(9) = cfg_dep('Realign: Estimate & Reslice: Realigned Images (Sess 9)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{9}, '.','cfiles'));
matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
matlabbatch{4}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
matlabbatch{4}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
matlabbatch{4}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                          78 76 85];
matlabbatch{4}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
matlabbatch{4}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{4}.spm.spatial.normalise.write.woptions.prefix = 'w';

% matlabbatch{5}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
% matlabbatch{5}.spm.spatial.smooth.fwhm = [8 8 8]; % standard gaussian smoothing with kernel size of 8 mm
% matlabbatch{5}.spm.spatial.smooth.dtype = 0; % no file extension convertion between input and output
% matlabbatch{5}.spm.spatial.smooth.im = 0; % applying no implicit mask
% matlabbatch{5}.spm.spatial.smooth.prefix = 's';
