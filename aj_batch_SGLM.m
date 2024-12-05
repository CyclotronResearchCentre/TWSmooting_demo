
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');

% List of open inputs
nrun = 1; % enter the number of runs here
jobfile = {'C:\Users\antoi\Documents\master_thesis\MATLAB\agingdata\aj_batch_SGLM_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
