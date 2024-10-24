% Create SPM's multiple conditions files
% Codes: Multimodal integration of M/EEG and f/MRI data in SPM12 of gllmflndn ; https://github.com/spm/MultimodalScripts/tree/master
%--------------------------------------------------------------------------
clear;clc;

base_dir = 'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117';
work_bids_root = fullfile(base_dir, 'work_copies', 'openneuro.org', 'ds000117');

runs = spm_BIDS(work_bids_root,'runs', 'modality','func', 'type','bold', 'task','facerecognition'); 
nrun = numel(runs);

subs = spm_BIDS(work_bids_root, 'subjects');
nsub = numel(subs);

trialtypes = {'Famous','Unfamiliar','Scrambled'}; % impose order

outpth = 'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\work_copies\openneuro.org\ds000117\derivatives\preprocessing';

for s = 1:1%nsub
    for r = 1:nrun  
            d = spm_load(char(spm_BIDS(work_bids_root,'data',... 
    'modality','func','type','events','sub',subs{s},'run',runs{r}))); 
            clear conds 
        for t = 1:numel(trialtypes) 
                    conds.names{t} = trialtypes{t}; 
                    conds.durations{t} = 0; 
                    conds.onsets{t} = d.onset(strcmpi(d.stim_type,trialtypes{t}));  
        end 
            save(fullfile(outpth,...
                sprintf('sub-%s',subs{s}),'ses-mri','func',...
                sprintf('sub-%s_run-%s_spmdef.mat',subs{s},runs{r})),'-struct','conds'); 
    end 
end 