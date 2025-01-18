% Toy script to test Bland-Altman Plot on fMRI smoothed contrast maps.
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------
%% Bland-Atlman for smoothed ds000117
close all; clear;clc;
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\agingdata\reprod_article');

gs_path = 'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\work_copies\openneuro.org\ds000117\derivatives\AJ-GS\sub-01\ses-mri\func\spm_smooth_con_0002.nii';
tws_GM_path = 'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\work_copies\openneuro.org\ds000117\derivatives\AJ-TWS\sub-01\ses-mri\func\TWS_GMw_con_0002.nii';
tspoon_GM_path = 'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\work_copies\openneuro.org\ds000117\derivatives\AJ-TSPOON\sub-01\ses-mri\func\TSPOON_GM_con_0002.nii';

gs = spm_read_vols(spm_vol(gs_path));
tws_GM = spm_read_vols(spm_vol(tws_GM_path));
tspoon_GM = spm_read_vols(spm_vol(tspoon_GM_path));

flag.drawPlot = 1;
flag.savePlot = 1;
pth_out = 'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\work_copies\openneuro.org\ds000117\derivatives\AJ-StatData\BlandAltman\TWS-TSPOON_GM_con0002.png';
plot_title = 'TWS-TSPOON GM';
[all_mean_diff,all_std_diff] = aj_BlandAltman(tws_GM, tspoon_GM, flag, pth_out, plot_title);
