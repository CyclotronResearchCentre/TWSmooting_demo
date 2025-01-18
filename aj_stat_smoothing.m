% Toy script to better understand difference and normalized difference
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------
%% Quick tests
clear;clc;
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');

diff_path = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\AJ-Diff_TWSTSPOON\sub-001\diff_GM_sub-S001_space-MNI_MTsat.nii';
normdiff_path = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\AJ-Diff_normTWSTSPOON\sub-001\diff_GM_sub-S001_space-MNI_MTsat.nii';

diff = spm_read_vols(spm_vol(diff_path));
normdiff = spm_read_vols(spm_vol(normdiff_path));

flag.drawPlot = 1;
flag.savePlot = 0;
[mean_diff,std_diff] = aj_BlandAltman(diff, normdiff, flag);
