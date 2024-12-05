% Draft script to better understand TWS and TSPOON smoothing
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------
%% Profile of TWS and TSPOON in 1D for one subject and combinaison
clear;clc;
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');

out_dir = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\AJ-Diff_TWSTSPOON';

ds_dir = 'D:\Master_Thesis\Data\BIDS_AgingData';
smoothing_methods = {'TWS', 'TSPOON'};
TWS_dir = fullfile(ds_dir,'derivatives',sprintf('AJ-%s',smoothing_methods{1}));
TSPOON_dir = fullfile(ds_dir,'derivatives',sprintf('AJ-%s',smoothing_methods{2}));

qMRI_params = {'MTsat', 'PDmap', 'R1map', 'R2starmap'};
TCs = {'GM', 'WM'};

% Combining parameters (Cartesian product)
[TCs_idx, qMRI_idx] = ndgrid(1:numel(TCs), 1:numel(qMRI_params));

% Creating TWS and TSPOON patterns
TWS_patterns = arrayfun(@(i, j) sprintf('^%s_%s.*%s\\.nii$', ...
                        smoothing_methods{1}, TCs{i}, qMRI_params{j}), ...
                        TCs_idx(:), qMRI_idx(:), 'UniformOutput', false);
TSPOON_patterns = arrayfun(@(i, j) sprintf('^%s_%s.*%s\\.nii$', ...
                        smoothing_methods{2}, TCs{i}, qMRI_params{j}), ...
                        TCs_idx(:), qMRI_idx(:), 'UniformOutput', false);
% Getting the list of TWS and TSPOON paths 
TWS_paths = cellfun(@(pattern) spm_select('FPListRec', TWS_dir, pattern), ...
                    TWS_patterns, 'UniformOutput', false);
TSPOON_paths = cellfun(@(pattern) spm_select('FPListRec', TSPOON_dir, pattern), ...
                       TSPOON_patterns, 'UniformOutput', false);
TWS_GM_MTsat_path = TWS_paths{1}(1,:);
TSPOON_GM_MTsat_path = TSPOON_paths{1}(1,:);
GMseg_path = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\AJ-TWS\sub-001\anat\gs_sub-S001_MTsat_space-MNI_desc-mod_label-GM_probseg.nii';
WMseg_path = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\AJ-TWS\sub-001\anat\gs_sub-S001_MTsat_space-MNI_desc-mod_label-WM_probseg.nii';

TWS_GM_MTsat = spm_read_vols(spm_vol(TWS_GM_MTsat_path));
TSPOON_GM_MTsat = spm_read_vols(spm_vol(TSPOON_GM_MTsat_path));
GMseg = spm_read_vols(spm_vol(GMseg_path));
WMseg = spm_read_vols(spm_vol(WMseg_path));

TWS_middleLine = TWS_GM_MTsat(90,110,:);
TSPOON_middleLine = TSPOON_GM_MTsat(90,110,:);
GMseg_middleLine = GMseg(90,110,:);
WMseg_middleLine = WMseg(90,110,:);

% Convertir les NaN en 0 pour chaque vecteur
TWS_middleLine(isnan(TWS_middleLine)) = 0;
TSPOON_middleLine(isnan(TSPOON_middleLine)) = 0;
GMseg_middleLine(isnan(GMseg_middleLine)) = 0;
WMseg_middleLine(isnan(WMseg_middleLine)) = 0;

% Extraire les vecteurs comme des lignes 1D
TWS_middleLine_1D = squeeze(TWS_middleLine);
TSPOON_middleLine_1D = squeeze(TSPOON_middleLine);
GMseg_middleLine_1D = squeeze(GMseg_middleLine);
WMseg_middleLine_1D = squeeze(WMseg_middleLine);

% Création du graphe
figure; % Ouvre une nouvelle figure
hold on; % Permet de superposer plusieurs tracés

% Tracer chaque vecteur
plot(TWS_middleLine_1D, '-r', 'LineWidth', 1.5, 'DisplayName', 'TWS GM MTsat');
plot(TSPOON_middleLine_1D, '-b', 'LineWidth', 1.5, 'DisplayName', 'TSPOON GM MTsat');
plot(GMseg_middleLine_1D, '-g', 'LineWidth', 1.5, 'DisplayName', 'GM segmentation');
plot(WMseg_middleLine_1D, '-k', 'LineWidth', 1.5, 'DisplayName', 'WM segmentation');

% Ajouter des légendes, titres et axes
legend('Location', 'best'); % Légende positionnée automatiquement
xlabel('Voxel Index');
ylabel('Intensity');
title('Middle Line Comparison');
grid on; % Ajouter une grille pour une meilleure visualisation

hold off; % Fin de la superposition

%% Compute mean of smoothed TCprobseg for all TWS subjects
dir = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\AJ-TWS';
GMprobseg_paths = spm_select('FPListRec',dir,'^gs_sub-.*CSF_probseg.nii$');
data = zeros(181,217,181,138);
for i = 1:138
    data(:,:,:,i) = spm_read_vols(spm_vol(GMprobseg_paths(i,:)));
end
mean_data = mean(data, 4);

% Save the mean as a nifti file
mean_vol = spm_vol(GMprobseg_paths(1,:));
mean_vol.fname = fullfile(dir, 'mean_sCSFprobseg.nii');
spm_write_vol(mean_vol, mean_data);

%% Histogramm of diff_GM vs. smwc1
diff_path = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\AJ-Diff_TWSTSPOON_GLM\diff_GM_MTsat\spmF_0001.nii';
dir = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\AJ-TWS';
GMseg_path = fullfile(dir, 'mean_GMprobseg.nii');

diff = spm_read_vols(spm_vol(diff_path));
GMseg = spm_read_vols(spm_vol(GMseg_path));

valid_idx = ~isnan(diff) & diff ~= 0 & ~isnan(GMseg) & GMseg ~= 0;
GMseg_values = GMseg(valid_idx);
reverse_GMseg_values = 1./GMseg_values;
diff_values = diff(valid_idx);

figure;
scatter(reverse_GMseg_values, diff_values, 5, 'filled');
xlabel('GMseg Values ^{-1}');
ylabel('Diff Values');
title('Scatter Plot of GMseg vs Diff');
grid on;

%% In Silico Simulation
close all;clear;clc;

GM_d = [1 1 1 1 1 1 1 1 0.75 0.75 0.6 0.5 0.25 0 0 0 0 0 0 0 0 0 0 0];

n = length(GM_d);
TC_d = zeros(3,n);
TC_d(1,:) = GM_d;
TC_d(2,:) = ones(1,n) - GM_d; % WM_d
TC_d(3,:) = zeros(1,n); % CSF_d

GM_s = zeros(1,n);
GM_s(:) = 20;
WM_s = zeros(1,n);
WM_s(:) = 10;

signal = TC_d(1,:).*GM_s + TC_d(2,:).*WM_s;

% figure;
% hold on;
% plot(signal, '-r', 'LineWidth', 1.5, 'DisplayName', 'signal');
% legend('Location', 'best');
% xlabel('Voxel Index');
% ylabel('Intensity');
% title('Signal');
% ylim([0 25]);
% grid on;
% hold off;

% Apply standard smoothing
% Y = filtfilt(B, A, X): The length of the input X must be more than three
% times the filter order, defined as max(length(B)-1,length(A)-1).
sm_kern = 4;
wg = gausswin(sm_kern);
wg = wg/sum(wg); % normalize
wg2 = gausswin(2*sm_kern); % double width smoothing
wg2 = wg2/sum(wg2); % normalize

gs_signal = filtfilt(wg,1,signal);
gs_TC_d = filtfilt(wg,1,TC_d')';

% Apply tissue-weighted smoothing
% Appplying the smoothing as implemented for VBQ,
% assuming the TPMs are like the tissue probability but smoothed with a
% kernel twice the size, for simplicity.
TWS_signal = zeros(2,n);
for ii=1:2
    tmp1 = signal .* TC_d(ii,:) .* (filtfilt(wg2,1,TC_d(ii,:))>.05); % Like the TPM masking
    tmp2 = gs_TC_d(ii,:) .* (gs_TC_d(ii,:)>.05); % masking from smoothed tissue
    TWS_signal(ii,:) = filtfilt(wg,1,tmp1) ./ tmp2;
end

% figure;
% hold on;
% plot(signal, '-r', 'LineWidth', 1.5, 'DisplayName', 'signal');
% plot(TWS_signal(1,:), '-b', 'LineWidth', 1.5, 'DisplayName', 'GM TWS signal');
% plot(TWS_signal(2,:), '-g', 'LineWidth', 1.5, 'DisplayName', 'WM TWS signal');
% legend('Location', 'best');
% xlabel('Voxel Index');
% ylabel('Intensity');
% title('TWS');
% ylim([0 25]);
% grid on;
% hold off;

% Explicit mask
% majority and above 20%
exMask = [  ...
    TC_d(1,:)>TC_d(2,:) & ... % GM>WM
    TC_d(1,:)>TC_d(3,:) & ... % GM>CSF
    TC_d(1,:)>.2 ; ...                % GM>.2
    TC_d(2,:)>TC_d(1,:) & ... % WM>GM
    TC_d(2,:)>TC_d(3,:) & ... % WM>CSF
    TC_d(2,:)>.2 ] ;                  % WM>.2

% figure;
% hold on;
% plot(signal, '-r', 'LineWidth', 1.5, 'DisplayName', 'signal');
% plot(TWS_signal(1,:).*exMask(1,:), '-b', 'LineWidth', 1.5, 'DisplayName', 'GM_TWS_signal*GM_exMask');
% plot(TWS_signal(2,:).*exMask(2,:), '-g', 'LineWidth', 1.5, 'DisplayName', 'WM_TWS_signal*WM_exMask');
% legend('Location', 'best');
% xlabel('Voxel Index');
% ylabel('Intensity');
% title('TWS*exMask');
% ylim([0 25]);
% grid on;
% hold off;

% Apply TSPOON smoothing
gs_exMask = filtfilt(wg,1,double(exMask)')';
TSPOON_signal = zeros(2,n);
for ii=1:2
    TSPOON_signal(ii,:) = filtfilt(wg,1,exMask(ii,:).*signal) ./ gs_exMask(ii,:);
end

% figure;
% hold on;
% plot(signal, '-r', 'LineWidth', 1.5, 'DisplayName', 'signal');
% plot(TSPOON_signal(1,:), '-b', 'LineWidth', 1.5, 'DisplayName', 'GM TSPOON signal');
% plot(TSPOON_signal(2,:), '-g', 'LineWidth', 1.5, 'DisplayName', 'WM TSPOON signal');
% legend('Location', 'best');
% xlabel('Voxel Index');
% ylabel('Intensity');
% title('TSPOON');
% ylim([0 25]);
% grid on;
% hold off;

figure;
hold on;
plot(signal, '-k', 'LineWidth', 1.5, 'DisplayName', 'signal');

plot(TWS_signal(1,:), '-b', 'LineWidth', 1.5, 'DisplayName', 'GM TWS signal');
plot(TWS_signal(2,:), '-c', 'LineWidth', 1.5, 'DisplayName', 'WM TWS signal');

plot(TSPOON_signal(1,:), '-r', 'LineWidth', 1.5, 'DisplayName', 'GM TSPOON signal');
plot(TSPOON_signal(2,:), '-y', 'LineWidth', 1.5, 'DisplayName', 'WM TSPOON signal');
legend('Location', 'best');
xlabel('Voxel Index');
ylabel('Intensity');
title('TWS + TSPOON');
ylim([0 25]);
grid on;
hold off;

