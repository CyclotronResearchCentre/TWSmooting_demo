function [gs_imgaussfilt3_paths, gs_spm_paths,...
    tws_paths, smwTC_paths,...
    tspoon_paths] = ...
    aj_smoothing(MPM_paths, TCseg_paths, param, flag, dim)
%--------------------------------------------------------------------------
% Function to apply different smoothing techniques (Gaussian, TWS, TSPOON)
% to the given data, based on the flags provided.
%
% INPUTS
% con_info:     A structure array containing image volume information
%               (from spm_vol) of data to be smoothed
% mwTC_info:    A structure array containing image volume information
%               (from spm_vol) of modulated warped tissue classes
% param:        Smoothing parameters (from default)
% flag:         Flags to determine which smoothing techniques to apply
%               (from default)
% dim:          Dimension of the data (1D, 2D or 3D)
%
% OUTPUTS
% gsP_signal:       Gaussian smoothed signal from imgaussfilt3
% gs_path:          Gaussian smoothed file path from spm_smooth
% twsP_signal:      Tissue-weighted smoothed signal
% tspoon_paths:     Cell array of filenames (nifti files) of the Tissue-
%                   SPecific smOOthing compeNsated (TSPOON) for GM & WM
% gs_exMask_paths:  Cell array of filenames (nifti files) of the individual
%                   explicit mask
%--------------------------------------------------------------------------
% FUTURE DEV
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------
gs_imgaussfilt3_paths = 0;
gs_spm_paths = 0;
smwTC_paths = 0;
tws_paths = 0;
tspoon_paths = 0;

%% Apply standard Gaussian smoothing
if flag.gaussian
    disp('Executing standard Gaussian smoothing...');
    
    % Automatically create a new smoothing-specific derivatives directory for the subject
    pth_out = regexprep(spm_file(MPM_paths(1,:), 'fpath'), param.outDerivName, 'AJ-GS');
    if ~exist(pth_out, 'dir')
        mkdir(pth_out);
    else
        % If the directory already exists, delete all the files inside it
        file_list = dir(fullfile(pth_out, '*'));
        for k = 1:length(file_list)
            % Ignore the special directories '.' and '..'
            if ~file_list(k).isdir
                delete(fullfile(pth_out, file_list(k).name));
            end
        end
    end
    
    [gs_imgaussfilt3_paths, gs_spm_paths] = aj_smooth_gaussian(MPM_paths, pth_out, param.fwhm_gs, dim);
end

%% Apply tissue-weighted smoothing (TWS)
if flag.tws
    disp('Executing tissue-weighted smoothing (TWS)...');
    
    % Taking the tissue probability map paths from TPM.nii
    nTC = size(TCseg_paths,1);
    TPM_paths = cell(1,nTC);
    for i = 1:nTC
        TPM_paths{i} = fullfile(spm('Dir'), 'tpm', ['TPM.nii,' num2str(i)]);
    end
    
    % Automatically create a new smoothing-specific derivatives directory for the subject
    pth_out = regexprep(spm_file(MPM_paths(1,:), 'fpath'), param.outDerivName, 'AJ-TWS');
    if ~exist(pth_out, 'dir')
        mkdir(pth_out);
    else
        % If the directory already exists, delete all the files inside it
        file_list = dir(fullfile(pth_out, '*'));
        for k = 1:length(file_list)
            % Ignore the special directories '.' and '..'
            if ~file_list(k).isdir
                delete(fullfile(pth_out, file_list(k).name));
            end
        end
    end
    
    [tws_paths, smwTC_paths] = hmri_proc_MPMsmooth(char(MPM_paths),...
        char(TCseg_paths), char(TPM_paths), param.fwhm_tws,param.l_TC, pth_out);
end

%% Apply tissue-specific smoothing compensated (TSPOON)
if flag.tspoon
    disp('Executing tissue-specific smoothing compensated (TSPOON)...');
    
    % Automatically create a new smoothing-specific derivatives directory for the subject
    pth_out = regexprep(spm_file(MPM_paths(1,:), 'fpath'), param.outDerivName, 'AJ-TSPOON');
    if ~exist(pth_out, 'dir')
        mkdir(pth_out);
    else
        % If the directory already exists, delete all the files inside it
        file_list = dir(fullfile(pth_out, '*'));
        for k = 1:length(file_list)
            % Ignore the special directories '.' and '..'
            if ~file_list(k).isdir
                delete(fullfile(pth_out, file_list(k).name));
            end
        end
    end
    
    [tspoon_paths, smwTC_paths] = aj_proc_MPMTPSOON(char(MPM_paths), ...
        char(TCseg_paths), param.fwhm_tspoon,param.l_TC, pth_out);
end

%% Plot results
if flag.plot_fig && dim == 3 && false
    for i = 1:1 %length(con_info) % only for con0002
        con = spm_read_vols(MPM_paths{i});
        
        [~, name, ~] = fileparts(MPM_paths{i}.fname);
        figure('Name', sprintf('Original image for file: %s', name));
        
        subplot(1, 2, 1);
        imagesc(con(:, :, round(size(con,3)/2)));
        axis image;
        colorbar;
        title('2D view');
        xlabel('Position X'); ylabel('Position Y');

        subplot(1, 2, 2);
        slice(double(con), size(con, 1) / 2, size(con, 2) / 2, size(con, 3) / 2);
        title('3D view');
        xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis');
        colorbar;
    end
    
    if flag.gaussian
        for i = 1:1 %length(gs_imgaussfilt3_paths) % only for con0002
            gs_imgaussfilt3_info = spm_vol(gs_imgaussfilt3_paths{i});
            gs = spm_read_vols(gs_imgaussfilt3_info);
            
            [~, name, ~] = fileparts(gs_imgaussfilt3_paths{i});
            figure('Name', sprintf('GS for file: %s', name));
            
            subplot(1, 2, 1);
            imagesc(gs(:, :, round(size(gs,3)/2)));
            axis image;
            colorbar;
            title('2D view');
            xlabel('Position X'); ylabel('Position Y');

            subplot(1, 2, 2);
            slice(double(gs), size(gs,1)/2, size(gs,2)/2, size(gs,3)/2);
            title('3D view');
            xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis');
            colorbar;
        end
        
        for i = 1:1 %length(gs_spm_paths) % only for con0002
            gs = spm_read_vols(spm_vol(gs_spm_paths{i}));
            
            name = spm_file(gs_spm_paths{i}, 'filename');
            figure('Name', sprintf('GS for file: %s', name));
            
            subplot(1, 2, 1);
            imagesc(gs(:, :, round(size(gs,3)/2)));
            axis image;
            colorbar;
            title('2D view');
            xlabel('Position X'); ylabel('Position Y');

            subplot(1, 2, 2);
            slice(double(gs), size(gs,1)/2, size(gs,2)/2, size(gs,3)/2);
            title('3D view');
            xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis');
            colorbar;
        end
    end
    
    if flag.tws
        for i = 1:1 % size(tws_paths,2) % only for con0002
            twsGM = spm_read_vols(spm_vol(tws_paths{1}(1,:)));
            twsWM = spm_read_vols(spm_vol(tws_paths{1}(2,:)));
            
            name = [spm_file(tws_paths{1}(1,:), 'filename'), ' and ', spm_file(tws_paths{1}(2,:), 'filename')];
            figure('Name', sprintf('TWS for file: %s', name));
            
            subplot(2, 2, 1);
            imagesc(twsGM(:, :, round(size(twsGM, 3) / 2)));
            axis image;
            colorbar;
            title('2D view');
            xlabel('Position X'); ylabel('Position Y');
            
            subplot(2, 2, 2);
            imagesc(twsWM(:, :, round(size(twsWM, 3) / 2)));
            axis image;
            colorbar;
            title('2D view');
            xlabel('Position X'); ylabel('Position Y');
            
            subplot(2, 2, 3);
            slice(double(twsGM), size(twsGM,1)/2, size(twsGM,2)/2, size(twsGM,3)/2);
            title('3D view');
            xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis');
            colorbar;
            
            subplot(2, 2, 4);
            slice(double(twsWM), size(twsWM,1)/2, size(twsWM,2)/2, size(twsWM,3)/2);
            title('3D view');
            xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis');
            colorbar;
        end
    end
    
    if flag.tspoon
        for i = 1:1 % size(tspoon_paths,1) % only for con0002
            tspoonGM = spm_read_vols(spm_vol(tspoon_paths{1}(1,:)));
            tspoonWM = spm_read_vols(spm_vol(tspoon_paths{1}(2,:)));
            
            name = [spm_file(tspoon_paths{1}(1,:), 'filename'), ' and ', spm_file(tspoon_paths{1}(2,:), 'filename')];
            figure('Name', sprintf('TSPOON for file: %s', name));
            
            subplot(2, 2, 1);
            imagesc(tspoonGM(:, :, round(size(tspoonGM, 3) / 2)));
            axis image;
            colorbar;
            title('2D view');
            xlabel('Position X'); ylabel('Position Y');
            
            subplot(2, 2, 2);
            imagesc(tspoonWM(:, :, round(size(tspoonWM, 3) / 2)));
            axis image;
            colorbar;
            title('2D view');
            xlabel('Position X'); ylabel('Position Y');
            
            subplot(2, 2, 3);
            slice(double(tspoonGM), size(tspoonGM,1)/2, size(tspoonGM,2)/2, size(tspoonGM,3)/2);
            title('3D view');
            xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis');
            colorbar;
            
            subplot(2, 2, 4);
            slice(double(tspoonWM), size(tspoonWM,1)/2, size(tspoonWM,2)/2, size(tspoonWM,3)/2);
            title('3D view');
            xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis');
            colorbar;
        end
    end
end

end
