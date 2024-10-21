function [twsP_signal, f_twsP_signal] = aj_smooth_TWS_3D(ph, tissue_proba, param, pth_out)
% Function to apply tissue-weighted smoothing to 3D data (using spm_imcalc and spm_smooth)
%
% INPUT:
% ph           : 3D Data to smooth (Nifti file or matrix)
% tissue_proba : 4D tissue probability maps (N x N x N x nb_tissue)
% param        : Smoothing parameters (kernel size, etc.)
% pth_out      : Output path for the resulting files (optional)
%
% OUTPUT:
% twsP_signal  : Smoothed signal per tissue
% f_twsP_signal: Final signal after applying the explicit mask
%
% Draganski et al, 2011, doi:10.1016/j.neuroimage.2011.01.052
    
    %% Check input and prepare
    if nargin < 4, pth_out = ''; end
    fwhm = param.sm_kern_tws;  % Smoothing kernel width
    nb_tissue = size(tissue_proba, 4);  % Number of tissue classes
    [nx, ny, nz] = size(ph);  % Dimensions of the 3D data
    
    %% Initialize output
    twsP_signal = cell(nb_tissue, 1);
    
    % Create output paths for intermediate and final results
    if isempty(pth_out)
        pth_out = pwd;  % Default to current working directory if no output path provided
    end
    
    % Output filenames
    fn_out = fullfile(pth_out, 'final_smoothed_signal.nii');
    
    %% Check if `ph` is a file path or already a loaded matrix
    if ischar(ph) || isstring(ph)
        vol_ph = spm_vol(ph);  % Load NIfTI file header
        ph_data = spm_read_vols(vol_ph);  % Read the volume data
    else
        ph_data = ph;  % `ph` is already a 3D matrix
        vol_ph = [];  % No volume info needed if it's in memory
    end
    
    %% Loop over each tissue class and smooth using spm_smooth and spm_imcalc
    for ii = 1:nb_tissue
        % Create temporary filenames for intermediate results
        fn_weighted_tissue = fullfile(pth_out, sprintf('weighted_tissue_%d.nii', ii));
        fn_smooth_tissue = fullfile(pth_out, sprintf('smoothed_tissue_%d.nii', ii));
        
        % Prepare the tissue probability map
        tmp_proba = tissue_proba(:, :, :, ii);  % Extract tissue map
        
        % Check if the dimensions match between ph_data and tmp_proba
        if ~isequal(size(ph_data), size(tmp_proba))
            % Use spm_imcalc to resample tmp_proba to match ph_data
            vol_proba = spm_vol(tmp_proba);  % Volume info for tissue probability map
            
            % Write out the tissue probability map to a temp file and resample
            spm_imcalc({vol_proba, vol_ph}, fn_weighted_tissue, 'i1', struct('interp', 1)); % Use trilinear interpolation
            
            % Read the resampled tissue probability map
            tmp_proba = spm_read_vols(spm_vol(fn_weighted_tissue));
        end
        
        % Apply tissue weight to the data
        tmp_ph = ph_data .* tmp_proba;
        
        % Write the weighted tissue image to a temporary file if needed
        if ~isempty(vol_ph)
            spm_vol_fn = spm_vol(ph);  % Prepare volume structure for spm_imcalc
            spm_write_vol(spm_vol_fn, tmp_ph);
        end
        
        % Smooth the tissue-specific image using spm_smooth
        spm_smooth(fn_weighted_tissue, fn_smooth_tissue, fwhm);
        
        % Store the smoothed result
        twsP_signal{ii} = spm_read_vols(spm_vol(fn_smooth_tissue));
    end
    
    %% Combine tissue-weighted smoothed signals into the final output
    f_twsP_signal = zeros(nx, ny, nz);
    for ii = 1:nb_tissue
        f_twsP_signal = f_twsP_signal + twsP_signal{ii};
    end
    
    %% Save the final combined smoothed signal
    if ~isempty(vol_ph)
        spm_vol_out = vol_ph;  % Use the original data header for the final output
        spm_vol_out.fname = fn_out;  % Set the output filename
        spm_write_vol(spm_vol_out, f_twsP_signal);  % Write the smoothed output
    else
        % If no volume information, save data as .mat or output as needed
        save(fn_out, 'f_twsP_signal');
    end
end
