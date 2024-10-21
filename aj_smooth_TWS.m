function [twsP_signal, f_twsP_signal] = aj_smooth_TWS(ph, tissue_proba, param, dim)
    % Function to apply tissue-weighted smoothing to data that can be in 1D, 2D or 3D
    % (with tissue probabilities: GM, WM, CSF, Sculpt).
    %
    % INPUT:
    % ph: Data to smooth
    % tissue_proba: Tissue probabilities matrix
    % param: Smoothing parameters (including kernel size)
    % dim: Dimension of the data (1D, 2D or 3D)
    %
    % OUTPUT:
    % twsP_signal: Smoothed signal per tissue
    % final_signal: Final signal after applying the explicit mask
    %
    % Draganski et al, 2011, doi:10.1016/j.neuroimage.2011.01.052
   
    % Generate Gaussian kernels
    wg = gausswin(param.sm_kern_tws);
    wg = wg / sum(wg);
    wg2 = gausswin(2 * param.sm_kern_tws);
    wg2 = wg2 / sum(wg2);

    nb_tissue = min(size(tissue_proba));
    nb_pt = max(size(tissue_proba));  % same number as max(size(ph))

    switch dim
        case 1 % 1D Smoothing
            % INPUT DIMENSIONS : ph [1 x nb_pt] and tissue_proba [nb_tissue x nb_pt]
            % OUTPUT DIMENSIONS : twsP_signal [nb_tissue x nb_pt] and f_twsP_signal [1 x nb_pt]
            
            % Initialize the variable for the smoothed signal
            twsP_signal = zeros([nb_tissue, nb_pt]);

            % Create the explicit mask using the function aj_create_explicit_mask
            iexMask = aj_create_explicit_mask(tissue_proba, 1);
            
            % Ensure ph and tissue_proba have compatible dimensions 
            tissue_proba = tissue_proba'; % tissue_proba is now [nb_pt x nb_tissue]
            
            % Apply smoothing for each tissue
            for ii = 1:nb_tissue
                filt_data2 = filtfilt(wg2, 1, tissue_proba(:, ii));  % Smoothing on the corresponding tissue probability
                tissue_mask2 = (filt_data2 > 0.05);  % Create binary mask
                tmp1 = ph .* tissue_proba(:, ii) .* tissue_mask2;

                % CP
%                 filt_data1 = filtfilt(wg, 1, tissue_proba(:, ii));
%                 tissue_mask1 = (filt_data1 > 0.05);
                
                % V1
%                 tissue_mask1 = filtfilt(wg, 1, tissue_proba(:, ii));
%                 tissue_mask1(tissue_mask1 <= 0) = NaN;  % Prevent division by zero or negative values

                % V2
                tissue_mask1 = filtfilt(wg, 1, tissue_proba(:, ii));
                tissue_mask1(tissue_mask1 <= 0.05) = NaN;

                smoothed_signal = filtfilt(wg, 1, tmp1) ./ tissue_mask1;
                smoothed_signal(isnan(smoothed_signal) | isinf(smoothed_signal)) = 0;  % Handle NaN/Inf
                twsP_signal(ii, :) = smoothed_signal';
            end
            
            % Apply the explicit mask to combine tissue signals into the final signal
            f_twsP_signal = sum(twsP_signal .* iexMask, 1);  % Element-wise multiplication and summation across tissues
            f_twsP_signal(isnan(f_twsP_signal) | isinf(f_twsP_signal)) = 0;  % Handle NaN/Inf in the final signal
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 2 % 2D Smoothing
            % INPUT DIMENSIONS : ph [nb_pt x nb_pt] and tissue_proba [nb_tissue x nb_pt x nb_pt]
            % OUTPUT DIMENSIONS : twsP_signal [nb_tissue x nb_pt x nb_pt] and f_twsP_signal [nb_pt x nb_pt]
            
            % Initialize the variable for the smoothed signal
            twsP_signal = zeros([nb_tissue, nb_pt, nb_pt]);

            % Create the explicit mask using the function aj_create_explicit_mask
            iexMask = aj_create_explicit_mask(tissue_proba, 2);
            
            % Ensure tissue_proba has dimensions [N x N x nb_tissue]
            dims = size(tissue_proba);
            [nb_tissue, idx_min] = min(dims);
            permute_order = [setdiff(1:3, idx_min), idx_min];
            if idx_min ~= 3
                tissue_proba = permute(tissue_proba, permute_order);
            end
            
            for ii = 1:nb_tissue
                tissue_mask2 = imgaussfilt(tissue_proba(:,:,ii), 2*param.sm_kern_tws) > 0.05;
                tmp1 = ph .* tissue_proba(:,:,ii) .* tissue_mask2;
                
                % CP
                tissue_mask1 = imgaussfilt(tissue_proba(:,:,ii), param.sm_kern_tws) > 0.05;

                % V1 - Handle NaN/Inf after filtering
%                 tissue_mask1 = imgaussfilt(tissue_proba(:,:,ii), param.sm_kern_tws);
%                 tissue_mask1(tissue_mask1 <= 0) = NaN;
                
                % V2
%                 tissue_mask1 = imgaussfilt(tissue_proba(:,:,ii), param.sm_kern_tws);
%                 tissue_mask1(tissue_mask1 <= 0.05) = NaN;

                smoothed_signal = imgaussfilt(tmp1, param.sm_kern_tws) ./ tissue_mask1;
                smoothed_signal(isnan(smoothed_signal) | isinf(smoothed_signal)) = 0;
                twsP_signal(ii, :, :) = smoothed_signal;
            end

            % Apply the explicit mask to combine tissue signals into the final signal
            f_twsP_signal = zeros(nb_pt, nb_pt);
            for ii = 1:nb_tissue
                % Combine the smoothed signals across tissues using the explicit mask
                f_twsP_signal = f_twsP_signal + squeeze(twsP_signal(ii, :, :)) .* squeeze(iexMask(ii, :, :));
            end

            % Handle NaN/Inf in the final signal
            f_twsP_signal(isnan(f_twsP_signal) | isinf(f_twsP_signal)) = 0;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 3 % 3D Smoothing
            % INPUT DIMENSIONS : ph [nb_pt x nb_pt x nb_pt] and tissue_proba [nb_tissue x N x N x N]
            % OUTPUT DIMENSIONS : twsP_signal [nb_tissue x N x N x N] and f_twsP_signal [nb_pt x nb_pt x nb_pt]
            
            % Initialize the variable for the smoothed signal
            twsP_signal = zeros([nb_tissue, nb_pt, nb_pt, nb_pt]);

            % Create the explicit mask using the function aj_create_explicit_mask_1D
            iexMask = aj_create_explicit_mask(tissue_proba, dim);
                        
            % Ensure tissue_proba has dimensions [N x N x N x nb_tissue]
            dims = size(tissue_proba);
            [nb_tissue, idx_min] = min(dims);
            permute_order = [setdiff(1:4, idx_min), idx_min];
            if idx_min ~= 4
                tissue_proba = permute(tissue_proba, permute_order);
            end
            
            for ii = 1:nb_tissue
                tissue_mask2 = imgaussfilt3(tissue_proba(:,:,:,ii), 2*param.sm_kern_tws) > 0.05;
                tmp1 = ph .* tissue_proba(:,:,:,ii) .* tissue_mask2;

                % CP
                tissue_mask1 = imgaussfilt(tissue_proba(:,:,:,ii), param.sm_kern_tws) > 0.05;

                % V1 - Handle NaN/Inf after filtering
%                 tissue_mask1 = imgaussfilt3(tissue_proba(:,:,:,ii), param.sm_kern_tws);
%                 tissue_mask1(tissue_mask1 <= 0) = NaN;
                
                % V2
%                 tissue_mask1 = imgaussfilt(tissue_proba(:,:,:,ii), param.sm_kern_tws);
%                 tissue_mask1(tissue_mask1 <= 0.05) = NaN;

                smoothed_signal = imgaussfilt3(tmp1, param.sm_kern_tws) ./ tissue_mask1;
                smoothed_signal(isnan(smoothed_signal) | isinf(smoothed_signal)) = 0;
                
                twsP_signal(ii, :, :, :) = smoothed_signal;
            end
            
            % Apply the explicit mask to combine tissue signals into the final signal
            f_twsP_signal = zeros(nb_pt, nb_pt, nb_pt);
            for ii = 1:nb_tissue
                % Combine the smoothed signals across tissues using the explicit mask
                f_twsP_signal = f_twsP_signal + squeeze(twsP_signal(ii, :, :, :)) .* squeeze(iexMask(ii, :, :, :));
            end

            % Handle NaN/Inf in the final signal
            f_twsP_signal(isnan(f_twsP_signal) | isinf(f_twsP_signal)) = 0;

        otherwise
            error('Data dimension is neither 1D, 2D nor 3D');
    end
end