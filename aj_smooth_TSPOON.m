function [tosP_signal, f_tosP_signal] = aj_smooth_TSPOON(ph, tissue_proba, param, dim)
% Applies TSPOON (Tissue-SPecific smOOthing compeNsated) on 1D, 2D or 3D data.
%
% INPUT:
% ph: Data to smooth
% tissue_proba: Tissue probability data
% param: Smoothing parameters (including kernel size)
% dim: Dimension of the data (1D, 2D or 3D)
%
% OUTPUT:
% tosP_signal: Smoothed signal per tissue
% f_tosP_signal: Final signal after applying the explicit mask

    % Create Gaussian kernel
    wg = gausswin(param.sm_kern_tspoon);
    wg = wg / sum(wg);

    % Apply TSPOON smoothing based on the data dimension
    switch dim
        case 1  % 1D Smoothing
            % INPUT DIMENSIONS : ph [1 x nb_pt] and tissue_proba [nb_tissue x nb_pt]
            % OUTPUT DIMENSIONS : tosP_signal [nb_tissue x nb_pt] and f_tosP_signal [1 x nb_pt]
            
            % Ensure ph and tissue_proba have compatible dimensions
            ph = ph'; % now [nb_pt x 1]
            tissue_proba = tissue_proba'; % now [nb_pt x nb_tissue]
            
            % Generate the explicit mask using tissue probabilities and standard Gaussian smoothing
            gsP_GmWmCsfSculpt = filtfilt_gauss_proba(tissue_proba, wg, dim);
            iexMask = aj_create_explicit_mask(gsP_GmWmCsfSculpt, 1);

            % Smooth the explicit mask
            gsP_iexMask = filtfilt_gauss_mask(double(iexMask), wg, dim);

            % Initialize tosP_signal
            tosP_signal = zeros(size(iexMask));

            for ii = 1:size(iexMask, 1)
                tosP_signal(ii, :) = filtfilt(wg, 1, iexMask(ii, :) .* ph) ./ gsP_iexMask(ii, :);
            end
            
            tosP_signal(isnan(tosP_signal) | isinf(tosP_signal)) = 0;  % Handle NaN/Inf

            % Combine tissue signals into the final smoothed signal using the explicit mask
            f_tosP_signal = sum(tosP_signal .* iexMask, 1);  % Element-wise multiplication and summation across tissues
            f_tosP_signal(isnan(f_tosP_signal) | isinf(f_tosP_signal)) = 0;  % Handle NaN/Inf in the final signal
            
        case 2 % 2D Smoothing
            % Ensure dimensions match: ph should be [1 x nb_pt_x x nb_pt_y]
            ph = reshape(ph, [1, size(ph)]);  % Reshape to ensure 3D for consistency

            % Generate the explicit mask
            gsP_GmWmCsfSculpt = filtfilt_gauss_proba(tissue_proba, wg, dim);
            iexMask = aj_create_explicit_mask(gsP_GmWmCsfSculpt, 2);

            % Smooth the explicit mask
            gsP_iexMask = filtfilt_gauss_mask(double(iexMask), wg, dim);

            % Initialize tosP_signal
            tosP_signal = zeros(size(iexMask, 1), size(iexMask, 2), size(iexMask, 3));

            for ii = 1:size(iexMask, 1)  % Iterate over each tissue type
                tosP_signal(ii, :, :) = filtfilt_gauss_mask(iexMask(ii, :, :) .* ph, wg, dim) ./ gsP_iexMask(ii, :, :);
            end
            
            % Handle NaN/Inf values
            tosP_signal(isnan(tosP_signal) | isinf(tosP_signal)) = 0;  
            % Combine tissue signals into the final smoothed signal
            f_tosP_signal = sum(tosP_signal .* iexMask, 1);  % Summing over the first dimension

        case 3 % 3D Smoothing
            % Ensure dimensions match: ph should be [1 x nb_pt_x x nb_pt_y x nb_pt_z]
            ph = reshape(ph, [1, size(ph)]);  % Reshape to ensure 4D for consistency

            disp(size(tissue_proba));
            % Generate the explicit mask
            gsP_GmWmCsfSculpt = filtfilt_gauss_proba(tissue_proba, wg, dim);
            iexMask = aj_create_explicit_mask(gsP_GmWmCsfSculpt, 3);
            disp(size(gsP_GmWmCsfSculpt));
            disp(size(iexMask));

            % Smooth the explicit mask
            gsP_iexMask = filtfilt_gauss_mask(double(iexMask), wg, dim);
            disp(size(gsP_iexMask));

            % Initialize tosP_signal
            tosP_signal = zeros(size(iexMask));

            for ii = 1:size(iexMask, 1)  % Iterate over each tissue type
                tosP_signal(:, :, :, ii) = filtfilt_gauss(iexMask(:, :, :, ii) .* ph, wg, dim) ./ gsP_iexMask(:, :, :, ii);
            end

            % Handle NaN/Inf values
            tosP_signal(isnan(tosP_signal) | isinf(tosP_signal)) = 0;  
            % Combine tissue signals into the final smoothed signal
            f_tosP_signal = sum(tosP_signal .* iexMask, 4);  % Summing over the first dimension

        otherwise
            error('Data dimension is neither 1D, 2D, nor 3D');
    end
end

%% Function to apply filtering based on dimension (1D, 2D, 3D)
function gsP_iexMask = filtfilt_gauss_proba(data, wg, num_dims)
    if num_dims == 1  % 1D
        gsP_iexMask = filtfilt(wg, 1, data)';
    elseif num_dims == 2  % 2D
        sigma = std(wg);  % Using standard deviation of the Gaussian kernel
        gsP_iexMask = imgaussfilt(data, sigma);  % 2D Gaussian smoothing
    elseif num_dims == 3  % 3D
        gsP_iexMask = smooth3(data, 'gaussian', [wg wg wg]);  % 3D Gaussian smoothing
    else
        error('Unsupported dimensions for filtfilt_gauss_proba.');
    end
end

%% Function to apply filtering based on dimension (1D, 2D, 3D)
function gsP_iexMask = filtfilt_gauss_mask(data, wg, num_dims)
    if num_dims == 1  % 1D
        gsP_iexMask = filtfilt(wg, 1, data')';
    elseif num_dims == 2  % 2D
        sigma = std(wg);  % Using standard deviation of the Gaussian kernel
        gsP_iexMask = imgaussfilt(data, sigma);  % 2D Gaussian smoothing
    elseif num_dims == 3  % 3D
        gsP_iexMask = smooth3(data, 'gaussian', [wg wg wg]);  % 3D Gaussian smoothing
    else
        error('Unsupported dimensions for filtfilt_gauss_mask.');
    end
end
