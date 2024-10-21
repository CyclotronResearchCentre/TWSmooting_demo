function [gsP_signal] = aj_smooth_gaussian(ph, param, dim)
% Function to apply standard Gaussian smoothing to data
% that can be in 1D, 2D, or 3D.
%
% INPUT:
% ph: Data to smooth
% param: Smoothing parameters (including kernel size)
% dim: Dimension of the data (1D, 2D or 3D)
%
% OUTPUT:
% gsP_signal: Smoothed signal

    % Generate Gaussian kernel
    wg = gausswin(param.sm_kern_gaussian);
    wg = wg / sum(wg);

    % Apply smoothing based on data dimensions
    switch dim
        case 1  % 1D Smoothing
            % INPUT DIMENSIONS : ph [1 x nb_pt]
            % OUTPUT DIMENSIONS : gsP_signal [nb_tissue x nb_pt]
            
            gsP_signal = filtfilt(wg, 1, ph);  % Apply 1D Gaussian smoothing
        
        case 2  % 2D Smoothing
            % INPUT DIMENSIONS : ph [nb_pt x nb_pt]
            % OUTPUT DIMENSIONS : gsP_signal [nb_pt x nb_pt]
            
            % Apply Gaussian smoothing on each dimension
            gsP_signal = ph;  % Initialize with the original data
            for i = 1:2
                gsP_signal = imgaussfilt(gsP_signal, param.sm_kern_gaussian);  % 2D Gaussian smoothing
            end

        case 3  % 3D Smoothing
            % INPUT DIMENSIONS : ph [nb_pt x nb_pt x nb_pt]
            % OUTPUT DIMENSIONS : gsP_signal [nb_pt x nb_pt x nb_pt]
            
            % Apply 3D Gaussian smoothing
            gsP_signal = imgaussfilt3(ph, param.sm_kern_gaussian);  % 3D Gaussian smoothing
            
            figure;
            imagesc(gsP_signal(:,:,64));
            axis image;
            colorbar;
            title('Gaussian Smoothing');
            xlabel('Position X');
            ylabel('Position Y');
            
            figure;
            slice(double(gsP_signal), size(gsP_signal,1)/2, size(gsP_signal,2)/2, size(gsP_signal,3)/2);
            title('Gaussian Smoothed Signal (gsP\_signal)');
            xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis');
            colorbar;
%             axis equal;
%             shading interp;

        otherwise
            error('Data dimension is neither 1D, 2D, nor 3D');
    end
end

%% Identifier la dimension des données
%     data_size = size(ph);
%     if isvector(ph)
%         num_dims = 1;
%     elseif ismatrix(ph)
%         num_dims = 2;
%     else
%         num_dims = length(data_size);
%         if num_dims ~= 3
%             error('Argument ph de mauvaises dimensions');
%         end
%     end
