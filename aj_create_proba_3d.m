function [data_GmWmCsfSculpt_3D] = aj_create_proba_3d(phantom_var, param, plot_fig)
    % AJ_CREATE_PROBA_3D Calculate tissue probability for 3D image (volume)
    %
    % This function takes a 3D volume of intensity values (data_3D)
    % and calculates the probability of each voxel intensity corresponding to
    % four types of tissue: Bone Sculpt, CSF, WM, and GM. The probabilities 
    % are computed using Gaussian distributions.
    
    [rows, cols, slices] = size(phantom_var);
    
    % Initialize the output matrix for storing probabilities (Bone, CSF, WM, GM)
    data_GmWmCsfSculpt_3D = zeros(rows, cols, slices, 4);

    % Calculate probabilities for each voxel in the 3D volume
    for i = 1:rows
        for j = 1:cols
            for k = 1:slices
                intensity = phantom_var(i, j, k);

                % Calculate Gaussian probabilities
                prob_gm = normpdf(intensity, param.ref_gm, param.sigma_gm);
                prob_wm = normpdf(intensity, param.ref_wm, param.sigma_wm);
                prob_csf = normpdf(intensity, param.ref_csf, param.sigma_csf);
                prob_bone = normpdf(intensity, param.ref_bone, param.sigma_bone);

                % Sum of probabilities
                total_prob = prob_gm + prob_wm + prob_csf + prob_bone;

                % Normalize to get probabilities
                data_GmWmCsfSculpt_3D(i, j, k, 1) = prob_gm / total_prob;   % GM
                data_GmWmCsfSculpt_3D(i, j, k, 2) = prob_wm / total_prob;   % WM
                data_GmWmCsfSculpt_3D(i, j, k, 3) = prob_csf / total_prob;  % CSF
                data_GmWmCsfSculpt_3D(i, j, k, 4) = prob_bone / total_prob; % Bone
            end
        end
    end

    % Visualization of the 3D volume
    if plot_fig
        figure;
        num_slices_to_show = length(param.slice_indices); % Number of slices to display
        for idx = 1:num_slices_to_show
            k = param.slice_indices(idx); % Get the current slice index
            subplot(ceil(sqrt(num_slices_to_show)), ceil(sqrt(num_slices_to_show)), idx);
            hold on;

            % Create a color map for the visualization
            colors = zeros(rows, cols, 3); % Initialize color array for RGB

            for i = 1:rows
                for j = 1:cols
                    % Assign colors based on the tissue type with the highest probability
                    if data_GmWmCsfSculpt_3D(i, j, k, 1) > max([data_GmWmCsfSculpt_3D(i, j, k, 2), data_GmWmCsfSculpt_3D(i, j, k, 3), data_GmWmCsfSculpt_3D(i, j, k, 4)])
                        colors(i, j, :) = [0.5, 0.5, 0.5]; % GM - Gray
                    elseif data_GmWmCsfSculpt_3D(i, j, k, 2) > max([data_GmWmCsfSculpt_3D(i, j, k, 1), data_GmWmCsfSculpt_3D(i, j, k, 3), data_GmWmCsfSculpt_3D(i, j, k, 4)])
                        colors(i, j, :) = [0, 1, 0]; % WM - Green
                    elseif data_GmWmCsfSculpt_3D(i, j, k, 3) > max([data_GmWmCsfSculpt_3D(i, j, k, 1), data_GmWmCsfSculpt_3D(i, j, k, 2), data_GmWmCsfSculpt_3D(i, j, k, 4)])
                        colors(i, j, :) = [0, 0, 1]; % CSF - Blue
                    else
                        colors(i, j, :) = [1, 0, 0]; % Bone - Red
                    end
                end
            end

            % Display the 2D slice with tissue colors
            imagesc(colors);
            axis image; % Keep aspect ratio
            title(['Slice ' num2str(k)]);
            xlabel('X-axis (Pixels)');
            ylabel('Y-axis (Pixels)');
            grid on;
            hold off;
        end
    end
end
