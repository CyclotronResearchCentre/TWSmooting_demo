function [data_GmWmCsfSculpt_2D] = aj_create_proba_2d(data_2D, param, plot_fig)
    % AJ_CREATE_PROBA_2D Calculate tissue probability for 2D image
    %
    % This function takes a 2D image of intensity values (data_2D)
    % and calculates the probability of each pixel intensity corresponding to
    % four types of tissue: Bone Sculpt, CSF, WM, and GM. The probabilities 
    % are computed using Gaussian distributions.
    
    [rows, cols] = size(data_2D);
    
    % Initialize the output matrix for storing probabilities (Bone, CSF, WM, GM)
    data_GmWmCsfSculpt_2D = zeros(rows, cols, 4);

    % Calculate probabilities for each pixel in the 2D image
    for i = 1:rows
        for j = 1:cols
            intensity = data_2D(i, j);

            % Calculate Gaussian probabilities
            prob_gm = normpdf(intensity, param.ref_gm, param.sigma_gm);
            prob_wm = normpdf(intensity, param.ref_wm, param.sigma_wm);
            prob_csf = normpdf(intensity, param.ref_csf, param.sigma_csf);
            prob_bone = normpdf(intensity, param.ref_bone, param.sigma_bone);

            % Sum of probabilities
            total_prob = prob_gm + prob_wm + prob_csf + prob_bone;

            % Normalize to get probabilities
            data_GmWmCsfSculpt_2D(i, j, 1) = prob_gm / total_prob;   % GM
            data_GmWmCsfSculpt_2D(i, j, 2) = prob_wm / total_prob;   % WM
            data_GmWmCsfSculpt_2D(i, j, 3) = prob_csf / total_prob;  % CSF
            data_GmWmCsfSculpt_2D(i, j, 4) = prob_bone / total_prob; % Bone Sculpt
        end
    end

    % Visualization of the 2D image
    if plot_fig
        figure;
        hold on;

        % Create a color map for the visualization
        colors = zeros(rows, cols, 3); % Initialize color array for RGB
        for i = 1:rows
            for j = 1:cols
                % Assign colors based on the tissue type with the highest probability
                if data_GmWmCsfSculpt_2D(i, j, 1) > max([data_GmWmCsfSculpt_2D(i, j, 2), data_GmWmCsfSculpt_2D(i, j, 3), data_GmWmCsfSculpt_2D(i, j, 4)])
                    colors(i, j, :) = [0.5, 0.5, 0.5]; % GM - Gray
                elseif data_GmWmCsfSculpt_2D(i, j, 2) > max([data_GmWmCsfSculpt_2D(i, j, 1), data_GmWmCsfSculpt_2D(i, j, 3), data_GmWmCsfSculpt_2D(i, j, 4)])
                    colors(i, j, :) = [0, 1, 0]; % WM - Green
                elseif data_GmWmCsfSculpt_2D(i, j, 3) > max([data_GmWmCsfSculpt_2D(i, j, 1), data_GmWmCsfSculpt_2D(i, j, 2), data_GmWmCsfSculpt_2D(i, j, 4)])
                    colors(i, j, :) = [0, 0, 1]; % CSF - Blue
                else
                    colors(i, j, :) = [1, 0, 0]; % Bone - Red
                end
            end
        end

        % Display the 2D image with tissue colors
        imagesc(colors);
        axis image; % Keep aspect ratio
        title('Tissue Probability Map in 2D Image');
        xlabel('X-axis (Pixels)');
        ylabel('Y-axis (Pixels)');
        grid on;
        hold off;
    end
end
