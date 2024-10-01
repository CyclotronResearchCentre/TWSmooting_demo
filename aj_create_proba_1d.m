function [data_GmWmCsfSculpt_1D] = aj_create_proba_1d(data_1D, param, plot_fig)
    % AJ_CREATE_PROBA_1D Calculate tissue probability from intensity values
    %
    % This function takes a 1D array of intensity values (data_1D)
    % and calculates the probability of each intensity corresponding to
    % four types of tissue: Bone Sculpt, CSF (Cerebrospinal Fluid), WM (White Matter),
    % and GM (Gray Matter). The probabilities are computed using Gaussian
    % distributions centered around reference intensity values for each 
    % tissue type. The choice of using a PDF allows us to model the 
    % likelihood of observing each intensity, considering its proximity 
    % to the expected intensity for each tissue type. This method provides
    % a more nuanced representation of tissue probability by assigning
    % higher probabilities to intensities close to the reference values.

    % Initialize the output matrix
    data_GmWmCsfSculpt_1D = zeros(length(data_1D), 4); % 4 columns for GM, WM, CSF, Bone Sculpt

    % Calculate probabilities using Gaussian distributions
    for i = 1:length(data_1D)
        intensity = data_1D(i);
        
        % Calculate Gaussian probabilities
        prob_gm = normpdf(intensity, param.ref_gm, param.sigma_gm);
        prob_wm = normpdf(intensity, param.ref_wm, param.sigma_wm);
        prob_csf = normpdf(intensity, param.ref_csf, param.sigma_csf);
        prob_bone = normpdf(intensity, param.ref_bone, param.sigma_bone);
        
        % Sum of probabilities
        total_prob = prob_gm + prob_wm + prob_csf + prob_bone;

        % Normalize to get probabilities
        data_GmWmCsfSculpt_1D(i, 1) = prob_gm / total_prob;   % GM
        data_GmWmCsfSculpt_1D(i, 2) = prob_wm / total_prob;   % WM
        data_GmWmCsfSculpt_1D(i, 3) = prob_csf / total_prob;  % CSF
        data_GmWmCsfSculpt_1D(i, 4) = prob_bone / total_prob; % Bone Sculpt
    end
    
    % Visualization
    if plot_fig
        figure;
        hold on;

        % Assign colors based on tissue type probabilities
        colors = zeros(length(data_1D),3); % Initialize color array
        for i = 1:length(data_1D)
            if data_GmWmCsfSculpt_1D(i, 1) > max([data_GmWmCsfSculpt_1D(i, 2), data_GmWmCsfSculpt_1D(i, 3), data_GmWmCsfSculpt_1D(i, 4)])
                colors(i, :) = [0.5, 0.5, 0.5]; % GM - Gray
            elseif data_GmWmCsfSculpt_1D(i, 2) > max([data_GmWmCsfSculpt_1D(i, 1), data_GmWmCsfSculpt_1D(i, 3), data_GmWmCsfSculpt_1D(i, 4)])
                colors(i, :) = [0, 1, 0]; % WM - Green
            elseif data_GmWmCsfSculpt_1D(i, 3) > max([data_GmWmCsfSculpt_1D(i, 1), data_GmWmCsfSculpt_1D(i, 2), data_GmWmCsfSculpt_1D(i, 4)])
                colors(i, :) = [0, 0, 1]; % CSF - Blue
            else
                colors(i, :) = [1, 0, 0]; % Bone - Red
            end
        end

        % Create scatter plot based on index
        scatter(1:length(data_1D), data_1D, 36, colors, 'filled');

        % Formatting the plot
        xlabel('Index');
        ylabel('Intensity Values');
        title('Tissue Probability by Intensity');
        grid on;
        hold off;
    end
end
