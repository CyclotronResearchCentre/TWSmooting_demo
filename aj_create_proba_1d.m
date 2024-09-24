function [data_GmWmCsf] = aj_create_proba_1d(data_1D)
    % AJ_CREATE_PROBA_1D Calculate tissue probability from intensity values
    %
    % This function takes a 1D array of intensity values (data_1D)
    % and calculates the probability of each intensity corresponding to
    % four types of tissue: Bone, CSF (Cerebrospinal Fluid), WM (White Matter),
    % and GM (Gray Matter). The probabilities are computed using Gaussian
    % distributions centered around reference intensity values for each 
    % tissue type. The choice of using a PDF allows us to model the 
    % likelihood of observing each intensity, considering its proximity 
    % to the expected intensity for each tissue type. This method provides
    % a more nuanced representation of tissue probability by assigning
    % higher probabilities to intensities close to the reference values.

    % Initialize the output matrix
    data_GmWmCsf = zeros(length(data_1D), 4); % 4 columns for Bone, CSF, WM, GM

    % Define the reference intensities for each tissue type
    ref_bone = 0.5;
    ref_csf = 0.2;
    ref_wm = 0.3;
    ref_gm = 0;
    
    % Define standard deviations for the normal distributions
    sigma_bone = 0.1; % Adjust as needed
    sigma_csf = 0.05; % Adjust as needed
    sigma_wm = 0.05;  % Adjust as needed
    sigma_gm = 0.05;  % Adjust as needed

    % Calculate probabilities using Gaussian distributions
    for i = 1:length(data_1D)
        intensity = data_1D(i);
        
        % Calculate Gaussian probabilities
        prob_bone = normpdf(intensity, ref_bone, sigma_bone);
        prob_csf = normpdf(intensity, ref_csf, sigma_csf);
        prob_wm = normpdf(intensity, ref_wm, sigma_wm);
        prob_gm = normpdf(intensity, ref_gm, sigma_gm);
        
        % Sum of probabilities
        total_prob = prob_bone + prob_csf + prob_wm + prob_gm;

        % Normalize to get probabilities
        data_GmWmCsf(i, 1) = prob_bone / total_prob; % Bone
        data_GmWmCsf(i, 2) = prob_csf / total_prob;  % CSF
        data_GmWmCsf(i, 3) = prob_wm / total_prob;   % WM
        data_GmWmCsf(i, 4) = prob_gm / total_prob;   % GM
    end
    
    % Visualization
    figure;
    hold on;
    
    % Assign colors based on tissue type probabilities
    colors = zeros(length(data_1D), 3); % Initialize color array
    for i = 1:length(data_1D)
        if data_GmWmCsf(i, 1) > max([data_GmWmCsf(i, 2), data_GmWmCsf(i, 3), data_GmWmCsf(i, 4)])
            colors(i, :) = [1, 0, 0]; % Bone - Red
        elseif data_GmWmCsf(i, 2) > max([data_GmWmCsf(i, 1), data_GmWmCsf(i, 3), data_GmWmCsf(i, 4)])
            colors(i, :) = [0, 0, 1]; % CSF - Blue
        elseif data_GmWmCsf(i, 3) > max([data_GmWmCsf(i, 1), data_GmWmCsf(i, 2), data_GmWmCsf(i, 4)])
            colors(i, :) = [0, 1, 0]; % WM - Green
        else
            colors(i, :) = [0.5, 0.5, 0.5]; % GM - Gray
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
