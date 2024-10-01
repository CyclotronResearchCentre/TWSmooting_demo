function [param, flag] = aj_proba_default()
    % AJ_PROBA_DEFAULT Returns default parameters for tissue probabilities
    %
    % This function defines and returns the default reference intensities
    % and standard deviations for each tissue type (Bone, CSF, WM, GM).
    
    % Parameters Definitions
    % Define the reference intensities for each tissue type
    param.ref_bone = 0.6;
    param.ref_csf = 0.2;
    param.ref_wm = 0.3;
    param.ref_gm = 0;

    % Define standard deviations for the normal distributions
    param.sigma_bone = 0.1;
    param.sigma_csf = 0.05;
    param.sigma_wm = 0.05;
    param.sigma_gm = 0.05;
    
    % Slices to plot in the 3D phantom (warning: flag has to be true)
    param.slice_indices = [15, 45, 64, 85, 110];
    
    % Flag: Used to manage additional behaviors
    flag.plot_fig = true;  % Flag to plot figures (set to true for plotting)
end
