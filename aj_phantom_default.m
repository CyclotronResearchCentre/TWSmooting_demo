function [model, param, flag] = aj_phantom_default()

    % Model Definitions: Available phantom models
    model.shepp_logan = @shepp_logan;
    model.modified_shepp_logan = @modified_shepp_logan;
    model.yu_ye_wang = @yu_ye_wang;

    % Parameters Definitions: Default parameters for generating phantoms
    param.n = 1;  % Default number of phantoms to generate
    param.grid_size = 128;  % Default grid size for 3D phantom
    param.FOV = 256;  % Field of View in mm (this is the physical size of the phantom in mm)
    param.jitter_range = [0, 0.1];
    param.jitter_factor = 0.05;  % Default jitter factor for anatomical variability
    param.prenoise_level = 0.01; % Noise to apply before smoothing if flag.noise_before_smoothing is true
    param.sm_kern = 4; % taille du noyau gaussien pour le lissage pour rendre les données + réalistes
    param.noise_range = [0, 0.05];
    param.noise_level = 0.01;  % Default noise level to add to phantoms

    % Flag: Used to manage additional behaviors
    flag.noise_before_smoothing = false; % Flag to apply noise before smoothing (and still apply noise after)
    flag.plot_fig = true;  % Flag to plot figures (set to true for plotting)
    flag.verbose = false;  % Set to true for verbose output (optional)

end

% Default ellipsoids for Shepp-Logan phantom
function e = shepp_logan()
    e = modified_shepp_logan();  % Use the modified Shepp-Logan as base
    e(:, 1) = [1, -0.98, -0.02, -0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01];  % Intensity values
end

% Default ellipsoids for Modified Shepp-Logan phantom
function e = modified_shepp_logan()
    e = [1, 0.69, 0.92, 0.81, 0, 0, 0, 0, 0, 0;
         -0.8, 0.6624, 0.874, 0.78, 0, -0.0184, 0, 0, 0, 0;
         -0.2, 0.11, 0.31, 0.22, 0.22, 0, 0, -18, 0, 10;
         -0.2, 0.16, 0.41, 0.28, -0.22, 0, 0, 18, 0, 10;
         0.1, 0.21, 0.25, 0.41, 0, 0.35, -0.15, 0, 0, 0;
         0.1, 0.046, 0.046, 0.05, 0, 0.1, 0.25, 0, 0, 0;
         0.1, 0.046, 0.046, 0.05, 0, -0.1, 0.25, 0, 0, 0;
         0.1, 0.046, 0.023, 0.05, -0.08, -0.605, 0, 0, 0, 0;
         0.1, 0.023, 0.023, 0.02, 0, -0.606, 0, 0, 0, 0;
         0.1, 0.023, 0.046, 0.02, 0.06, -0.605, 0, 0, 0, 0];
end

% Default ellipsoids for Yu-Ye-Wang phantom
function e = yu_ye_wang()
    e = [1, 0.69, 0.92, 0.9, 0, 0, 0, 0, 0, 0;
         -0.8, 0.6624, 0.874, 0.88, 0, 0, 0, 0, 0, 0;
         -0.2, 0.41, 0.16, 0.21, -0.22, 0, -0.25, 108, 0, 0;
         -0.2, 0.31, 0.11, 0.22, 0.22, 0, -0.25, 72, 0, 0;
         0.2, 0.25, 0.21, 0.31, 0, 0.35, -0.25, 0, 0, 0;
         0.2, 0.046, 0.046, 0.05, 0, 0.1, 0.625, 0, 0, 0;
         0.2, 0.046, 0.046, 0.05, 0, -0.1, 0.625, 0, 0, 0;
         0.2, 0.023, 0.046, 0.05, -0.08, -0.605, 0, 0, 0, 0;
         0.2, 0.023, 0.023, 0.02, 0, -0.606, 0, 0, 0, 0;
         0.2, 0.023, 0.046, 0.02, 0.06, -0.605, 0, 0, 0, 0];
end
