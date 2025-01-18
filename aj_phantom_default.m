function [model, param, flag] = aj_phantom_default()
%--------------------------------------------------------------------------
% Function to get default parameters, flags and model
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------
% Model Definitions: Available phantom models
model.shepp_logan = @shepp_logan;
model.modified_shepp_logan = @modified_shepp_logan;
model.yu_ye_wang = @yu_ye_wang;

% Define default model type
param.model_type = 'modified_shepp_logan';  % Set default model

% Parameters Definitions: Default parameters for generating phantoms
param.n = 1;                    % Default number of phantoms to generate
param.grid_size = 128;          % Default grid size for 3D phantom
param.voxreal_res = 0.68;          % Resolution factor between the voxel size and the real size [mm/voxel]
param.modelreal_res = 1;        % Resolution factor between the model

param.jitter_range = [0, 0.1];
param.jitter_factor = 0.05;     % Default jitter factor for anatomical variability
param.prenoise_level = 0.01;    % Noise to apply before smoothing if flag.noise_before_smoothing is true
param.fwhm = 5;                 % gaussian kernel size for smoothing (to get more realistic images)
param.noise_range = [0, 0.05];
param.noise_level = 0.02;       % Default noise level to add to phantoms

param.proba_noise = 0.1;        % Noise level for probability maps

% Flag: Used to manage additional behaviors
flag.noise_before_smoothing = false;    % Flag to apply noise before smoothing (and still apply noise after)
flag.plot_fig = false;                   % Flag to plot figures (set to true for plotting)
flag.ph_stat = false;                    % Flag to compute RMSE between phantoms (anat) if n>1
flag.plot_ph_stat = false;              % Flag to plot histogram and Q-Q plot of RMSE values (only if flag.ph_stat == true)

flag.del_prevResults = true;
flag.data_1D2D = false;
end

% A: Intensity (positive for adding intensity, negative for subtracting).
% a, b, c: Semi-axis lengths in the x, y and z directions, respectively.
% x0, y0, z0: Center coordinates of the ellipsoid.
% phi, theta, psi: Rotation angles in degrees, corresponding to rotations around the x, y and z axes.

% Default ellipsoids for Shepp-Logan phantom
function e = shepp_logan()
    e = modified_shepp_logan();  % Use the modified Shepp-Logan as base
    % The intensity values are modified here for the Shepp-Logan phantom
    e(:, 1) = [1, -0.98, -0.02, -0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01];  
end

% Default ellipsoids for Modified Shepp-Logan phantom
function e = modified_shepp_logan()
    % Each row represents one ellipsoid with the following parameters:
    % [A, a, b, c, x0, y0, z0, phi, theta, psi]
    e = [ ...
         1,    0.69,   0.92,   0.81,  0,     0,      0,     0,    0,   0;   % Ellipsoid 1
        -0.8,  0.6624, 0.874,  0.78,  0,    -0.0184, 0,     0,    0,   0;   % Ellipsoid 2
        -0.2,  0.11,   0.31,   0.22,  0.22,  0,      0,    -18,   0,  10;   % Ellipsoid 3
        -0.2,  0.16,   0.41,   0.28, -0.22,  0,      0,     18,   0,  10;   % Ellipsoid 4
         0.1,  0.21,   0.25,   0.41,  0,     0.35,  -0.15,  0,    0,   0;   % Ellipsoid 5
         0.1,  0.046,  0.046,  0.05,  0,     0.1,    0.25,  0,    0,   0;   % Ellipsoid 6
         0.1,  0.046,  0.046,  0.05,  0,    -0.1,    0.25,  0,    0,   0;   % Ellipsoid 7
         0.1,  0.046,  0.023,  0.05, -0.08, -0.605,  0,     0,    0,   0;   % Ellipsoid 8
         0.1,  0.023,  0.023,  0.02,  0,    -0.606,  0,     0,    0,   0;   % Ellipsoid 9
         0.1,  0.023,  0.046,  0.02,  0.06, -0.605,  0,     0,    0,   0;   % Ellipsoid 10
    ];
end

% Default ellipsoids for Yu-Ye-Wang phantom
function e = yu_ye_wang()
    % Each row represents one ellipsoid with the following parameters:
    % [A, a, b, c, x0, y0, z0, phi, theta, psi]
    e = [ ...
         1,    0.69,   0.92,   0.9,    0,     0,      0,     0,    0,   0;   % Ellipsoid 1
        -0.8,  0.6624, 0.874,  0.88,   0,     0,      0,     0,    0,   0;   % Ellipsoid 2
        -0.2,  0.41,   0.16,   0.21,  -0.22,  0,     -0.25, 108,   0,   0;   % Ellipsoid 3
        -0.2,  0.31,   0.11,   0.22,   0.22,  0,     -0.25,  72,   0,   0;   % Ellipsoid 4
         0.2,  0.25,   0.21,   0.31,   0,     0.35,  -0.25,  0,    0,   0;   % Ellipsoid 5
         0.2,  0.046,  0.046,  0.05,   0,     0.1,    0.625, 0,    0,   0;   % Ellipsoid 6
         0.2,  0.046,  0.046,  0.05,   0,    -0.1,    0.625, 0,    0,   0;   % Ellipsoid 7
         0.2,  0.023,  0.046,  0.05,  -0.08, -0.605,  0,     0,    0,   0;   % Ellipsoid 8
         0.2,  0.023,  0.023,  0.02,   0,    -0.606,  0,     0,    0,   0;   % Ellipsoid 9
         0.2,  0.023,  0.046,  0.02,   0.06, -0.605,  0,     0,    0,   0;   % Ellipsoid 10
    ];
end

