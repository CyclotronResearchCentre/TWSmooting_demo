function [ph_anat, ph_fn] = ...
    aj_phantom_variability(ph_orig, model_orig, param, flag)
%--------------------------------------------------------------------------
% aj_phantom_variability Adds anatomical variability and noise to a 3D aj_create_phantom_3d
%
% INPUTS
% ph_orig:      Original phantom volume (3D matrix)
% model_orig:   Original ellipsoid parameters from the original model
% param:        structure of parameters from aj_phantom_default
% flag:         structure of flags from aj_phantom_default
%
% OUTPUTS
% ph_anat:      Phantom with anatomical variability
% ph_noise:     Phantom with noise
% ph_fn:        Phantom with anatomical variability and noise
% model_anat:   Ellipsoid parameters model with anatomical variability
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------

%% Dealing with inputs
% Jitter factor
if isempty(param.jitter_factor) || param.jitter_factor < param.jitter_range(1) || param.jitter_factor > param.jitter_range(2)
    param.jitter_factor = mean(param.jitter_range); % default value = 0.05
    disp(['Jitter factor out of recommended range, set to default value of ', num2str(param.jitter_factor)]);
end

% Noise level
if isempty(param.noise_level) || param.noise_level < param.noise_range(1) || param.noise_level > param.noise_range(2)
    param.noise_level = mean(param.noise_range);
    disp(['Noise level out of recommended range, set to default value of ', num2str(param.noise_level)]);
end

%% Do the job
% Applying a anatomical variation to the original phantom
model_anat = model_orig;
x_ellipse = 3; % keep the first two ellipses constant
nb_ellipse = size(model_orig);
subellipse = model_orig(x_ellipse:nb_ellipse(2), 2:10); % x_ellipse to 10th ellipses, 2nd to 10th parameters
subellipse_anatomical = subellipse; % Create a copy to apply changes

for l = 1:size(subellipse, 1) % Iterate over ellipses (rows)
    % Calculate ellipse size (e.g., the mean of the dimensions)
    ellipse_size = mean(subellipse(l, 2:4)); % Assuming columns 2-4 contain radii

    % Adjust jitter based on the size of the ellipse
    adjusted_jitter = param.jitter_factor * (ellipse_size / max(subellipse(:, 2:4), [], 'all'));

    for k = 1:size(subellipse, 2) % Iterate over parameters (columns)
        if subellipse(l, k) == 0
            subellipse_anatomical(l, k) = adjusted_jitter * randn(1);  % Gaussian noise centered on 0
        else
            subellipse_anatomical(l, k) = subellipse(l, k) + ...
                (2 * randi([0, 1]) - 1) * adjusted_jitter * subellipse(l, k) * randn(1);
        end
    end
end

model_anat(x_ellipse:nb_ellipse(2), 2:10) = subellipse_anatomical;
ph_anat = aj_phantom_create3D(model_anat, param);

% Adding prenoise before smoothing (option)
if flag.noise_before_smoothing
    ph_noisy_anat = ph_anat + param.prenoise_level * randn(size(ph_anat)); % noise with a normal distribution
else
    ph_noisy_anat = ph_anat;
end

% Applying 3D Gaussian smoothing 
smooth_anat_ph = imgaussfilt3(ph_noisy_anat, 'FilterSize', param.fwhm);

% Adding noise with a normal distribution to the original phantom
ph_noise = ph_orig + param.noise_level * randn(size(ph_orig));

% Smoothed anatomical variation + noise
ph_fn = smooth_anat_ph + param.noise_level * randn(size(smooth_anat_ph));

%% Result visualization
if flag.plot_fig
    figure;
    
    slice_num = round(size(ph_orig, 1) / 2);  % Get the middle slice

    min_val = min(ph_orig(:));
    max_val = max(ph_orig(:));

    % Conversion of the axes to millimeters
    voxel_size = param.voxreal_res;  % Resolution factor [mm/voxel]
    x_axis_mm = ((1:size(ph_orig, 2)) - size(ph_orig, 2)/2) * voxel_size;  % X-axis in mm
    y_axis_mm = ((1:size(ph_orig, 1)) - size(ph_orig, 1)/2) * voxel_size;  % Y-axis in mm

    % Plot subplots with shared color scale
    subplot(2, 2, 1);
    imagesc(x_axis_mm, y_axis_mm, squeeze(ph_orig(:, :, slice_num)), [min_val, max_val]);
    title('Original Phantom (OP)');
    colormap gray;
    axis image;
    set(gca, 'YDir', 'normal'); % imagesc: By default, displays the matrix so that the first row is at the top.
    xlabel('X [mm]');
    ylabel('Y [mm]');

    subplot(2, 2, 2);
    imagesc(x_axis_mm, y_axis_mm, squeeze(ph_anat(:, :, slice_num)), [min_val, max_val]);
    title('OP with anatomical variability');
    colormap gray;
    axis image;
    set(gca, 'YDir', 'normal');
    xlabel('X [mm]');
    ylabel('Y [mm]');

    subplot(2, 2, 3);
    imagesc(x_axis_mm, y_axis_mm, squeeze(ph_noise(:, :, slice_num)), [min_val, max_val]);
    title('OP with noise');
    colormap gray;
    axis image;
    set(gca, 'YDir', 'normal');
    xlabel('X [mm]');
    ylabel('Y [mm]');

    subplot(2, 2, 4);
    imagesc(x_axis_mm, y_axis_mm, squeeze(ph_fn(:, :, slice_num)), [min_val, max_val]);
    title('Smoothed P with anatomical variability & noise');
    colormap gray;
    axis image;
    set(gca, 'YDir', 'normal');
    xlabel('X [mm]');
    ylabel('Y [mm]');
end

end
