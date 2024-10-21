function [ph_anatomical, ph_noise, ph_var, ellipse_anatomical] = ...
    aj_phantom_variability(ph, ellipse, param, flag)
% aj_phantom_variability Adds anatomical variability and noise to a 3D aj_create_phantom_3d
%
% INPUTS:
% ph           : Original phantom volume
% ellipse      : Original ellipsoid parameters
% param        : see aj_phantom_default
% flag         : see aj_phantom_default
%
% OUTPUTS:
% ph_anatomical     : Phantom with anatomical variability
% ph_noise          : Phantom with noise
% ph_var            : Phantom with anatomical variability and noise
% ellipse_anatomical: Ellipsoid parameters with variability

%% Jitter factor
if isempty(param.jitter_factor) || param.jitter_factor < param.jitter_range(1) || param.jitter_factor > param.jitter_range(2)
    param.jitter_factor = mean(param.jitter_range); % default value = 0.05
    disp(['Jitter factor out of recommended range, set to default value of ', num2str(param.jitter_factor)]);
end

%% Noise level
if isempty(param.noise_level) || param.noise_level < param.noise_range(1) || param.noise_level > param.noise_range(2)
    param.noise_level = mean(param.noise_range);
    disp(['Noise level out of recommended range, set to default value of ', num2str(param.noise_level)]);
end

%% Original phantom with added anatomical variability
ellipse_anatomical = ellipse; % Make a copy of the original ellipse parameters

x_ellipse = 3; % keep the first two ellipses constant

nb_ellipse = size(ellipse);
subellipse = ellipse(x_ellipse:nb_ellipse(2), 2:10); % x_ellipse to 10th ellipses, 2nd to 10th parameters
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

ellipse_anatomical(x_ellipse:nb_ellipse(2), 2:10) = subellipse_anatomical;

[ph_anatomical, ~] = aj_phantom_create3D(ellipse_anatomical, param); 
% resulting ellipse is the same as the input ellipse_anatomical

%% Original phantom with added noise
ph_noise = ph + param.noise_level * randn(size(ph));

%% Smoothed phantom with anatomical variability & noise
% Apply noise before smoothing
if flag.noise_before_smoothing
    ph_noisy_anat = ph_anatomical + param.prenoise_level * randn(size(ph_anatomical));
else
    ph_noisy_anat = ph_anatomical;
end

% Apply Gaussian smoothing
wg = gausswin(param.sm_kern); % Create Gaussian kernel
wg = wg / sum(wg); % Normalize the kernel

% Apply 3D smoothing
smoothed_ph = convn(ph_noisy_anat, reshape(wg, [length(wg), 1, 1]), 'same');
smoothed_ph = convn(smoothed_ph, reshape(wg, [1, length(wg), 1]), 'same');
smoothed_ph = convn(smoothed_ph, reshape(wg, [1, 1, length(wg)]), 'same');

ph_var = smoothed_ph + param.noise_level * randn(size(smoothed_ph));

%% Optional: Visualize the results
if flag.plot_fig
    figure;
    
    slice_num = round(size(ph, 1) / 2);  % Get the middle slice

    min_val = min(ph(:));
    max_val = max(ph(:));

    % Conversion of the axes to millimeters
    voxel_size = param.voxreal_res;  % Resolution factor [mm/voxel]
    x_axis_mm = ((1:size(ph, 2)) - size(ph, 2)/2) * voxel_size;  % X-axis in mm
    y_axis_mm = ((1:size(ph, 1)) - size(ph, 1)/2) * voxel_size;  % Y-axis in mm

    % Plot subplots with shared color scale
    subplot(2, 2, 1);
    imagesc(x_axis_mm, y_axis_mm, squeeze(ph(:, :, slice_num)), [min_val, max_val]);
    title('Original Phantom (OP)');
    colormap gray;
    axis image;
    set(gca, 'YDir', 'normal'); % imagesc: By default, displays the matrix so that the first row is at the top.
    xlabel('X [mm]');
    ylabel('Y [mm]');

    subplot(2, 2, 2);
    imagesc(x_axis_mm, y_axis_mm, squeeze(ph_anatomical(:, :, slice_num)), [min_val, max_val]);
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
    imagesc(x_axis_mm, y_axis_mm, squeeze(ph_var(:, :, slice_num)), [min_val, max_val]);
    title('Smoothed P with anatomical variability & noise');
    colormap gray;
    axis image;
    set(gca, 'YDir', 'normal');
    xlabel('X [mm]');
    ylabel('Y [mm]');
end

end
