function [data_1D, data_2D] = aj_phantom_extract1D2D(ph, param, plot_fig)
%--------------------------------------------------------------------------
% Function to extract 1D and 2D data from a 3D phantom
%
% INPUTS
% ph          : 3D phantom volume with noise and anatomical variability
% param       : Structure containing default parameters
% plot_fig    : Flag to plot figures for visual comparison
%
% OUTPUTS
% data_1D     : 1D profile of the phantom (intensity values along a line)
% data_2D     : 2D slice of the phantom (intensity values on a 2D plane)
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------

%% Dealing with inputs
if nargin < 2, plot_fig = false; end

%% Extracting 1D & 2D data
dim_y = size(ph, 2);
dim_z = size(ph, 3);

center_y = round(dim_y / 2);
center_z = round(dim_z / 2);

data_1D = squeeze(ph(:, center_y, center_z))';  % 1D profile along the x-axis
data_2D = ph(:, :, center_z);  % 2D slice at the center of the Z-axis volume

%% Optional: Visualize the results
if plot_fig
    figure;
    
    % Conversion factor from voxel to mm
    voxel_size = param.voxreal_res;  % mm/voxel
    
    % X-axis in mm for 1D profile
    x_axis_mm_1D = ((1:length(data_1D)) - length(data_1D)/2) * voxel_size;
    
    % X and Y axes in mm for 2D slice
    x_axis_mm_2D = ((1:size(data_2D, 2)) - size(data_2D, 2)/2) * voxel_size;
    y_axis_mm_2D = ((1:size(data_2D, 1)) - size(data_2D, 1)/2) * voxel_size;
    
    % Display the original 1D profile in mm
    subplot(1, 3, 1);
    plot(x_axis_mm_1D, data_1D);
    title('1D Profile');
    xlabel('X [mm]');
    ylabel('Intensity');
    
    % Display the original 2D slice in mm
    subplot(1, 3, 2);
    imagesc(x_axis_mm_2D, y_axis_mm_2D, data_2D);
    title('2D Slice');
    colormap gray;
    axis image;
    xlabel('X [mm]');
    ylabel('Y [mm]');
    
    % Display the 2D slice with the red line representing the position of the 1D profile
    subplot(1, 3, 3);
    imagesc(x_axis_mm_2D, y_axis_mm_2D, data_2D);
    hold on;
    % The vertical red line is also scaled to match mm
    plot([x_axis_mm_2D(center_y), x_axis_mm_2D(center_y)], [y_axis_mm_2D(1), y_axis_mm_2D(end)], 'r', 'LineWidth', 2);
    hold off;
    title('2D Slice with 1D Profile Line');
    colormap gray;
    axis image;
    xlabel('X [mm]');
    ylabel('Y [mm]');
end


end