function [data_1D, data_2D] = aj_extract_data_from_phantom(ph_var, plot_fig)
% aj_extract_data_from_phantom Extract 1D and 2D data from a 3D phantom
% 
% INPUTS:
% noisy_p      : Phantom 3D volume with noise and anatomical variability
% plot_fig     : Flag to plot figures for visual comparison
%
% OUTPUTS:
% data_1D          : 1D profile of the phantom (intensity values along a line)
% data_2D          : 2D slice of the phantom (intensity values on a 2D plane)

%% Default parameters
if nargin < 2, plot_fig = false; end

%% Extraction des données 1D & 2D
dim_y = size(ph_var, 2);
dim_z = size(ph_var, 3);

center_y = round(dim_y / 2);
center_z = round(dim_z / 2);

data_1D = squeeze(ph_var(:, center_y, center_z));  % Profil 1D sur l'axe x
data_2D = ph_var(:, :, center_z);  % Coupe 2D au centre du volume d'axe Z

%% Affichage optionnel
if plot_fig
    figure;
    
    % Afficher le profil 1D original
    subplot(1, 3, 1);
    plot(data_1D);
    title('1D Profile');
    xlabel('X-axis index');
    ylabel('Intensity');
    
    % Afficher la coupe 2D originale
    subplot(1, 3, 2);
    imagesc(data_2D);
    title('2D Slice');
    colormap gray;
    axis image;
    
    % Afficher la coupe 2D avec la ligne rouge représentant la position du profil 1D
    subplot(1, 3, 3);
    imagesc(data_2D);
    hold on;
    plot([center_y, center_y], [1, size(data_2D, 1)], 'r', 'LineWidth', 2);  % Ligne rouge verticale
    hold off;
    title('2D Slice with 1D Profile Line');
    colormap gray;
    axis image;
end

end