function [P_signal, P_GmWmCsf, T_names] = aj_create_data_3D(r_jitter, T_noise, plot_fig, n)
    % Function to create some synthetic 3D data:
    % - synthetic signal profile for a 3D Modified Shepp-Logan phantom
    % - some tissue probabilities for GM, WM & CSF
    
    % INPUT:
    % r_jitter : random jitter to add variability to ellipsoids, default 0
    % T_noise  : noise variance for the signal, default [5 5 10]
    % plot_fig : flag to plot 3D slices or not, default 0
    % n        : grid size for the 3D phantom, default 128
    
    if nargin < 4, n = 128; end
    if nargin < 3, plot_fig = 0; end
    if nargin < 2, T_noise = [5 5 10]; end
    if nargin < 1, r_jitter = 0; end
    
    % Define tissue names and default intensity values
    T_names = {'GM', 'WM', 'CSF'};
    T_signal = [50, 100, 5];
    
    % Create a 3D Modified Shepp-Logan phantom
    [P_signal, ellipsoids] = aj_phantom3d('modified shepp-logan', n);
    
    % Add jitter (variability) to ellipsoids if requested
    if r_jitter
        for k = 1:size(ellipsoids, 1)
            ellipsoids(k, 5:7) = ellipsoids(k, 5:7) + rand(1, 3) * r_jitter;
        end
        P_signal = aj_phantom3d(ellipsoids, n); % Recompute with jitter
    end
    
    % Add noise to the signal based on T_noise
    P_signal = P_signal + randn(size(P_signal)) .* reshape(T_noise, 1, 1, 1, length(T_noise));
    
    % Placeholder for tissue probability (simplified in this example)
    P_GmWmCsf = zeros(3, n, n, n); % For GM, WM, CSF
    
%     % Simulate tissue probabilities (simplified, could be adjusted based on ellipsoids)
%     P_GmWmCsf(1, :, :, :) = P_signal > 0.5;  % GM
%     P_GmWmCsf(2, :, :, :) = P_signal <= 0.5 & P_signal > 0.2;  % WM
%     P_GmWmCsf(3, :, :, :) = P_signal <= 0.2;  % CSF

    % Boucle sur les segments (ou ellipses)
    for ii = 1:Nsegm
        % Ici vous devez définir les indices dans les trois dimensions pour
        % les segments de tissus et affecter les valeurs de probabilité.
        % Par exemple, pour le premier tissu (GM), vous pourriez avoir :

        % L'indice GM (par exemple) sera rempli de manière spécifique dans les
        % dimensions x, y, z.

        P_GmWmCsf(1, :, :, :) = P_signal > 0.5;  % Exemple simple pour GM
        P_GmWmCsf(2, :, :, :) = P_signal > 0.2 & P_signal <= 0.5;  % Exemple pour WM
        P_GmWmCsf(3, :, :, :) = P_signal <= 0.2;  % Exemple pour CSF
    end
    
    % Plot figures if requested
    if plot_fig
        figure;
        slice_view = squeeze(P_signal(n/2, :, :));
        imagesc(slice_view);
        title('Middle Slice of the 3D Phantom');
        colormap gray;
        colorbar;
    end
end
