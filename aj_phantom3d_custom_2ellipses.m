function [p, ellipse] = aj_phantom3d_custom_2ellipses(params_ellipse1, params_ellipse2, n)
    % PHANTOM3D_CUSTOM_2ELLIPSES Génère un phantom 3D avec 2 ellipses personnalisées
    %   [P, ELLIPSE] = AJ_PHANTOM3D_CUSTOM_2ELLIPSES(PARAMS_ELLIPSE1, PARAMS_ELLIPSE2, N)
    %   génère un phantom avec deux ellipses, dont les paramètres sont spécifiés.
    %
    %   PARAMS_ELLIPSE1 et PARAMS_ELLIPSE2 sont des vecteurs contenant les
    %   10 paramètres de chaque ellipse :
    %     [A, a, b, c, x0, y0, z0, phi, theta, psi].
    %
    %   N spécifie la taille de la grille 3D.
    %
    % Exemple d'utilisation :
    %   params1 = [1, 0.5, 0.3, 0.2, 0, 0, 0, 0, 0, 0];  % Paramètres ellipse 1
    %   params2 = [0.8, 0.4, 0.2, 0.3, 0.3, 0, 0, 30, 0, 45];  % Paramètres ellipse 2
    %   [p, e] = aj_phantom3d_custom_2ellipses(params1, params2, 128);

    % Combine the parameters for both ellipses into a matrix
    ellipse = [params_ellipse1];% params_ellipse2];

    % Call the original phantom generation function
    p = aj_phantom3d(ellipse, n);

    % Display the central slice of the phantom for visualization
    figure;
    imshow(squeeze(p(round(n/2), :, :)), []);
    title('Central Slice of 3D Phantom with 2 Ellipses');
end
