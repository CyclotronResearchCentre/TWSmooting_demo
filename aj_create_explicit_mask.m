function [iexMask] = aj_create_explicit_mask(tissue_proba)
    % Fonction pour créer un masque explicite pour la substance grise (GM),
    % la substance blanche (WM), et le CSF à partir des probabilités tissulaires lissées.
    % Compatible avec les données 1D, 2D, et 3D.
    %
    % INPUT :
    % tissue_proba : Cell array {3 x 1} contenant les probabilités pour GM, WM, CSF,
    %                ou une matrice [3 x Npt] en 1D, [3 x Nx x Ny] en 2D, [3 x Nx x Ny x Nz] en 3D.
    %
    % OUTPUT :
    % iexMask : Masque explicite pour GM, WM et CSF. Dimension [3 x ...] avec la même taille que tissue_proba.
    %           Le masque contient 3 couches : GM, WM, et CSF.
    %
    % Conditions :
    % - GM : Probabilité de GM > WM et CSF, et probabilité > 0.2
    % - WM : Probabilité de WM > GM et CSF, et probabilité > 0.2
    % - CSF : Probabilité de CSF > GM et WM, et probabilité > 0.2

    % Vérification de l'entrée : Si ce n'est pas un cell array, on le convertit.
    if ~iscell(tissue_proba)
        tmp = tissue_proba;
        tissue_proba = cell(3,1);
        for ii = 1:3
            tissue_proba{ii} = tmp(ii, :);  % Distribution des probabilités entre GM, WM, CSF
        end
        need2avg = 0;  % Pas besoin de moyenne si un seul sujet
    else
        need2avg = 1;  % Moyenne requise si plusieurs sujets
    end
    
    % Identifier la dimension des données
    data_size = size(tissue_proba{1});
    if isvector(tissue_proba{1})
        disp('Ceci est un vecteur.');
        data_dims = 1;
    elseif ismatrix(tissue_proba{1})
        disp('Ceci est une matrice.');
        data_dims = 2;
    else
        data_dims = length(data_size);  % Nombre de dimensions des données (1, 2 ou 3)
        disp('eeee');
        if data_dims == 3
            disp('Ceci est une matrice 3D.');
        end
    end

    % Moyenne des probabilités tissulaires si plusieurs sujets
    if need2avg
        avgP_GmWmCsf = zeros([3, data_size]);  % Initialisation
        for ii = 1:3
            avgP_GmWmCsf(ii, :) = mean(tissue_proba{ii}, 'all');  % Moyenne des probabilités sur tous les sujets
        end
    else
        avgP_GmWmCsf = zeros([3, data_size]);
        for ii = 1:3
            avgP_GmWmCsf(ii, :) = tissue_proba{ii};  % Utilisation des probabilités uniques si un seul sujet
        end
    end

    % Initialisation du masque explicite avec les mêmes dimensions
    iexMask = false([3, data_size]);

    % Calcul du masque explicite pour GM, WM et CSF en fonction des dimensions
    if data_dims == 2  % 2D
        % GM
        iexMask(1, :, :) = avgP_GmWmCsf(1, :, :) > avgP_GmWmCsf(2, :, :) & ...
                           avgP_GmWmCsf(1, :, :) > avgP_GmWmCsf(3, :, :) & ...
                           avgP_GmWmCsf(1, :, :) > 0.2;
        % WM
        iexMask(2, :, :) = avgP_GmWmCsf(2, :, :) > avgP_GmWmCsf(1, :, :) & ...
                           avgP_GmWmCsf(2, :, :) > avgP_GmWmCsf(3, :, :) & ...
                           avgP_GmWmCsf(2, :, :) > 0.2;
        % CSF
        iexMask(3, :, :) = avgP_GmWmCsf(3, :, :) > avgP_GmWmCsf(1, :, :) & ...
                           avgP_GmWmCsf(3, :, :) > avgP_GmWmCsf(2, :, :) & ...
                           avgP_GmWmCsf(3, :, :) > 0.2;

    elseif data_dims == 3  % 3D
        % GM
        iexMask(1, :, :, :) = avgP_GmWmCsf(1, :, :, :) > avgP_GmWmCsf(2, :, :, :) & ...
                              avgP_GmWmCsf(1, :, :, :) > avgP_GmWmCsf(3, :, :, :) & ...
                              avgP_GmWmCsf(1, :, :, :) > 0.2;
        % WM
        iexMask(2, :, :, :) = avgP_GmWmCsf(2, :, :, :) > avgP_GmWmCsf(1, :, :, :) & ...
                              avgP_GmWmCsf(2, :, :, :) > avgP_GmWmCsf(3, :, :, :) & ...
                              avgP_GmWmCsf(2, :, :, :) > 0.2;
        % CSF
        iexMask(3, :, :, :) = avgP_GmWmCsf(3, :, :, :) > avgP_GmWmCsf(1, :, :, :) & ...
                              avgP_GmWmCsf(3, :, :, :) > avgP_GmWmCsf(2, :, :, :) & ...
                              avgP_GmWmCsf(3, :, :, :) > 0.2;

    else  % 1D (par défaut)
        % GM
        iexMask(1, :) = avgP_GmWmCsf(1, :) > avgP_GmWmCsf(2, :) & ...
                        avgP_GmWmCsf(1, :) > avgP_GmWmCsf(3, :) & ...
                        avgP_GmWmCsf(1, :) > 0.2;
        % WM
        iexMask(2, :) = avgP_GmWmCsf(2, :) > avgP_GmWmCsf(1, :) & ...
                        avgP_GmWmCsf(2, :) > avgP_GmWmCsf(3, :) & ...
                        avgP_GmWmCsf(2, :) > 0.2;
        % CSF
        iexMask(3, :) = avgP_GmWmCsf(3, :) > avgP_GmWmCsf(1, :) & ...
                        avgP_GmWmCsf(3, :) > avgP_GmWmCsf(2, :) & ...
                        avgP_GmWmCsf(3, :) > 0.2;
    end
end
