function [iexMask] = aj_create_explicit_mask(tissue_proba)
    % Fonction pour créer un masque explicite pour la substance grise (GM),
    % la substance blanche (WM), et le CSF à partir des probabilités tissulaires lissées.
    %
    % INPUT :
    % tissue_proba : Matrice [4 x Npt] (ou [4 x Nx x Ny] pour 2D, ou [4 x Nx x Ny x Nz] pour 3D)
    %                contenant les probabilités pour GM, WM, CSF et BoneSculpt.
    %
    % OUTPUT :
    % iexMask : Masque explicite pour GM, WM et CSF [3 x Npt] (ou [3 x Nx x Ny] ou [3 x Nx x Ny x Nz]).
    %
    % Conditions :
    % - GM : Probabilité de GM > WM et CSF, et probabilité > 0.2
    % - WM : Probabilité de WM > GM et CSF, et probabilité > 0.2
    % - CSF : Probabilité de CSF > GM et WM, et probabilité > 0.2

    % Vérification de l'entrée : Si ce n'est pas un cell array, le convertir.
    if ~iscell(tissue_proba)
        tissue_proba = num2cell(tissue_proba, 1);
    end

    % Nombre de points (Npt) dans les données
    Npt = size(tissue_proba{1}, 2);
    
    % Nombre de dimensions
    data_dims = ndims(tissue_proba{1});

    % Moyenne des probabilités tissulaires
    avgP_GmWmCsf = zeros(3, Npt);
    for ii = 1:3
        avgP_GmWmCsf(ii, :) = mean(tissue_proba{ii}, 2);  % Moyenne sur la dimension des sujets
    end

    % Initialiser iexMask en fonction des dimensions
    iexMask = false(3, Npt);

    % Calcul du masque explicite pour GM, WM et CSF
    iexMask(1, :) = avgP_GmWmCsf(1, :) > avgP_GmWmCsf(2, :) & ...  % GM > WM
                    avgP_GmWmCsf(1, :) > avgP_GmWmCsf(3, :) & ...  % GM > CSF
                    avgP_GmWmCsf(1, :) > 0.2;                        % GM > 0.2

    iexMask(2, :) = avgP_GmWmCsf(2, :) > avgP_GmWmCsf(1, :) & ...  % WM > GM
                    avgP_GmWmCsf(2, :) > avgP_GmWmCsf(3, :) & ...  % WM > CSF
                    avgP_GmWmCsf(2, :) > 0.2;                        % WM > 0.2

    iexMask(3, :) = avgP_GmWmCsf(3, :) > avgP_GmWmCsf(1, :) & ...  % CSF > GM
                    avgP_GmWmCsf(3, :) > avgP_GmWmCsf(2, :) & ...  % CSF > WM
                    avgP_GmWmCsf(3, :) > 0.2;                        % CSF > 0.2

    % Ajustement de la taille pour 2D et 3D
    if data_dims == 2
        iexMask = reshape(iexMask, [3, size(tissue_proba{1}, 2)]);  % 3 x Nx x Ny
    elseif data_dims == 3
        iexMask = reshape(iexMask, [3, size(tissue_proba{1}, 2), size(tissue_proba{1}, 3)]);  % 3 x Nx x Ny x Nz
    end
end
