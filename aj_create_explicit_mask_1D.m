function [iexMask] = aj_create_explicit_mask_1D(tissue_proba)
    % Fonction pour créer un masque explicite pour la substance grise (GM),
    % la substance blanche (WM), et le CSF à partir des probabilités tissulaires lissées.
    %
    % INPUT :
    % tissue_proba : Matrice [4 x Npt] contenant les probabilités pour GM, WM, CSF et BoneSculpt,
    %                ou un cell array de {4 x 1} avec des probabilités pour plusieurs sujets.
    %
    % OUTPUT :
    % iexMask : Masque explicite pour GM, WM et CSF [3 x Npt].
    %           Le masque contient 3 lignes : GM, WM, et CSF.
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
            tissue_proba{ii} = tmp(ii,:);
        end
        need2avg = 0;  % Pas besoin de moyenne si un seul sujet
    else
        need2avg = 1;  % Moyenne requise si plusieurs sujets
    end

    % Nombre de points (Npt) dans les données
    Npt = size(tissue_proba{1}, 2);

    % Moyenne des probabilités tissulaires si plusieurs sujets
    if need2avg
        avgP_GmWmCsf = zeros(3, Npt);
        for ii = 1:3
            avgP_GmWmCsf(ii, :) = mean(tissue_proba{ii});
        end
    else
        avgP_GmWmCsf = zeros(3, Npt);
        for ii = 1:3
            avgP_GmWmCsf(ii, :) = tissue_proba{ii};
        end
    end

    % Calcul du masque explicite pour GM, WM et CSF
    iexMask = [ ...
        avgP_GmWmCsf(1,:) > avgP_GmWmCsf(2,:) & ...  % GM > WM
        avgP_GmWmCsf(1,:) > avgP_GmWmCsf(3,:) & ...  % GM > CSF
        avgP_GmWmCsf(1,:) > 0.2 ; ...               % GM > 0.2
        avgP_GmWmCsf(2,:) > avgP_GmWmCsf(1,:) & ...  % WM > GM
        avgP_GmWmCsf(2,:) > avgP_GmWmCsf(3,:) & ...  % WM > CSF
        avgP_GmWmCsf(2,:) > 0.2 ; ...               % WM > 0.2
        avgP_GmWmCsf(3,:) > avgP_GmWmCsf(1,:) & ...  % CSF > GM
        avgP_GmWmCsf(3,:) > avgP_GmWmCsf(2,:) & ...  % CSF > WM
        avgP_GmWmCsf(3,:) > 0.2 ];                  % CSF > 0.2
end
