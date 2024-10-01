function [tosP_signal] = aj_smooth_TSPOON(data, data_GmWmCsfSculpt, param)
% Applique un lissage TSPOON sur des données 1D, 2D ou 3D

    % Identifier la dimension des données
    data_size = size(data);
    if isvector(data)
        disp('Ceci est un vecteur.');
        data_dims = 1;
    elseif ismatrix(data)
        disp('Ceci est une matrice.');
        data_dims = 2;
    else
        data_dims = length(data_size);  % Nombre de dimensions des données (1, 2 ou 3)
        disp('eeee');
        if data_dims == 3
            disp('Ceci est une matrice 3D.');
        end
    end

    % Sélectionner le noyau gaussien en fonction de la dimension
    wg = gausswin(param.sm_kern_tspoon);
    wg = wg / sum(wg);  % Normalisation

    % Générer le masque explicite à partir des probabilités tissulaires obtenues via le lissage gaussien standard
    gsP_GmWmCsfSculpt = filtfilt_gauss_proba(data_GmWmCsfSculpt, wg, data_dims);
    iexMask = aj_create_explicit_mask_1D(gsP_GmWmCsfSculpt);

    % Lissage du masque explicite
    gsP_iexMask = filtfilt_gauss_mask(double(iexMask), wg, data_dims);

    % Appliquer le lissage TSPOON pour chaque tissu
    tosP_signal = zeros(size(iexMask));

    % Appliquer un lissage adapté en fonction de la dimension
    if data_dims == 1  % 1D
        % Vérifier que data_1D est un vecteur colonne et le transposer si nécessaire
        if size(data, 2) == 1
            data = data';
        end
        for ii = 1:size(iexMask, 1)
            tosP_signal(ii, :) = filtfilt(wg, 1, iexMask(ii, :) .* data) ./ gsP_iexMask(ii, :);
        end
        
    elseif data_dims == 2  % 2D
        % Redimensionner data pour correspondre à la taille de iexMask si nécessaire
        if size(data, 1) ~= size(iexMask, 1) || size(data, 2) ~= size(iexMask, 2)
            data_resized = imresize(data, [size(iexMask, 1), size(iexMask, 2)]);
        else
            data_resized = data;
        end
        
        for ii = 1:size(iexMask, 3)  % Parcourir chaque plan 2D (chaque tissu)
            tosP_signal(:, :, ii) = filtfilt_gauss_mask(iexMask(:, :, ii) .* data_resized, wg, data_dims) ./ gsP_iexMask(:, :, ii);
            % tosP_signal(:, :, ii) = filtfilt_gauss_mask(iexMask(:, :, ii) .* data, wg, data_dims) ./ gsP_iexMask(:, :, ii);
        end
        
    elseif data_dims == 3  % 3D
        for ii = 1:size(iexMask, 4)  % Parcourir chaque volume 3D (chaque tissu)
            tosP_signal(:, :, :, ii) = filtfilt_gauss(iexMask(:, :, :, ii) .* data, wg, data_dims) ./ gsP_iexMask(:, :, :, ii);
        end
    else
        error('Dimensions non supportées.');
    end

%     % Afficher les résultats (optionnel)
%     if flag.plot_fig
%         visualize_results(data, tosP_signal, data_dims);
%     end
% 
%     % Sauvegarder les résultats
%     if flag.save_data
%         save(sprintf('smoothed_tspoon_%dD.mat', data_dims), 'tosP_signal');
%     end
end

% % Fonction pour appliquer filtfilt en fonction de la dimension (1D, 2D, 3D)
% function gsP_iexMask = filtfilt_gauss_proba(data, wg, data_dims)
%     if data_dims == 1  % 1D
%         gsP_iexMask = filtfilt(wg, 1, data)';
%     elseif data_dims == 2  % 2D
%         gsP_iexMask = imgaussfilt(data, wg);  % Lissage Gaussien 2D
%     elseif data_dims == 3  % 3D
%         gsP_iexMask = smooth3(data, 'gaussian', [wg wg wg]);  % Lissage Gaussien 3D
%     else
%         error('Dimensions non supportées pour filtfilt_gauss.');
%     end
% end

% Fonction pour appliquer filtfilt en fonction de la dimension (1D, 2D, 3D)
function gsP_iexMask = filtfilt_gauss_proba(data, wg, data_dims)
    if data_dims == 1  % 1D
        gsP_iexMask = filtfilt(wg, 1, data)';
    elseif data_dims == 2  % 2D
        sigma = std(wg);  % Utilisation de l'écart-type du noyau gaussien pour obtenir un sigma scalaire
        gsP_iexMask = imgaussfilt(data, sigma);  % Lissage Gaussien 2D avec un sigma scalaire
    elseif data_dims == 3  % 3D
        gsP_iexMask = smooth3(data, 'gaussian', [wg wg wg]);  % Lissage Gaussien 3D
    else
        error('Dimensions non supportées pour filtfilt_gauss_proba.');
    end
end

% % Fonction pour appliquer filtfilt en fonction de la dimension (1D, 2D, 3D)
% function gsP_iexMask = filtfilt_gauss_mask(data, wg, data_dims)
%     if data_dims == 1  % 1D
%         gsP_iexMask = filtfilt(wg, 1, data')';
%     elseif data_dims == 2  % 2D
%         gsP_iexMask = imgaussfilt(data, wg);  % Lissage Gaussien 2D
%     elseif data_dims == 3  % 3D
%         gsP_iexMask = smooth3(data, 'gaussian', [wg wg wg]);  % Lissage Gaussien 3D
%     else
%         error('Dimensions non supportées pour filtfilt_gauss.');
%     end
% end

% Fonction pour appliquer filtfilt en fonction de la dimension (1D, 2D, 3D)
function gsP_iexMask = filtfilt_gauss_mask(data, wg, data_dims)
    if data_dims == 1  % 1D
        gsP_iexMask = filtfilt(wg, 1, data')';
    elseif data_dims == 2  % 2D
        sigma = std(wg);  % Utilisation de l'écart-type pour obtenir un sigma scalaire
        gsP_iexMask = imgaussfilt(data, sigma);  % Lissage Gaussien 2D avec un sigma scalaire
    elseif data_dims == 3  % 3D
        gsP_iexMask = smooth3(data, 'gaussian', [wg wg wg]);  % Lissage Gaussien 3D
    else
        error('Dimensions non supportées pour filtfilt_gauss_mask.');
    end
end

