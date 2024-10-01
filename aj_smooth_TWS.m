function [twsP_signal] = aj_smooth_TWS(data, data_GmWmCsfSculpt, param)
    % Fonction pour appliquer un lissage pondéré par les tissus sur des données
    % qui peuvent être en 1D, 2D ou 3D (avec probabilités tissulaires GM, WM, CSF, Sculpt).
    
    % Identifier la dimension des données
    data_size = size(data);
    if isvector(data)
        disp('Ceci est un vecteur.');
        num_dims = 1;
    elseif ismatrix(data)
        disp('Ceci est une matrice.');
        num_dims = 2;
    else
        num_dims = length(data_size);  % Nombre de dimensions des données (1, 2 ou 3)
        if num_dims == 3
            disp('Ceci est une matrice 3D.');
        end
    end

    % Générer les noyaux gaussiens
    wg = gausswin(param.sm_kern_tws);
    wg = wg / sum(wg);  % Normalisation
    wg2 = gausswin(2 * param.sm_kern_tws);  % Noyau de lissage plus large
    wg2 = wg2 / sum(wg2);

    % Initialiser la variable pour le signal lissé
    twsP_signal = zeros([4, data_size]);  % 4 canaux : GM, WM, CSF, Sculpt

    % Appliquer le lissage pondéré par les tissus selon la dimension des données
    switch num_dims
        case 1  % Lissage 1D
            % Vérifier si les dimensions sont correctes
            if size(data_GmWmCsfSculpt, 1) ~= length(data)
                error('Les dimensions des données 1D et des probabilités tissulaires ne correspondent pas.');
            end
            % Vérifier que data_1D est un vecteur ligne et le transposer si nécessaire
            if size(data, 1) == 1
                data = data';
            end
            for ii = 1:4
                % Vérifier que les dimensions de `data` et `data_GmWmCsfSculpt` correspondent
                tmp1 = data .* data_GmWmCsfSculpt(:, ii) .* (filtfilt(wg2, 1, data_GmWmCsfSculpt(:, ii)) > 0.05);
                twsP_signal(ii, :) = filtfilt(wg, 1, tmp1) ./ (filtfilt(wg, 1, data_GmWmCsfSculpt(:, ii)) > 0.05);
            end
        
        case 2  % Lissage 2D
            for ii = 1:4
                tissue_mask = imgaussfilt(data_GmWmCsfSculpt(:,:,ii), param.sm_kern_tws) > 0.05;
                tmp1 = data .* data_GmWmCsfSculpt(:,:,ii) .* tissue_mask;
                twsP_signal(ii, :, :) = imgaussfilt(tmp1, param.sm_kern_tws) ./ imgaussfilt(data_GmWmCsfSculpt(:,:,ii), param.sm_kern_tws);
            end
        
        case 3  % Lissage 3D
            for ii = 1:4
                tissue_mask = imgaussfilt3(data_GmWmCsfSculpt(:,:,:,ii), param.sm_kern_tws) > 0.05;
                tmp1 = data .* data_GmWmCsfSculpt(:,:,:,ii) .* tissue_mask;
                twsP_signal(ii, :, :, :) = imgaussfilt3(tmp1, param.sm_kern_tws) ./ imgaussfilt3(data_GmWmCsfSculpt(:,:,:,ii), param.sm_kern_tws);
            end

        otherwise
            error('La dimension des données n''est ni 1D, ni 2D, ni 3D');
    end

%     % Option d'affichage des résultats
%     if flag.plot_fig
%         switch num_dims
%             case 1  % Affichage pour les données 1D
%                 figure;
%                 plot(data, 'b-', 'DisplayName', 'Données brutes');
%                 hold on;
%                 plot(twsP_signal(1,:), 'r-', 'DisplayName', 'Lissage Pondéré (GM)');
%                 plot(twsP_signal(2,:), 'g-', 'DisplayName', 'Lissage Pondéré (WM)');
%                 plot(twsP_signal(3,:), 'Color', [1, 0, 1], 'DisplayName', 'Lissage Pondéré (CSF)');
%                 plot(twsP_signal(4,:), 'Color', [1, 0.5, 0], 'DisplayName', 'Lissage Pondéré (SCULPT)');
%                 legend;
%                 xlabel('Position');
%                 ylabel('Valeur du signal');
%                 title('Lissage Pondéré par les Tissus sur Données 1D');
%                 hold off;
% 
%             case 2  % Affichage pour les données 2D
%                 figure;
%                 subplot(1, 2, 1);
%                 imagesc(data);
%                 title('Données brutes 2D');
%                 colorbar;
% 
%                 subplot(1, 2, 2);
%                 imagesc(twsP_signal(1,:,:));  % Affichage pour GM, par exemple
%                 title('Lissage Pondéré (GM) 2D');
%                 colorbar;
% 
%             case 3  % Affichage pour les données 3D
%                 slice_idx = round(data_size(3) / 2);  % Afficher une coupe au milieu
%                 figure;
%                 subplot(1, 2, 1);
%                 imagesc(data(:,:,slice_idx));
%                 title('Données brutes 3D (Coupe)');
%                 colorbar;
% 
%                 subplot(1, 2, 2);
%                 imagesc(twsP_signal(1,:,:,slice_idx));  % Affichage pour GM, par exemple
%                 title('Lissage Pondéré (GM) 3D (Coupe)');
%                 colorbar;
%         end
%     end
% 
%     % Option de sauvegarde des résultats
%     if flag.save_data
%         output_filename = sprintf('smoothed_tissue_weighted_%dD.mat', num_dims);
%         save(output_filename, 'twsP_signal');
%     end
end
