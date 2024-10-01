function [tosP_signal] = aj_smooth_TSPOON_1D(data_1D, data_GmWmCsfSculpt_1D, param)
% Script pour appliquer un lissage TSPOON sur des données 1D, incluant GM, WM, et CSF.

    % Vérifier que data_1D est un vecteur colonne et le transposer si nécessaire
    if size(data_1D, 2) == 1
        data_1D = data_1D';
    end

    wg = gausswin(param.sm_kern_tspoon);
    wg = wg / sum(wg);  % Normalisation
    
    % Générer le masque explicite à partir des probabilités tissulaires obtenue via le lissage gaussien standard
    gsP_GmWmCsfSculpt_1D = filtfilt(wg,1,data_GmWmCsfSculpt_1D)';
    iexMask = aj_create_explicit_mask_(gsP_GmWmCsfSculpt_1D);

    % Lissage du masque explicite
    gsP_iexMask = filtfilt(wg, 1, double(iexMask)')';

    % Appliquer le lissage TSPOON pour chaque tissu
    tosP_signal = zeros(size(iexMask));  % Taille ajustée selon iexMask
    for ii = 1:size(iexMask, 1)
        tosP_signal(ii, :) = filtfilt(wg, 1, iexMask(ii, :) .* data_1D) ./ gsP_iexMask(ii, :);
    end
    
%     % Afficher les résultats
%     if flag.plot_fig
%         figure;
%         plot(data_1D, 'b-', 'DisplayName', 'Données brutes');
%         hold on;
%         plot(tosP_signal(1,:), 'r-', 'DisplayName', 'Lissage TSPOON (GM)');
%         plot(tosP_signal(2,:), 'g-', 'DisplayName', 'Lissage TSPOON (WM)');
%         plot(tosP_signal(3,:), 'm-', 'DisplayName', 'Lissage TSPOON (CSF)');
%         legend;
%         xlabel('Position');
%         ylabel('Valeur du signal');
%         title('Lissage TSPOON sur Données 1D');
%         hold off;
%     end
% 
%     % Sauvegarder le signal lissé
%     if flag.save_data
%         save('smoothed_tspoon_1D.mat', 'tosP_signal');
%     end
end
