function [twsP_signal] = aj_smooth_TWS_1D(data_1D, data_GmWmCsfSculpt_1D, param)
% Script pour appliquer un lissage pondéré par les tissus sur des données 1D

    wg = gausswin(param.sm_kern_tws);
    wg = wg / sum(wg);  % Normalisation
    wg2 = gausswin(2 * param.sm_kern_tws);  % Noyau de lissage plus large
    wg2 = wg2 / sum(wg2);

    % Vérifier que data_1D est un vecteur ligne et le transposer si nécessaire
    if size(data_1D, 1) == 1
        data_1D = data_1D';
    end

    % Initialiser la variable pour le signal lissé
    twsP_signal = zeros(4, length(data_1D));

    % Appliquer le lissage pondéré par les tissus
    for ii = 1:4
        tmp1 = data_1D .* data_GmWmCsfSculpt_1D(:,ii) .* (filtfilt(wg2, 1, data_GmWmCsfSculpt_1D(:,ii)) > 0.05);
        twsP_signal(ii,:) = filtfilt(wg, 1, tmp1) ./ (filtfilt(wg, 1, data_GmWmCsfSculpt_1D(:,ii)) > 0.05);
    end

%     % Afficher les résultats
%     if flag.plot_fig
%         figure;
%         plot(data_1D, 'b-', 'DisplayName', 'Données brutes');
%         hold on;
%         plot(twsP_signal(1,:), 'r-', 'DisplayName', 'Lissage Pondéré (GM)');
%         plot(twsP_signal(2,:), 'g-', 'DisplayName', 'Lissage Pondéré (WM)');
%         plot(twsP_signal(3,:), 'Color', [1, 0, 1], 'DisplayName', 'Lissage Pondéré (CSF)');
%         plot(twsP_signal(4,:), 'Color', [1, 0.5, 0], 'DisplayName', 'Lissage Pondéré (SCULPT)');
%         legend;
%         xlabel('Position');
%         ylabel('Valeur du signal');
%         title('Lissage Pondéré par les Tissus sur Données 1D');
%         hold off;
%     end
% 
%     % Sauvegarder le signal lissé
%     if flag.save_data
%         save('smoothed_tissue_weighted_1D.mat', 'twsP_signal');
%     end
end