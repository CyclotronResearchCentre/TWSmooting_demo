function [gsP_signal] = aj_smooth_gaussian_1D(data_1D, param)
    % Script pour appliquer un lissage gaussien standard sur des données 1D

    % Générer le noyau gaussien
    wg = gausswin(param.sm_kern_gaussian);
    wg = wg / sum(wg);  % Normalisation

    % Appliquer le lissage
    gsP_signal = filtfilt(wg, 1, data_1D)'; % transpose to get 1x128 vector

%     % Afficher les résultats
%     if flag.plot_fig
%         figure;
%         plot(data_1D, 'b-', 'DisplayName', 'Données brutes');
%         hold on;
%         plot(gsP_signal, 'r-', 'DisplayName', 'Lissage Gaussien');
%         legend;
%         xlabel('Position');
%         ylabel('Valeur du signal');
%         title('Lissage Gaussien sur Données 1D');
%         hold off;
%     end
% 
%     % Sauvegarder le signal lissé
%     if flag.save_data
%         save('smoothed_gaussian_1D.mat', 'gsP_signal');
%     end
end
