function [gsP_signal] = aj_smooth_gaussian(data, param)
    % Fonction pour appliquer un lissage gaussien standard sur des données
    % qui peuvent être en 1D, 2D, ou 3D.
    
    % Identifier la dimension des données
    if isvector(data)
        disp('Ceci est un vecteur.');
        num_dims = 1;
    elseif ismatrix(data)
        disp('Ceci est une matrice.');
        num_dims = 2;
    else
        data_size = size(data);
        num_dims = length(data_size);  % Nombre de dimensions des données (1, 2 ou 3)
        if num_dims == 3
            disp('Ceci est une matrice 3D.');
        end
    end

    % Générer le noyau gaussien
    wg = gausswin(param.sm_kern_gaussian);  % Noyau Gaussien pour le lissage
    wg = wg / sum(wg);  % Normalisation

    % Appliquer le lissage selon les dimensions des données
    switch num_dims
        case 1  % Lissage 1D
            gsP_signal = filtfilt(wg, 1, data)';  % transpose to get 1xN vector
        
        case 2  % Lissage 2D
            % Appliquer le lissage successivement sur chaque dimension
            gsP_signal = data;  % Initialiser avec les données d'origine
            for i = 1:2
                gsP_signal = imgaussfilt(gsP_signal, param.sm_kern_gaussian);  % Lissage Gaussien 2D
            end

        case 3  % Lissage 3D
            % Appliquer le lissage Gaussien 3D
            gsP_signal = imgaussfilt3(data, param.sm_kern_gaussian);  % Lissage Gaussien 3D

        otherwise
            error('La dimension des données n''est ni 1D, ni 2D, ni 3D');
    end
end
%     % Option d'affichage des résultats
%     if flag.plot_fig
%         switch num_dims
%             case 1
%                 figure;
%                 plot(data, 'b-', 'DisplayName', 'Données brutes');
%                 hold on;
%                 plot(gsP_signal, 'r-', 'DisplayName', 'Lissage Gaussien');
%                 legend;
%                 xlabel('Position');
%                 ylabel('Valeur du signal');
%                 title('Lissage Gaussien sur Données 1D');
%                 hold off;
% 
%             case 2
%                 figure;
%                 subplot(1, 2, 1);
%                 imagesc(data);
%                 title('Données brutes 2D');
%                 colorbar;
% 
%                 subplot(1, 2, 2);
%                 imagesc(gsP_signal);
%                 title('Lissage Gaussien 2D');
%                 colorbar;
% 
%             case 3
%                 figure;
%                 slice_idx = round(data_size(3) / 2);  % Afficher une coupe au milieu
%                 subplot(1, 2, 1);
%                 imagesc(data(:, :, slice_idx));
%                 title('Données brutes 3D (Coupe)');
%                 colorbar;
% 
%                 subplot(1, 2, 2);
%                 imagesc(gsP_signal(:, :, slice_idx));
%                 title('Lissage Gaussien 3D (Coupe)');
%                 colorbar;
%         end
%     end
% 
%     % Option de sauvegarde des résultats
%     if flag.save_data
%         output_filename = sprintf('smoothed_gaussian_%dD.mat', num_dims);
%         save(output_filename, 'gsP_signal');
%     end
