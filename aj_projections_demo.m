function aj_projections_demo(phantom_3D, opt)
% Fonction pour effectuer des projections 1D et 2D à partir d'un phantom
% 3D bruité.
%
% Entrée:
%   - phantom_3D : Phantom 3D bruité (matrice 3D).
%   - opt        : Options pour la visualisation et le traitement (structure).
%       - opt.plot_all   : Afficher les figures multiples ou non (1 pour oui, 0 pour non, par défaut 0).
%       - opt.save_fig   : Sauvegarder les figures en PNG (1 pour oui, 0 pour non, par défaut 0).
%       - opt.fn_data    : Nom du fichier pour sauvegarder les résultats (facultatif).
%
%__________________________________________________________________________
% Copyright (C) 2024 Cyclotron Research Centre

%% 0/ Options
% Définir les options par défaut
opt_def = struct('plot_all', 1, 'save_fig', 0, 'fn_data', 'phantom_projections.mat');

% Vérifier les options fournies et les compléter avec les valeurs par défaut si nécessaire
if nargin < 2, opt = []; end
opt = crc_check_flag(opt_def, opt);

% Récupérer les options
plot_all  = opt.plot_all;
save_fig  = opt.save_fig;
fn_data   = opt.fn_data;

% Dimensions du phantom 3D
[dim_x, dim_y, dim_z] = size(phantom_3D);

%% 1/ Projections 1D et 2D
% Créer des projections 1D et 2D depuis le phantom 3D

% Projection 1D (somme sur l'axe Z)
projection_1D = squeeze(sum(phantom_3D, 3)); % Projeter le phantom sur un plan 2D (somme sur la 3ème dimension)

% Projection 2D (plan XY, YZ, et XZ)
projection_2D_XY = squeeze(sum(phantom_3D, 3)); % Projection sur le plan XY
projection_2D_XZ = squeeze(sum(phantom_3D, 2)); % Projection sur le plan XZ
projection_2D_YZ = squeeze(sum(phantom_3D, 1))'; % Projection sur le plan YZ

%% 2/ Visualisation et sauvegarde

% Affichage des projections 1D et 2D si plot_all est activé
if plot_all
    figure;
    
    % Projection 1D
    subplot(2,2,1);
    plot(projection_1D);
    title('Projection 1D (somme sur Z)');
    
    % Projection 2D XY
    subplot(2,2,2);
    imagesc(projection_2D_XY);
    axis image;
    title('Projection 2D (XY)');
    
    % Projection 2D XZ
    subplot(2,2,3);
    imagesc(projection_2D_XZ);
    axis image;
    title('Projection 2D (XZ)');
    
    % Projection 2D YZ
    subplot(2,2,4);
    imagesc(projection_2D_YZ);
    axis image;
    title('Projection 2D (YZ)');
    
    set(gcf,'Position',[500 150 1000 800]);
    
    % Sauvegarde des figures si save_fig est activé
    if save_fig
        saveas(gcf, 'phantom_projections.png');
    end
end

% Sauvegarde des résultats si demandé
save(fn_data, 'projection_1D', 'projection_2D_XY', 'projection_2D_XZ', 'projection_2D_YZ');

end

% Fonction pour vérifier et ajuster les options si elles ne sont pas spécifiées
function opt = crc_check_flag(opt_def, opt)
    fields = fieldnames(opt_def);
    for i = 1:length(fields)
        if ~isfield(opt, fields{i}) || isempty(opt.(fields{i}))
            opt.(fields{i}) = opt_def.(fields{i});
        end
    end
end

