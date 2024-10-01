function [param, flag] = aj_smooth_default()

% Parameters Definitions
param.sm_kern_gaussian = 4;  % Taille du noyau de lissage gaussien pour Gaussian
param.sm_kern_tws = 4;       % Taille du noyau de lissage gaussien pour TWS
param.sm_kern_tspoon = 4;    % Taille du noyau de lissage gaussien pour TSPOON

% fonction filtfilt (filtrage passe-bas par double filtrage) nécessite 
% que la longueur des données soit plus grande que trois fois la longueur du noyau de filtre

% Flag: Used to manage additional behaviors
flag.plot_fig = true;
flag.save_data = true;

end