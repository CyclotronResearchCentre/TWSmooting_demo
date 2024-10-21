function [param, flag] = aj_smooth_default()

% Parameters Definitions
param.sm_kern_gaussian = 3;  % Taille du noyau de lissage gaussien pour Gaussian
param.sm_kern_tws = 6;       % Taille du noyau de lissage gaussien pour TWS
param.sm_kern_tspoon = 6;    % Taille du noyau de lissage gaussien pour TSPOON

% fonction filtfilt (filtrage passe-bas par double filtrage) nécessite 
% que la longueur des données soit plus grande que trois fois la longueur du noyau de filtre

% Flag: Used to manage additional behaviors
flag.plot_fig = true;           % Flag to plot figures (set to true for plotting)
flag.save_data = true;          % Flag to save the results (set to true for saving)

flag.gaussian = true;           % Flag to compute gaussian smoothing
flag.tws = true;                % Flag to compute Tissue Weighted Smoothing
flag.tspoon = true;             % Flag to compute Tissue-SPecific smOOthing compeNsated

end