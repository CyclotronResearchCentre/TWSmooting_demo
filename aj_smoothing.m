function [gsP_signal, twsP_signal, f_twsP_signal, tosP_signal, f_tosP_signal] = ...
    aj_smoothing(ph, proba_map, param, flag, dim)
% Function to apply different smoothing techniques (Gaussian, TWS, TSPOON)
% to the given data, based on the flags provided.
%
% INPUT:
% ph:               Data to smooth                                      1D:[1 x nb_pt]
% proba_map:        Tissue probability maps matrix (for TWS and TSPOON) 1D:[nb_tissue x nb_pt]
% param:            Smoothing parameters
% flag:             Flags to determine which smoothing techniques to apply
% dim:              Dimension of the data (1D, 2D or 3D)
%
% OUTPUT:
% gsP_signal:       Gaussian smoothed signal                    1D:[1 x nb_pt]
% twsP_signal:      Tissue-weighted smoothed signal             1D:[nb_tissue x nb_pt]
% f_twsP_signal:    Final TWS signal after mask application     1D:[1 x nb_pt]
% tosP_signal:      TSPOON smoothed signal                      1D:[nb_tissue x nb_pt]

    % Apply standard Gaussian smoothing
    if flag.gaussian
        disp('Executing standard Gaussian smoothing...');
        gsP_signal = aj_smooth_gaussian(ph, param, dim);
    end

    % Apply tissue-weighted smoothing (TWS)
    if flag.tws
        disp('Executing tissue-weighted smoothing (TWS)...');
        [twsP_signal, f_twsP_signal] = aj_smooth_TWS_3D(ph, proba_map, param);
    end

    % Apply Tissue-Specific Smoothing (TSPOON)
    if flag.tspoon
        disp('Executing TSPOON smoothing...');
        [tosP_signal, f_tosP_signal] = aj_smooth_TSPOON(ph, proba_map, param, dim);
    end
end
