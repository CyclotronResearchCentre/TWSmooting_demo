function [iexMask] = aj_create_explicit_mask(tissue_proba, dim)
    % Function to create an explicit mask for different tissues
    % from the smoothed tissue probabilities in 1D, 2D, or 3D.
    %
    % INPUT:
    % tissue_proba: Matrix [nb_tissue x nb_pt_x ...] containing probabilities
    %               for different tissues in 1D, 2D, or 3D.
    %               Optionally, a cell array {nb_tissue x 1} with probabilities for multiple subjects.
    % dim: Dimension of the data (1 for 1D, 2 for 2D, 3 for 3D)
    %
    % OUTPUT:
    % iexMask: Explicit mask [nb_tissue x nb_pt_x ...] where the dimensions depend on the input `dim`.
    %          The mask contains nb_tissue slices for GM, WM, CSF, and Sculpt.
    %
    % Conditions:
    % - GM: GM probability > WM, CSF, and Sculpt probabilities, and probability > 0.2
    % - WM: WM probability > GM, CSF, and Sculpt probabilities, and probability > 0.2
    % - CSF: CSF probability > GM, WM, and Sculpt probabilities, and probability > 0.2
    % - Sculpt: Sculpt probability > GM, WM, and CSF probabilities, and probability > 0.2

    % Number of tissues
    nb_tissue = size(tissue_proba, 1);

    % Initialize avgP_GmWmCsf and iexMask based on dimension
    switch dim
        case 1  % 1D case
            nb_pt = size(tissue_proba, 2);  % Number of points in 1D
            avgP_GmWmCsf = tissue_proba;  % Tissue probabilities are already 1D
            iexMask = false(nb_tissue, nb_pt);  % Initialize the explicit mask

        case 2  % 2D case
            [nb_tissue, nb_pt_x, nb_pt_y] = size(tissue_proba);
            avgP_GmWmCsf = tissue_proba;  % Tissue probabilities are already 2D
            iexMask = false(nb_tissue, nb_pt_x, nb_pt_y);  % Initialize the explicit mask

        case 3  % 3D case
            [nb_tissue, nb_pt_x, nb_pt_y, nb_pt_z] = size(tissue_proba);
            avgP_GmWmCsf = tissue_proba;  % Tissue probabilities are already 3D
            iexMask = false(nb_tissue, nb_pt_x, nb_pt_y, nb_pt_z);  % Initialize the explicit mask

        otherwise
            error('Invalid dimension. dim must be 1, 2, or 3.');
    end

    % Handle the case where tissue_proba is a cell array for multiple subjects
    if iscell(tissue_proba)
        avgP_GmWmCsf = zeros(size(tissue_proba{1}));
        for ii = 1:nb_tissue
            avgP_GmWmCsf(ii, :, :, :) = mean(tissue_proba{ii}, dim + 1);  % Average over the last dimension (subjects)
        end
    end

    % Create the explicit mask for each tissue type
    for i = 1:nb_tissue
        % Start with a mask that is true for all points (1D, 2D, or 3D)
        switch dim
            case 1
                mask = true(1, nb_pt);  % Mask for 1D
            case 2
                mask = true(nb_pt_x, nb_pt_y);  % Mask for 2D
            case 3
                mask = true(nb_pt_x, nb_pt_y, nb_pt_z);  % Mask for 3D
        end

        % Compare tissue i with the other tissues
        for j = 1:nb_tissue
            if i ~= j
                switch dim
                    case 1
                        mask = mask & (avgP_GmWmCsf(i, :) > avgP_GmWmCsf(j, :));
                    case 2
                        mask = mask & (squeeze(avgP_GmWmCsf(i, :, :)) > squeeze(avgP_GmWmCsf(j, :, :)));
                    case 3
                        mask = mask & (squeeze(avgP_GmWmCsf(i, :, :, :)) > squeeze(avgP_GmWmCsf(j, :, :, :)));
                end
            end
        end

        % Apply the probability > 0.2 condition
        switch dim
            case 1
                mask = mask & (avgP_GmWmCsf(i, :) > 0.2);
            case 2
                mask = mask & (squeeze(avgP_GmWmCsf(i, :, :)) > 0.2);
            case 3
                mask = mask & (squeeze(avgP_GmWmCsf(i, :, :, :)) > 0.2);
        end

        % Assign the final mask for tissue i
        switch dim
            case 1
                iexMask(i, :) = mask;  % Correctly assign the mask for 1D
            case 2
                iexMask(i, :, :) = mask;  % Assign for 2D
            case 3
                iexMask(i, :, :, :) = mask;  % Assign for 3D
        end
    end
end

