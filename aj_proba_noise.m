function noisy_proba_map = aj_proba_noise(proba_map, noise_range)
% Function to add noise to the tissue probability maps
%
% INPUT
% proba_map:    3D matrix containing the tissue probability maps
% noise_range:  Scalar value which multiplies the added noise
%
% OUTPUT
% noisy_proba_map:  3D matrix containing the noisy tissue probability maps
%
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------
%% Dealing with inputs
% Determine the number of dimensions (1D, 2D or 3D)
dims = size(proba_map);
if length(dims) == 4
    [Nx, Ny, Nz, ~] = size(proba_map);  % For 3D input
elseif length(dims) == 3
    [Nx, Ny, ~] = size(proba_map);      % For 2D input
    Nz = 1;
elseif length(dims) == 2
    Nx = size(proba_map, 1);            % For 1D input
    Ny = 1;
    Nz = 1;
else
    error('Invalid input dimensions. Expected 1D, 2D or 3D data.');
end

%% Do the job
% Initialize noisy_proba_map with the same dimensions
noisy_proba_map = zeros(size(proba_map));

% Iterate through each voxel and add noise to its probability vector
for x = 1:Nx
    for y = 1:Ny
        for z = 1:Nz
            % Extract the probability vector for the current voxel
            if Nz > 1
                prob_vector = squeeze(proba_map(x, y, z, :));
            elseif Ny > 1
                prob_vector = squeeze(proba_map(x, y, :));
            else
                prob_vector = squeeze(proba_map(x, :));
            end

            % Skip if the prob_vector is all zeros
            if all(prob_vector == 0)
                continue;
            end

            % Add noise to the probability vector
            noise = (rand(size(prob_vector)) - 0.5) * noise_range;

            % Add noise to the probability vector
            noisy_vector = prob_vector + noise;

            % Ensure all values are non-negative
            noisy_vector(noisy_vector < 0) = 0;

            % Normalize the noisy vector so that its sum equals 1
            noisy_vector = noisy_vector / sum(noisy_vector);

            % Assign the noisy vector back to the noisy_proba_map
            if Nz > 1
                noisy_proba_map(x, y, z, :) = noisy_vector;
            elseif Ny > 1
                noisy_proba_map(x, y, :) = noisy_vector;
            else
                noisy_proba_map(x, :) = noisy_vector;
            end
        end
    end
end
end
