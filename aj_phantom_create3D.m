function ph = aj_phantom_create3D(model, param)
% aj_create_phantom_3d Three-dimensional phantom
%   P = aj_create_phantom_3d(DEF, N) generates a 3D phantom with size N and default types.
%   DEF can be 'Shepp-Logan', 'Modified Shepp-Logan', 'Yu-Ye-Wang', or custom ellipsoids matrix.
%
% INPUTS:
%   model - String specifying phantom type: 'Shepp-Logan' or 'Modified
%         Shepp-Logan' or 'Yu-Ye-Wang'.
%   param - Structure containing default parameters such as the grid size.
% 
% OUTPUT:
%   ph      - Generated 3D phantom volume.
%
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------

grid_size = param.grid_size;
ph = zeros([grid_size, grid_size, grid_size], 'double');
rng = ((0:grid_size-1) - (grid_size-1)/2)/ ((grid_size-1)/2); % Normalize to have grid_size values between [-1 1]
[x, y, z] = meshgrid(rng, rng, rng);

% Flatten the grids
coord = [x(:), y(:), z(:)]';  % 3 x N matrix for voxel coordinates
ph = ph(:);  % Flatten the phantom

for k = 1:size(model, 1)
    % Get ellipsoid parameters
    A = model(k, 1);
    asq = model(k, 2)^2;            % square of the semi-axis along the x-axis
    bsq = model(k, 3)^2;            % square of the semi-axis along the y-axis
    csq = model(k, 4)^2;            % square of the semi-axis along the z-axis
    x0 = model(k, 5);               % center x-coordinate
    y0 = model(k, 6);               % center y-coordinate
    z0 = model(k, 7);               % center z-coordinate
    phi = deg2rad(model(k, 8));     % x rotation angle 
    theta = deg2rad(model(k, 9));   % y rotation angle 
    psi = deg2rad(model(k, 10));    % z rotation angle 

    % Euler rotation matrix
    alpha = euler_rotation(phi, theta, psi);

    % Apply rotation to voxel coordinates
    coordp = alpha * coord;

    % Find points inside the ellipsoid
    idx = ((coordp(1, :) - x0).^2 / asq) + ((coordp(2, :) - y0).^2 / bsq) + ((coordp(3, :) - z0).^2 / csq) <= 1;
    ph(idx) = ph(idx) + A;
end

ph = reshape(ph, [grid_size, grid_size, grid_size]);  % Reshape to original 3D volume

end

%% Helper function: computes Euler rotation matrix
function alpha = euler_rotation(phi, theta, psi)
    cphi = cos(phi); sphi = sin(phi);
    ctheta = cos(theta); stheta = sin(theta);
    cpsi = cos(psi); spsi = sin(psi);

    alpha = [cpsi*cphi - ctheta*sphi*spsi,  cpsi*sphi + ctheta*cphi*spsi,  spsi*stheta;
             -spsi*cphi - ctheta*sphi*cpsi, -spsi*sphi + ctheta*cphi*cpsi, cpsi*stheta;
             stheta*sphi,                  -stheta*cphi,                  ctheta];
end

