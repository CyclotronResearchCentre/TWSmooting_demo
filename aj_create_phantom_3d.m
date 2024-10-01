function [ph, model] = aj_create_phantom_3d(model, FOV_size)
    % aj_create_phantom_3d Three-dimensional Shepp-Logan phantom for MRI simulation
    %   P = aj_create_phantom_3d(DEF, N) generates a 3D phantom with size N and default types.
    %   DEF can be 'Shepp-Logan', 'Modified Shepp-Logan', 'Yu-Ye-Wang', or custom ellipsoids matrix.
    %
    % INPUTS:
    %   model - String specifying phantom type: 'Shepp-Logan' or 'Modified
    %         Shepp-Logan' or 'Yu-Ye-Wang'.
    %   FOV_size   - Scalar specifying the grid size.
    % 
    % OUTPUT:
    %   P       - Generated 3D phantom volume.
    %   ELLIPSE - Matrix of ellipsoid parameters used to generate the phantom.
    %
    
    ph = zeros([FOV_size, FOV_size, FOV_size], 'double');  % Initialize 3D grid with zeros
    rng = ((0:FOV_size-1) - (FOV_size-1)/2) / ((FOV_size-1)/2);
    [x, y, z] = meshgrid(rng, rng, rng);
    
    % Flatten the grids
    coord = [x(:), y(:), z(:)]';  % 3 x N matrix for voxel coordinates
    ph = ph(:);  % Flatten the phantom
    
    for k = 1:size(model, 1)
        % Get ellipsoid parameters
        A = model(k, 1);
        asq = model(k, 2)^2;
        bsq = model(k, 3)^2;
        csq = model(k, 4)^2;
        x0 = model(k, 5);
        y0 = model(k, 6);
        z0 = model(k, 7);
        phi = deg2rad(model(k, 8));
        theta = deg2rad(model(k, 9));
        psi = deg2rad(model(k, 10));
        
        % Euler rotation matrix
        alpha = euler_rotation(phi, theta, psi);
        
        % Apply rotation to voxel coordinates
        coordp = alpha * coord;
        
        % Find points inside the ellipsoid
        idx = ((coordp(1, :) - x0).^2 / asq) + ((coordp(2, :) - y0).^2 / bsq) + ((coordp(3, :) - z0).^2 / csq) <= 1;
        ph(idx) = ph(idx) + A;
    end
    
    ph = reshape(ph, [FOV_size, FOV_size, FOV_size]);  % Reshape to original 3D volume
    
%     figure;
%     slice_num = round(size(ph, 1) / 2);
%     min_val = min(ph(:));
%     max_val = max(ph(:));
%     imagesc(squeeze(ph(:, :, slice_num)), [min_val, max_val]);
%     title('Original Phantom');
%     colormap gray;
%     axis image;    
end

% Helper function: computes Euler rotation matrix
function alpha = euler_rotation(phi, theta, psi)
    cphi = cos(phi); sphi = sin(phi);
    ctheta = cos(theta); stheta = sin(theta);
    cpsi = cos(psi); spsi = sin(psi);

    alpha = [cpsi*cphi - ctheta*sphi*spsi,  cpsi*sphi + ctheta*cphi*spsi,  spsi*stheta;
             -spsi*cphi - ctheta*sphi*cpsi, -spsi*sphi + ctheta*cphi*cpsi, cpsi*stheta;
             stheta*sphi,                  -stheta*cphi,                  ctheta];
end
