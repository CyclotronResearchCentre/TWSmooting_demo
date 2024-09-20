function [ph, ellipse] = aj_create_phantom_3d(varargin)
    % aj_create_phantom_3d Three-dimensional Shepp-Logan phantom for MRI simulation
    %   P = aj_create_phantom_3d(DEF, N) generates a 3D phantom with size N and default types.
    %   DEF can be 'Shepp-Logan', 'Modified Shepp-Logan', 'Yu-Ye-Wang', or custom ellipsoids matrix.
    %
    % INPUTS:
    %   DEF - String specifying phantom type: 'Shepp-Logan' or 'Modified
    %         Shepp-Logan' or 'Yu-Ye-Wang'. (default is 'Modified Shepp-Logan')
    %   N   - Scalar specifying the grid size (default is 128).
    % 
    % OUTPUT:
    %   P       - Generated 3D phantom volume.
    %   ELLIPSE - Matrix of ellipsoid parameters used to generate the phantom.
    %

    [ellipse, n] = parse_inputs(varargin{:});
    
    ph = zeros([n, n, n], 'double');  % Initialize 3D grid with zeros
    rng = ((0:n-1) - (n-1)/2) / ((n-1)/2);
    [x, y, z] = meshgrid(rng, rng, rng);
    
    % Flatten the grids
    coord = [x(:), y(:), z(:)]';  % 3 x N matrix for voxel coordinates
    ph = ph(:);  % Flatten the phantom
    
    for k = 1:size(ellipse, 1)
        % Get ellipsoid parameters
        A = ellipse(k, 1);
        asq = ellipse(k, 2)^2;
        bsq = ellipse(k, 3)^2;
        csq = ellipse(k, 4)^2;
        x0 = ellipse(k, 5);
        y0 = ellipse(k, 6);
        z0 = ellipse(k, 7);
        phi = deg2rad(ellipse(k, 8));
        theta = deg2rad(ellipse(k, 9));
        psi = deg2rad(ellipse(k, 10));
        
        % Euler rotation matrix
        alpha = euler_rotation(phi, theta, psi);
        
        % Apply rotation to voxel coordinates
        coordp = alpha * coord;
        
        % Find points inside the ellipsoid
        idx = ((coordp(1, :) - x0).^2 / asq) + ((coordp(2, :) - y0).^2 / bsq) + ((coordp(3, :) - z0).^2 / csq) <= 1;
        ph(idx) = ph(idx) + A;
    end
    
    ph = reshape(ph, [n, n, n]);  % Reshape to original 3D volume
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

% Helper function: parse input arguments
function [e, n] = parse_inputs(varargin)
    % Set defaults
    n = 128;  % Default grid size
    e = [];

    % Valid phantom types
    valid_phantoms = {'shepp-logan', 'modified shepp-logan', 'yu-ye-wang'};
    
    % Parse inputs
    for i = 1:nargin
        arg = varargin{i};
        if ischar(arg)
            def = lower(arg);
            idx = find(strncmp(def, valid_phantoms, length(def)), 1);
            if isempty(idx)
                error('Invalid phantom type: %s. Valid options are: %s.', def, strjoin(valid_phantoms, ', '));
            end
            switch valid_phantoms{idx}
                case 'shepp-logan'
                    e = shepp_logan();
                case 'modified shepp-logan'
                    e = modified_shepp_logan();
                case 'yu-ye-wang'
                    e = yu_ye_wang();
            end
        elseif isscalar(arg)
            n = arg;  % Grid size
        elseif ismatrix(arg) && size(arg, 2) == 10
            e = arg;  % Custom ellipsoid parameters
        else
            error('Invalid input argument format.');
        end
    end

    % Use modified Shepp-Logan phantom if no ellipsoid is specified
    if isempty(e)
        e = modified_shepp_logan();
    end
end

% Flatten a matrix into a row vector (for compatibility)
function out = flatten(in)
    out = in(:);
end

% Default ellipsoids for Shepp-Logan phantom
function e = shepp_logan()
    e = modified_shepp_logan();
    e(:, 1) = [1, -0.98, -0.02, -0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01];
end

% Default ellipsoids for Modified Shepp-Logan phantom
function e = modified_shepp_logan()
    e = [1, 0.69, 0.92, 0.81, 0, 0, 0, 0, 0, 0;
         -0.8, 0.6624, 0.874, 0.78, 0, -0.0184, 0, 0, 0, 0;
         -0.2, 0.11, 0.31, 0.22, 0.22, 0, 0, -18, 0, 10;
         -0.2, 0.16, 0.41, 0.28, -0.22, 0, 0, 18, 0, 10;
         0.1, 0.21, 0.25, 0.41, 0, 0.35, -0.15, 0, 0, 0;
         0.1, 0.046, 0.046, 0.05, 0, 0.1, 0.25, 0, 0, 0;
         0.1, 0.046, 0.046, 0.05, 0, -0.1, 0.25, 0, 0, 0;
         0.1, 0.046, 0.023, 0.05, -0.08, -0.605, 0, 0, 0, 0;
         0.1, 0.023, 0.023, 0.02, 0, -0.606, 0, 0, 0, 0;
         0.1, 0.023, 0.046, 0.02, 0.06, -0.605, 0, 0, 0, 0];
end

% Default ellipsoids for Yu-Ye-Wang phantom
function e = yu_ye_wang()
    e = [1, 0.69, 0.92, 0.9, 0, 0, 0, 0, 0, 0;
         -0.8, 0.6624, 0.874, 0.88, 0, 0, 0, 0, 0, 0;
         -0.2, 0.41, 0.16, 0.21, -0.22, 0, -0.25, 108, 0, 0;
         -0.2, 0.31, 0.11, 0.22, 0.22, 0, -0.25, 72, 0, 0;
         0.2, 0.25, 0.21, 0.31, 0, 0.35, -0.25, 0, 0, 0;
         0.2, 0.046, 0.046, 0.05, 0, 0.1, 0.625, 0, 0, 0;
         0.2, 0.046, 0.046, 0.05, 0, -0.1, 0.625, 0, 0, 0;
         0.2, 0.023, 0.046, 0.05, -0.08, -0.605, 0, 0, 0, 0;
         0.2, 0.023, 0.023, 0.02, 0, -0.606, 0, 0, 0, 0;
         0.2, 0.023, 0.046, 0.02, 0.06, -0.605, 0, 0, 0, 0];
end
