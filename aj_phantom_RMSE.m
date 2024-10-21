function rmse_value = aj_phantom_RMSE(ph_1, ph_2)
% aj_phantom_RMSE Computes the Root Mean Square Error (RMSE) between two 3D phantoms
%
% INPUTS:
% ph_1    : First 3D phantom volume
% ph_2    : Second 3D phantom volume
%
% OUTPUT:
% rmse_value : Computed RMSE value

% Ensure the input volumes are of the same size
if size(ph_1) ~= size(ph_2)
    error('Input volumes must have the same dimensions.');
end

% Calculate RMSE
rmse_value = sqrt(mean((ph_1(:) - ph_2(:)).^2));

end
