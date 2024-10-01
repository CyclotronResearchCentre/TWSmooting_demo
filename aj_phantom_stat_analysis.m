function aj_phantom_stat_analysis(phantoms_3D, param)
    % Function to perform statistical analysis on generated phantoms
    % Inputs:
    %   - phantoms_3D: Cell array containing the generated 3D phantoms
    %   - param: The parameter structure containing grid size, etc.
    
    num_phantoms = length(phantoms_3D);
    
    % Initialize arrays to store global variance and std for each phantom
    global_variances = zeros(num_phantoms, 1);
    global_stddevs = zeros(num_phantoms, 1);
    
    % Loop through each phantom and compute global variance and std dev
    for i = 1:num_phantoms
        current_phantom = phantoms_3D{i};
        
        % Compute global variance and standard deviation for the whole phantom
        global_variances(i) = var(current_phantom(:));
        global_stddevs(i) = std(current_phantom(:));
        
        fprintf('Phantom %d: Global variance = %.4f, Global std dev = %.4f\n', ...
                i, global_variances(i), global_stddevs(i));
    end
    
    % Now, compute mean and std dev across all phantoms
    mean_variance = mean(global_variances);
    mean_stddev = mean(global_stddevs);
    
    fprintf('Mean variance across all phantoms: %.4f\n', mean_variance);
    fprintf('Mean std dev across all phantoms: %.4f\n', mean_stddev);
    
    % Optional: Analyze variance per ellipse (if ellipses are available)
    if isfield(param, 'ellipsoids')
        num_ellipses = size(param.ellipsoids, 1);
        
        % Initialize containers for per-ellipse stats
        ellipse_variances = zeros(num_ellipses, num_phantoms);
        ellipse_stddevs = zeros(num_ellipses, num_phantoms);
        
        for i = 1:num_phantoms
            current_phantom = phantoms_3D{i};
            
            % Loop through each ellipse and compute variance/stddev
            for j = 1:num_ellipses
                % Assuming ellipses are stored in a logical mask (ellipses should be defined)
                % Example: mask = aj_extract_ellipse_mask(current_phantom, param.ellipsoids(j, :));
                % Compute variance and stddev inside the mask region
                mask = aj_extract_ellipse_mask(current_phantom, param.ellipsoids(j, :));
                ellipse_data = current_phantom(mask);
                
                ellipse_variances(j, i) = var(ellipse_data);
                ellipse_stddevs(j, i) = std(ellipse_data);
                
                fprintf('Phantom %d, Ellipse %d: Variance = %.4f, Std dev = %.4f\n', ...
                        i, j, ellipse_variances(j, i), ellipse_stddevs(j, i));
            end
        end
        
        % Compute mean and std dev per ellipse across all phantoms
        mean_ellipse_variance = mean(ellipse_variances, 2);
        mean_ellipse_stddev = mean(ellipse_stddevs, 2);
        
        for j = 1:num_ellipses
            fprintf('Ellipse %d: Mean variance = %.4f, Mean std dev = %.4f\n', ...
                    j, mean_ellipse_variance(j), mean_ellipse_stddev(j));
        end
    end
end
