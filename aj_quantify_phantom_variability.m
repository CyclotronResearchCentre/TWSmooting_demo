function aj_quantify_phantom_variability(phantoms_3D)
    % Quantify anatomical variability between generated phantoms.
    %
    % This function calculates different metrics to quantify the variability
    % between the 3D phantoms generated, including voxel-wise differences,
    % volume differences, and correlation coefficients.
    
    num_phantoms = length(phantoms_3D);
    
    if num_phantoms < 2
        error('At least two phantoms are required for comparison.');
    end
    
    % Initialize containers for the metrics
    voxelwise_diff = zeros(num_phantoms, num_phantoms);
    volume_diff = zeros(num_phantoms, num_phantoms);
    correlation_coeff = zeros(num_phantoms, num_phantoms);
    IoU = zeros(num_phantoms, num_phantoms);
    
    % Loop over all pairs of phantoms
    for i = 1:num_phantoms
        for j = i+1:num_phantoms
            % Retrieve the two phantoms
            phantom1 = phantoms_3D{i};
            phantom2 = phantoms_3D{j};
            
            % Ensure both phantoms have the same size
            if any(size(phantom1) ~= size(phantom2))
                error('Phantoms must have the same dimensions for comparison.');
            end
            
            % Metric 1: Voxel-wise difference (mean absolute difference)
            voxelwise_diff(i,j) = mean(abs(phantom1(:) - phantom2(:)));
            voxelwise_diff(j,i) = voxelwise_diff(i,j); % Symmetry
            
            % Metric 2: Volume difference (sum of voxel intensities)
            volume1 = sum(phantom1(:));
            volume2 = sum(phantom2(:));
            volume_diff(i,j) = abs(volume1 - volume2);
            volume_diff(j,i) = volume_diff(i,j); % Symmetry
            
            % Metric 3: Correlation coefficient between voxel intensities
            correlation_coeff(i,j) = corr(phantom1(:), phantom2(:));
            correlation_coeff(j,i) = correlation_coeff(i,j); % Symmetry
            
            % Metric 4: Intersection over Union (IoU)
            % Define a threshold for binary comparison (e.g., 0.1)
            threshold = 0.1;
            binary1 = phantom1 > threshold;
            binary2 = phantom2 > threshold;
            intersection = sum(binary1(:) & binary2(:));
            union = sum(binary1(:) | binary2(:));
            IoU(i,j) = intersection / union;
            IoU(j,i) = IoU(i,j); % Symmetry
        end
    end
    
    % Display the results
    fprintf('Voxel-wise Differences:\n');
    disp(voxelwise_diff);
    
    fprintf('Volume Differences:\n');
    disp(volume_diff);
    
    fprintf('Correlation Coefficients:\n');
    disp(correlation_coeff);
    
    fprintf('Intersection over Union (IoU):\n');
    disp(IoU);
    
    % Plot heatmaps for each metric
    figure;
    subplot(2,2,1);
    imagesc(voxelwise_diff);
    colorbar;
    title('Voxel-wise Difference (Mean Absolute)');
    xlabel('Phantom Index');
    ylabel('Phantom Index');
    
    subplot(2,2,2);
    imagesc(volume_diff);
    colorbar;
    title('Volume Difference');
    xlabel('Phantom Index');
    ylabel('Phantom Index');
    
    subplot(2,2,3);
    imagesc(correlation_coeff);
    colorbar;
    title('Correlation Coefficient');
    xlabel('Phantom Index');
    ylabel('Phantom Index');
    
    subplot(2,2,4);
    imagesc(IoU);
    colorbar;
    title('Intersection over Union (IoU)');
    xlabel('Phantom Index');
    ylabel('Phantom Index');
end