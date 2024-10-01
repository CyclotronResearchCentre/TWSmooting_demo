function aj_quantify_ellipse_variability(ellipses_all)
    % Quantify variability between ellipses across different phantoms
    %
    % ellipses_all: a cell array where each cell contains the ellipses
    % parameters for a given phantom (Nx7 matrix for each phantom)
    
    num_phantoms = length(ellipses_all);
    
    if num_phantoms < 2
        error('At least two phantoms are required for comparison.');
    end
    
    % Initialize containers for metrics
    center_diff = zeros(num_phantoms, num_phantoms);
    size_diff = zeros(num_phantoms, num_phantoms);
    orientation_diff = zeros(num_phantoms, num_phantoms);
    
    % Loop over all pairs of phantoms
    for i = 1:num_phantoms
        for j = i+1:num_phantoms
            % Retrieve the two sets of ellipses
            ellipses1 = ellipses_all{i};
            ellipses2 = ellipses_all{j};
            
            % Ensure both phantoms have the same number of ellipses
            if size(ellipses1, 1) ~= size(ellipses2, 1)
                error('Phantoms must have the same number of ellipses for comparison.');
            end
            
            % Metric 1: Center difference (Euclidean distance)
            centers1 = ellipses1(:, 2:4); % Extract centers from ellipses1
            centers2 = ellipses2(:, 2:4); % Extract centers from ellipses2
            center_diff(i,j) = mean(sqrt(sum((centers1 - centers2).^2, 2))); % Euclidean distance
            center_diff(j,i) = center_diff(i,j); % Symmetry
            
            % Metric 2: Size difference (absolute difference of axes lengths)
            sizes1 = ellipses1(:, 5:7); % Extract radii from ellipses1
            sizes2 = ellipses2(:, 5:7); % Extract radii from ellipses2
            size_diff(i,j) = mean(abs(sizes1(:) - sizes2(:))); % Mean absolute difference
            size_diff(j,i) = size_diff(i,j); % Symmetry
            
            % Metric 3: Orientation difference (absolute difference in angles)
            angles1 = ellipses1(:, 8); % Extract orientation angles from ellipses1
            angles2 = ellipses2(:, 8); % Extract orientation angles from ellipses2
            orientation_diff(i,j) = mean(abs(angles1 - angles2)); % Mean absolute difference
            orientation_diff(j,i) = orientation_diff(i,j); % Symmetry
        end
    end
    
    % Display the results
    fprintf('Center Differences (Euclidean distance):\n');
    disp(center_diff);
    
    fprintf('Size Differences (Mean Absolute Difference):\n');
    disp(size_diff);
    
    fprintf('Orientation Differences (Mean Absolute Difference in angles):\n');
    disp(orientation_diff);
    
    % Plot heatmaps for each metric
    figure;
    subplot(1,3,1);
    imagesc(center_diff);
    colorbar;
    title('Center Difference (Euclidean Distance)');
    xlabel('Phantom Index');
    ylabel('Phantom Index');
    
    subplot(1,3,2);
    imagesc(size_diff);
    colorbar;
    title('Size Difference (Mean Absolute)');
    xlabel('Phantom Index');
    ylabel('Phantom Index');
    
    subplot(1,3,3);
    imagesc(orientation_diff);
    colorbar;
    title('Orientation Difference (Mean Absolute)');
    xlabel('Phantom Index');
    ylabel('Phantom Index');
end
