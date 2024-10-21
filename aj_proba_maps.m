function proba_map = aj_proba_maps(anat_phantom, param, flag)
    % Initialize the probability map matrix based on the dimension of the input
    data_size = size(anat_phantom);
    if isvector(anat_phantom)
        dims = 1;
        Nx = length(anat_phantom);
        proba_map = zeros(Nx, 4); % Adjusted to 1D with 4 channels for each voxel
    elseif ismatrix(anat_phantom)
        dims = 2;
        [Nx, Ny] = size(anat_phantom);
        proba_map = zeros(Nx, Ny, 4); % Adjusted to 2D with 4 channels for each voxel
    else
        dims = length(data_size);
        if dims ~= 3
            error('Input must be 1D, 2D or 3D');
        else 
            [Nx, Ny, Nz] = size(anat_phantom);
            proba_map = zeros(Nx, Ny, Nz, 4);
        end
    end
    
    % Fill the probability map for each type of tissue
    % Following SPM classification order : GmWmCsfSculpt
    switch dims
        case 1
            proba_map(:,4) = round(anat_phantom, 3) == 1;           % Bone Sculpt
            proba_map(:,3) = round(anat_phantom, 3) == 0.3;         % CSF
            proba_map(:,2) = round(anat_phantom, 3) == 0.2;         % WM
            proba_map(:,1) = round(anat_phantom, 3) == 0;           % GM
        case 2
            proba_map(:,:,4) = round(anat_phantom, 3) == 1;         % Bone Sculpt
            proba_map(:,:,3) = round(anat_phantom, 3) == 0.3;       % CSF
            proba_map(:,:,2) = round(anat_phantom, 3) == 0.2;       % WM
            proba_map(:,:,1) = round(anat_phantom, 3) == 0;         % GM
            
        case 3
            proba_map(:,:,:,4) = round(anat_phantom, 3) == 1;       % Bone Sculpt
            proba_map(:,:,:,3) = round(anat_phantom, 3) == 0.3;     % CSF
            proba_map(:,:,:,2) = round(anat_phantom, 3) == 0.2;     % WM
            proba_map(:,:,:,1) = round(anat_phantom, 3) == 0;       % GM
    end
    
    if flag.plot_fig
        % Convert voxel indices to mm
        voxel_size = param.voxreal_res;  % mm/voxel
        x_axis_mm = ((1:Nx) - Nx/2) * voxel_size;

        % Display the middle slice in color (for 3D) or the full map for 1D/2D
        switch dims
            case 1
                show_colored_1D(anat_phantom, x_axis_mm);
            case 2
                y_axis_mm = ((1:Ny) - Ny/2) * voxel_size;
                show_colored_slice(anat_phantom, x_axis_mm, y_axis_mm);
            case 3
                y_axis_mm = ((1:Ny) - Ny/2) * voxel_size;
                middle_slice = round(Nz / 2);
                slice = anat_phantom(:,:,middle_slice);
                show_colored_slice(slice, x_axis_mm, y_axis_mm);
        end
    end
end

%% Helper function to display the colored slice for 2D and 3D inputs
function show_colored_slice(slice, x_axis_mm, y_axis_mm)
    % Define colors for each tissue type
    bone_color = [1, 1, 1];   % White for bone
    csf_color = [0, 0, 1];    % Blue for CSF
    wm_color = [0, 1, 0];     % Green for white matter
    gm_color = [1, 0, 0];     % Red for gray matter

    % Create a colored image based on the tissue type
    Nx = length(x_axis_mm);
    Ny = length(y_axis_mm);
    colored_slice = zeros(Nx, Ny, 3);
    for i = 1:Nx
        for j = 1:Ny
            if round(slice(i,j), 3) == 1
                colored_slice(i,j,:) = bone_color;  % Bone
            elseif round(slice(i,j), 3) == 0.3
                colored_slice(i,j,:) = csf_color;   % CSF
            elseif round(slice(i,j), 3) == 0.2
                colored_slice(i,j,:) = wm_color;    % WM
            elseif round(slice(i,j), 3) == 0
                colored_slice(i,j,:) = gm_color;    % GM
            end
        end
    end

    % Display the colored slice with axes in mm
    figure;
    imagesc(x_axis_mm, y_axis_mm, colored_slice);
    title('Colored slice with different tissues');
    xlabel('X [mm]');
    ylabel('Y [mm]');
    axis image;
    set(gca, 'YDir', 'normal'); % imagesc: By default, displays the matrix so that the first row is at the top.
end

%% Helper function to display 1D colored data with gray intensity
function show_colored_1D(slice, x_axis_mm)
    % Define colors for each tissue type
    bone_color = [0, 0, 0];   % Black for bone
    csf_color = [0, 0, 1];    % Blue for CSF
    wm_color = [0, 1, 0];     % Green for white matter
    gm_color = [1, 0, 0];     % Red for gray matter

    % Create a colored representation for 1D
    Nx = length(x_axis_mm);
    colored_slice = zeros(Nx, 3);
    
    % Assign colors based on tissue types
    for i = 1:Nx
        if round(slice(i), 3) == 1
            colored_slice(i,:) = bone_color;  % Bone
        elseif round(slice(i), 3) == 0.3
            colored_slice(i,:) = csf_color;   % CSF
        elseif round(slice(i), 3) == 0.2
            colored_slice(i,:) = wm_color;    % WM
        elseif round(slice(i), 3) == 0
            colored_slice(i,:) = gm_color;    % GM
        end
    end

    % Get the gray intensity values
    gray_intensity = slice; % Assuming 'slice' represents intensity values

    % Display the 1D tissue representation with intensity on y-axis
    figure;
    hold on;  % Hold on to overlay the colored representation

    % Plotting the gray intensities along with colored markers
    for i = 1:Nx
        % Plot gray intensity as a line
        plot(x_axis_mm(i), gray_intensity(i), 'o', 'Color', colored_slice(i,:), 'MarkerFaceColor', colored_slice(i,:));
    end
    hold off;
    
    % Title and labels
    title('Colored 1D Tissue Data with Gray Intensity');
    xlabel('X [mm]');
    ylabel('Gray Intensity');
    axis tight;  % Adjust axes limits to fit the data
    ylim([0, 1]); % Adjust according to the expected range of gray intensity
end
