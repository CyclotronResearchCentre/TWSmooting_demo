function all_files = aj_get_bids_files(N, image_type)
    % This function retrieves NIfTI files based on the number of subjects (N) and the image type (func, anat, dwi, fmap).
    % Arguments:
    %   N          : The number of subjects
    %   image_type : The type of image to retrieve ('func', 'anat', 'dwi', 'fmap')
    %
    % Returns:
    %   all_files  : A cell array of paths to the selected NIfTI files
    
    % Base directory paths
    base_path = pwd;
    data_path = fullfile(base_path, 'work_copies', 'openneuro.org', 'ds000117');
    
    % Automatically generate the list of subject IDs based on [1, N]
    subjects = arrayfun(@(x) sprintf('sub-%02d', x), 1:N, 'UniformOutput', false);
    disp('List of subjects:');
    disp(subjects);

    % Initialize an empty cell array to store file paths
    all_files = {};
    
    % Iterate over subjects
    for i = 1:length(subjects)
        subject_id = subjects{i};
        
        % Path to the images based on the specified type (BIDS format)
        image_dir = fullfile(data_path, subject_id, 'ses-mri', image_type);
        
        % Select NIfTI images (.nii files) in the specified directory
        f = spm_select('FPList', image_dir, '^.*\.nii$');
        
        % Check if files are found
        if isempty(f)
            fprintf('No %s images found for %s in %s\n', image_type, subject_id, image_dir);
        else
            fprintf('Found %s images for %s in %s\n', image_type, subject_id, image_dir);
            disp(f);  % Display the files found
            
            % Convert the character array to a cell array and store the files
            file_list = cellstr(f);
            all_files = [all_files; file_list]; % Concatenate the list of found files
        end
    end
    
    % Display available NIfTI files for user selection
    disp('Available NIfTI files for analysis:');
    for j = 1:length(all_files)
        [~, name, ext] = fileparts(all_files{j}); % Get file name and extension
        fprintf('%d: %s%s\n', j, name, ext); % Print file name and extension
    end
end
