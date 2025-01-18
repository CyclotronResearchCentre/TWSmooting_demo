% Main function to download BIDS files from the shell script and manage zip
% files
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------
% Name of the shell file
file_path = './ds000117-1.0.6.sh';

% Read the content of the shell file
file_content = fileread(file_path);

% Extract lines containing URLs
lines = strsplit(file_content, '\n');
image_urls = {};
image_types = {};

% Loop through each line and extract information about the images
for i = 1:length(lines)
    line = lines{i};
    if contains(line, 'curl') && contains(line, 'ses-mri')  % Ensure these are MRI files
        % Extract the URL and destination file
        tokens = regexp(line, 'https://\S+', 'match');
        if ~isempty(tokens)
            url = tokens{1};
            image_urls{end+1} = url;  % Store the URL

            % Identify the image type based on the subfolder in the URL
            if contains(url, 'ses-mri/anat')
                image_types{end+1} = 'Anatomical (anat)';
            elseif contains(url, 'ses-mri/dwi')
                image_types{end+1} = 'Diffusion (dwi)';
            elseif contains(url, 'ses-mri/fmap')
                image_types{end+1} = 'Magnetic Field Map (fmap)';
            elseif contains(url, 'ses-mri/func')
                image_types{end+1} = 'Functional (func)';
            else
                image_types{end+1} = 'Other';
            end
        end
    end
end

% Count the number of subjects based on 'sub-' in the URLs
subject_matches = regexp(image_urls, 'sub-\d+', 'match');
subject_ids = [subject_matches{:}];
unique_subjects = unique(subject_ids);
num_subjects = length(unique_subjects);
fprintf('Total number of subjects: %d\n', num_subjects);

% Display a selection menu to the user
unique_types = unique(image_types);
fprintf('Available image types for download:\n');
for i = 1:length(unique_types)
    fprintf('%d: %s\n', i, unique_types{i});
end
choice = input('Select the number corresponding to the image type to download: ');

% Check if the user made a valid choice
if choice > 0 && choice <= length(unique_types)
    selected_type = unique_types{choice};
    fprintf('Downloading images of type: %s\n', selected_type);

    % Download the selected images
    download_selected_images(image_urls, image_types, selected_type);
else
    fprintf('Invalid choice.\n');
end

%% HELP FUNCTIONS
% Function to download the selected images and manage zip copies
function download_selected_images(urls, types, selected_type)
    % Filter the URLs based on the selected type
    selected_urls = urls(strcmp(types, selected_type));

    % Total number of files to download
    num_files = length(selected_urls);

    % Check if a parallel pool is already running, and delete it if so
    pool = gcp('nocreate'); % Get current pool if it exists
    if ~isempty(pool)
        delete(pool); % Terminate the existing session
    end
    
    % Create a parallel pool
    parpool;  % Enable parallelism

    % Create a folder to store the working copies
    work_folder = 'work_copies';
    if ~exist(work_folder, 'dir')
        mkdir(work_folder);
    end

    % Download and process the files
    parfor i = 1:num_files
        url = selected_urls{i};
        file_path = urlsplit(url);  % Get the relative path from URL
        dest_file = fullfile('downloads', file_path);  % Original download path
        
        % Create necessary directories for the file
        dest_folder = fileparts(dest_file);
        if ~exist(dest_folder, 'dir')
            mkdir(dest_folder);
        end

        % Download the file if it does not already exist
        if ~exist(dest_file, 'file')
            try
                fprintf('Downloading %s...\n', url);
                websave(dest_file, url);  % Download the file
                fprintf('Successfully downloaded to %s\n', dest_file);
            catch ME
                fprintf('Error downloading %s: %s\n', url, ME.message);
            end
        else
            fprintf('File already exists: %s\n', dest_file);
        end

        % Copy the downloaded zip file to the working folder
        work_file = fullfile(work_folder, file_path);
        work_folder_path = fileparts(work_file);  % Create working folder structure
        if ~exist(work_folder_path, 'dir')
            mkdir(work_folder_path);
        end

        % Copy the zip file to the working directory
        copyfile(dest_file, work_file);
        fprintf('File copied to work directory: %s\n', work_file);

        % Unzip the copied file in the work folder
        if contains(work_file, '.gz')
            try
                fprintf('Unzipping %s...\n', work_file);
                gunzip(work_file, work_folder_path);
                fprintf('Successfully unzipped in work directory: %s\n', work_folder_path);
                
                % Delete the zip file after unzipping in work_folder
                delete(work_file);
                fprintf('Deleted the zip file: %s\n', work_file);
            catch ME
                fprintf('Error unzipping %s: %s\n', work_file, ME.message);
            end
        end
    end
    delete(pool); % Terminate the existing session
end

% Function to extract the relative path from a URL
function file_path = urlsplit(url)
    % Remove the parameters from the URL (everything after '?')
    base_url = strtok(url, '?');
    
    % Extract the part after the host (relative path)
    split_url = split(base_url, '/');
    
    % Use fullfile properly to join the elements
    file_path = fullfile(split_url{4:end}); % Ignore the first three elements (protocol + host)
end
