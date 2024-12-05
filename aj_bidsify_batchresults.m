function aj_bidsify_batchresults(orig_bids_root, work_bids_root, subject)
%--------------------------------------------------------------------------
% Script to BIDSify files after running a MATLAB batch process (e.g., 
% preprocessing steps). The script identifies and organizes new files 
% created during preprocessing (typically named with prefixes like r*, w*,
% mean*) into a derivative directory according to the BIDS standards.
%
% INPUTS
% orig_bids_root:   path to the download folder from which the files
%                   will be copied to restore the work_copies.
% work_bids_root:   path to the work_copies folder in which the
%                   subjects and derivatives (if exist) are.
% subjects:         a cell array which contains the list of subjects to be
%                   processed.
%
% OUTPUT
% None:             the work_bids_root is bidsified for the specified 
%                   subjects thanks to comparing to orig_bids_root
%
% PROCESS
% 1. The script compares the original BIDS dataset (`orig_bids_root`) with 
%    a working copy of the BIDS dataset (`work_bids_root`) to identify 
%    newly created files.
% 2. New files are detected by comparing the base names (filename without 
%    extensions) between the two datasets for each type of data (e.g., 
%    `func`, `anat`, `fmap`, `dwi`).
% 3. Once new files are identified, they are moved into a 
%    `derivatives/preprocessing` folder, following the BIDS structure. The
%    `derivatives` folder is created if it doesn't exist.
% 4. Each data type (e.g., `anat`, `func`, etc.) has its own sub-folder 
%    under `derivatives/preprocessing`, ensuring a clear organization of 
%    the newly created files.
% 5. After all new files are moved, the script checks if the number of 
%    moved files matches the number of new files detected.
% 6. If the number matches (i.e., all files were successfully moved), the 
%    script deletes The corresponding files from the working BIDS dataset 
%    (`work_bids_root`), cleaning up the workspace.
%
% LIMITATIONS
% - If two files have the same name but different extensions (e.g., 
%   `file.nii` and `file.nii.gz`), the current method of detecting new files
%   (by comparing base names) may not treat them as distinct files. As a 
%   result, such files may not be correctly identified as new files and may 
%   remain in the working dataset.
%
% - This limitation exists because the comparison only looks at the file 
%   name without extensions. In future updates, this issue could be 
%   addressed by extending the comparison to include both base names and 
%   extensions.
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------

% Load the BIDS structure using spm_BIDS
bids_info = spm_BIDS(work_bids_root);

% Identify the specific subject
subject_info = bids_info.subjects(strcmp({bids_info.subjects.name}, subject));

% If the subject does not exist, display an error message
if isempty(subject_info)
    error('Subject %s not found in the BIDS structure.', subject);
end

% List of valid data types
valid_data_types = {'anat', 'func', 'fmap', 'dwi'};

% Iterate over the valid data types
for i = 1:numel(valid_data_types)
    data_type = valid_data_types{i};

    % Check if this data type contains files (non-empty struct)
    if ~isempty(subject_info.(data_type))
        % Get the original and work files for this data type
        orig_content = dir(fullfile(orig_bids_root, subject, subject_info.session, data_type));
        work_content = dir(fullfile(work_bids_root, subject, subject_info.session, data_type));

        % Filter to get only files (ignore directories)
        orig_files = orig_content(~[orig_content.isdir]);
        work_files = work_content(~[work_content.isdir]);

        % Extract base names of files (ignoring extensions)
        orig_names = cellfun(@(x) get_base_name(x), {orig_files.name}, 'UniformOutput', false);
        work_names = cellfun(@(x) get_base_name(x), {work_files.name}, 'UniformOutput', false);

        % Find new files
        new_files = setdiff(work_names, orig_names);
        fprintf('Found %d new file(s) for data type %s\n', size(new_files, 2), data_type);

        % If new files are found, copy or move them to the appropriate BIDS structure
        if ~isempty(new_files)
            successful_move = move_files_to_bids_structure(work_bids_root, subject, subject_info.session, data_type, new_files);

            % After confirming all files were successfully moved, delete from work directory
            if successful_move
                delete_new_files_in_work(work_bids_root, subject, subject_info.session, data_type, new_files);
            end
        end
    else
        fprintf('No file for data type %s \n', data_type);
    end
end

disp('BIDSIFYING: DONE');
disp('_____________________________________________________________________');

end

%% Function to get the base name by ignoring everything after the first "."
function base_name = get_base_name(name)
    base_name = strtok(name, '.'); % Ignore everything after the first "."
end

%% Function to copy/move files to the correct BIDS structure
function success = move_files_to_bids_structure(work_bids_root, subject, session, data_type, new_files)
% Path to the derivatives/preprocessing directory for this data type
target_dir = fullfile(work_bids_root, 'derivatives', 'preprocessing', subject, session, data_type);

% Create the directory if it doesn't exist
if ~exist(target_dir, 'dir')
    mkdir(target_dir);
end

% Variable to track if the move was successful
files_moved = 0;

% Loop over new files and copy/move them
for i = 1:numel(new_files)
    file_name = new_files{i};

    % Find the full original file (with extension) in work_bids_root
    src_file = dir(fullfile(work_bids_root, subject, session, data_type, [file_name, '*'])); % Uses '*' to match varied extensions
    if ~isempty(src_file)
        % Full path to the source file
        src_path = fullfile(src_file(1).folder, src_file(1).name);

        % Full path to the target file
        target_file = fullfile(target_dir, src_file(1).name);

        % Copy the file
        copyfile(src_path, target_file);

        % Display the result
%         fprintf('File %s copied to %s\n', src_path, target_file);

        % Count the file as moved
        files_moved = files_moved + 1;
    else
        fprintf('Empty fullfile src_file %s \n', src_file);
    end
end

% Check if all files were moved
if files_moved == numel(new_files)
    success = true;
    fprintf('The %d %s files have been moved.', size(new_files, 2), data_type);
else
    success = false;
    fprintf('Error: Not all files were successfully moved for %s.\n', data_type);
end
end

%% Function to delete files from the work directory after they have been moved
function delete_new_files_in_work(work_bids_root, subject, session, data_type, new_files)
    % Loop over the new files to delete them
    for i = 1:numel(new_files)
        file_name = new_files{i};
        
        % Find the full path to the file in the work directory
        src_file = dir(fullfile(work_bids_root, subject, session, data_type, [file_name, '*']));
        if ~isempty(src_file)
            % Full path to the source file
            src_path = fullfile(src_file(1).folder, src_file(1).name);
            
            % Delete the file
            delete(src_path);
            
            % Display the result
%             fprintf('File %s deleted from %s\n', src_path, fullfile(work_bids_root, subject, session, data_type));
        end
    end
end
