% Function to clear the work_copies folder and restore files from downloads
function aj_restore_work_copies()

clc; clear;
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');

downloads_path='C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\downloads\openneuro.org\ds000117\sub-02\ses-mri\anat';
work_copies_path='C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\work_copies\openneuro.org\ds000117\sub-02\ses-mri\anat'; 

    % Check if the work_copies folder exists
    if ~exist(work_copies_path, 'dir')
        fprintf('The specified work_copies folder does not exist: %s\n', work_copies_path);
        return;
    end
    
    % Clear all files and subdirectories from work_copies
    try
        rmdir(work_copies_path, 's'); % Delete the folder and its contents
        fprintf('All files in work_copies have been deleted.\n');
    catch ME
        fprintf('Error deleting files in work_copies: %s\n', ME.message);
        return; % Exit if deletion fails
    end
    
    % Recreate the work_copies folder
    mkdir(work_copies_path);
    fprintf('Recreated work_copies folder.\n');

    % Check if the downloads folder exists
    if ~exist(downloads_path, 'dir')
        fprintf('The specified downloads folder does not exist: %s\n', downloads_path);
        return;
    end
    
    % Copy files from downloads to work_copies
    try
        copyfile(fullfile(downloads_path, '*'), work_copies_path); % Copy all files
        fprintf('Successfully copied files from downloads to work_copies.\n');
    catch ME
        fprintf('Error copying files: %s\n', ME.message);
        return;
    end

    % Get the list of compressed files in work_copies_path (both .zip and other types)
    archive_files = dir(fullfile(work_copies_path, '*.*'));
    
    % Process each file and unzip if it's a recognized archive type
    for k = 1:length(archive_files)
        [~, ~, ext] = fileparts(archive_files(k).name);
        file_path = fullfile(work_copies_path, archive_files(k).name);
        
        switch lower(ext)
            case '.zip'
                try
                    unzip(file_path, work_copies_path);
                    fprintf('Unzipped file: %s\n', archive_files(k).name);
                    delete(file_path); % Delete after unzipping
                    fprintf('Deleted zip file: %s\n', archive_files(k).name);
                catch ME
                    fprintf('Error unzipping file %s: %s\n', archive_files(k).name, ME.message);
                end
                
            case {'.tar', '.tar.gz'}
                try
                    untar(file_path, work_copies_path);
                    fprintf('Untarred file: %s\n', archive_files(k).name);
                    delete(file_path); % Delete after unzipping
                    fprintf('Deleted tar file: %s\n', archive_files(k).name);
                catch ME
                    fprintf('Error untarring file %s: %s\n', archive_files(k).name, ME.message);
                end

            case '.rar'
                try
                    % Use system call to unrar (requires 7-Zip or WinRAR installed)
                    system(['7z x "' file_path '" -o"' work_copies_path '"']);
                    fprintf('Unrarred file: %s\n', archive_files(k).name);
                    delete(file_path); % Delete after unzipping
                    fprintf('Deleted rar file: %s\n', archive_files(k).name);
                catch ME
                    fprintf('Error un-rarring file %s: %s\n', archive_files(k).name, ME.message);
                end
                
            case '.gz'
                gunzip(file_path, work_copies_path);
                fprintf('Successfully unzipped in work directory: %s\n', work_copies_path);
                delete(file_path);
                fprintf('Deleted the zip file: %s\n', archive_files(k).name);

            % Add more cases for other archive types if needed

            otherwise
                fprintf('Skipping non-archive file: %s\n', archive_files(k).name);
        end
    end

    fprintf('All archives unzipped and copied successfully.\n');
end
