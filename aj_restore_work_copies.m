function aj_restore_work_copies(subjects, orig_bids_root, work_bids_root)
%--------------------------------------------------------------------------
% Function to clear the work_copies folder by removing all derivatives 
% folders for a list of subjects and by restoring files from downloads for
% a list of subjects.
%
% INPUTS
% subjects:         a cell array which contains the list of subjects to be
%                   processed.
% orig_bids_root:   path to the download folder from which the files
%                   will be copied to restore the work_copies.
% work_bids_root:   path to the work_copies folder in which the
%                   subjects and derivatives (if exist) are.
%
% OUTPUT
% None:             the work_bids_root is a strict copy of orig_bids_root
%                   for specified subjects
%
% LIMITATIONS
% - All the data are bidsyfied (download folder, work folder and
%   derivatives folder).
% - SPM path has to be added.
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------

% Find all the different derivatives folders if they exist and take their paths in deriv_paths
disp('_____________________________________________________________________');
disp('Finding paths to all derivatives folders...');

if exist(work_bids_root, 'dir')
    work_dir = dir(work_bids_root);
    if ismember('derivatives',{work_dir.name})
        deriv_base_path = fullfile(work_bids_root, 'derivatives');
        if exist(deriv_base_path, 'dir')
            deriv_dir = dir(deriv_base_path);
            deriv_paths = cell(length(deriv_dir),1);
            for j = 1:length(deriv_dir)
                subdir_name = deriv_dir(j).name;
                if strcmp(subdir_name, '.') || strcmp(subdir_name, '..')
                    continue; % Skip '.' and '..'
                end
                deriv_paths{j} = fullfile(deriv_base_path, subdir_name);
                fprintf('Subdirectory called %s found in derivatives folder.\n', subdir_name);
            end
            emptyCells = cellfun(@isempty, deriv_paths); % delete empty cells (the ones corresponding to '.' and '..')
            deriv_paths(emptyCells) = [];
        else
            disp('No folder in the derivatives folder.');
        end
    else
        disp('No derivatives folder found.');
    end
else
    error('No folder found in the work directory (work_bids_root).');
end

disp('---------------------------------------------------------------------');

%--------------------------------------------------------------------------
% Look for subjects in the BIDS derivatives folder and remove their
% directory and subdirectories if they are contained in the input list
% subjects

disp('Removing derivatives folder for all subjects in the list...');

% Delete empty cells (if not all subjects of the database are selected)
emptyCells = cellfun(@isempty, subjects);
subjects(emptyCells) = [];

for i = 1:length(deriv_paths)
    deriv_BIDS = spm_BIDS(char(deriv_paths(i)));
    [~, lastWord, ~] = fileparts(deriv_paths(i));
    
    for ii = 1:length(subjects)
        if ismember(subjects{ii},{deriv_BIDS.subjects.name}) % using spm_BIDS(BIDS,'subjects') doesn't work: sub-01 and not 01
            deriv_sub_path = fullfile(deriv_BIDS.dir,subjects{ii});
            [SUCCESS,MESSAGE,MESSAGEID] = rmdir(deriv_sub_path, 's'); % directory and subdirectory tree will be removed recursively.
            if SUCCESS && isempty(MESSAGE)
                fprintf('Successfully deleted derivatives folder for subject %s in %s.\n', char(subjects{ii}), lastWord);
            else
                fprintf('Error deleting work folder for subject %s in %s: %s. MESSAGEID: %s.\n', char(subjects{ii}), lastWord, MESSAGE, MESSAGEID);
            end
        else
            fprintf('No folder found for subject %s in derivatives folder called %s.\n', char(subjects{ii}), lastWord);
        end
    end
end

disp('---------------------------------------------------------------------');

%--------------------------------------------------------------------------
% Look for subjects in the BIDS work folder and remove their
% directory and subdirectories if they are contained in the input list
% subjects

disp('Removing work folder for all subjects in the list...');

work_sub_path = cell(length(subjects),1);
work_BIDS = spm_BIDS(work_bids_root);
[~, lastWord, ~] = fileparts(work_bids_root);

for ii = 1:length(subjects)
    if ismember(subjects{ii},{work_BIDS.subjects.name}) % using spm_BIDS(BIDS,'subjects') doesn't work: sub-01 and not 01
        work_sub_path{ii} = fullfile(work_bids_root,subjects{ii});
        [SUCCESS,MESSAGE,MESSAGEID] = rmdir(work_sub_path{ii}, 's'); % directory and subdirectory tree will be removed recursively.
        if SUCCESS && isempty(MESSAGE)
            fprintf('Successfully deleted work folder for subject %s in %s.\n', char(subjects{ii}), lastWord);
        else
            fprintf('Error deleting work folder for subject %s in %s: %s. MESSAGEID: %s.\n', char(subjects{ii}), lastWord, MESSAGE, MESSAGEID);
        end
    else
        fprintf('No folder found for subject %s in work folder called %s.\n', char(subjects{ii}), lastWord);
    end
end

disp('---------------------------------------------------------------------');

%--------------------------------------------------------------------------
% Take the list of download paths for all subjects (orig_sub_path) from
% (orig_bids_root) if they exist

disp('Restoring work folder for all subjects in the list...');

orig_sub_path = cell(length(subjects),1);
if exist(orig_bids_root, 'dir')
    orig_dir = dir(orig_bids_root);
    for ii = 1:length(subjects)
        if ismember(subjects{ii},{orig_dir.name})
            orig_sub_path{ii} = fullfile(orig_bids_root,subjects{ii});
        else
            fprintf('No folder found for subject %s in download folder.\n', char(subjects{ii}));
        end
    end
else
    error('No folder found in the download directory (orig_bids_root).');
end

%--------------------------------------------------------------------------
% Copy the subject files from the download folder (orig_sub_path) to work
% folder (work_sub_path)
if eq(length(orig_sub_path),length(work_sub_path))
    for ii = 1:length(orig_sub_path)
        [~, orig_lastWord, ~] = fileparts(orig_sub_path{ii});
        [~, work_lastWord, ~] = fileparts(work_sub_path{ii});
        if eq(orig_lastWord,work_lastWord)
            [SUCCESS,MESSAGE,MESSAGEID] = copyfile(orig_sub_path{ii}, work_sub_path{ii}); % create DESTINATION if not exist as a directory and copy SOURCE to DESTINATION
            if SUCCESS && isempty(MESSAGE)
                fprintf('Successfully copied files from downloads to work_copies for %s.\n', char(subjects{ii}));
            else
                fprintf('Error copying files from downloads to work_copies for subject %s in %s: %s. MESSAGEID: %s.\n', char(subjects{ii}), lastWord, MESSAGE, MESSAGEID);
            end
        else
            fprintf('No match between removed subject in work folder %s and copied subject in download folder %s.\n', work_lastWord, orig_lastWord);
        end
    end
else
    fprintf('Number of removed subjects is not equal to the number of copied subjects.\n');
end

%--------------------------------------------------------------------------
% "Unarchive" archive files in the work folder for each subject in the
% input list
for ii = 1:length(work_sub_path)
    archive_files = dir(fullfile(work_sub_path{ii}, '**', '*.gz')); % get the list of compressed files for a specific subject ii
    if ~isempty(archive_files)
        % Process each file and unzip if it's a recognized archive type
        for k = 1:length(archive_files)
            file_path = fullfile(archive_files(k).folder, archive_files(k).name);
            try
                gunzip(file_path, archive_files(k).folder);
                delete(file_path);
            catch ME
                fprintf('Error unzipping GNU zip file (.gz) %s: %s\n', archive_files(k).name, ME.message);
            end
        end
        fprintf('Successfully unzipping GNU zip files (.gz) for subject %s\n', char(subjects{ii}))
    else
        fprintf('No archive file (.gz) for subject %s\n', char(subjects{ii}));
    end
end

disp('RESTORING: DONE');
disp('_____________________________________________________________________');

end
