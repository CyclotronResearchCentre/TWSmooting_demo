function file_paths = aj_BIDS_select(bids_info, varargin)
    % Parse input arguments
    p = inputParser;
    addParameter(p, 'sub', '', @ischar);
    addParameter(p, 'ses', '', @ischar);
    addParameter(p, 'modality', '', @ischar);
    addParameter(p, 'pattern', '', @ischar);
    addParameter(p, 'derivatives', '', @ischar);
    parse(p, varargin{:});
    
    subject = p.Results.sub;
    session = p.Results.ses;
    modality = p.Results.modality;
    pattern = p.Results.pattern;
    derivatives = p.Results.derivatives;

    % Get the subject info
    subject_info = bids_info.subjects(strcmp({bids_info.subjects.name}, subject));
    if isempty(subject_info)
        error('Subject %s not found in the BIDS structure.', subject);
    end
    
    % Validate session
    if ~isempty(session) && ~strcmp(subject_info.session, session)
        error('Session %s not found for subject %s.', session, subject);
    end
    
    % Check if the modality (func, anat, etc.) exists for this subject
    if ~isfield(subject_info, modality) || isempty(subject_info.(modality))
        error('Modality %s not found for subject %s.', modality, subject);
    end
    
    % Get the directory for the given modality in the derivatives folder
    modality_dir = fullfile(derivatives, subject, session, modality);
    if ~exist(modality_dir, 'dir')
        error('Directory %s does not exist.', modality_dir);
    end
    
%     % Display contents of the modality folder
%     fprintf('Searching in directory: %s\n', modality_dir);
%     dir_content = dir(modality_dir);
%     disp({dir_content.name});  % Display the list of files in the directory
    
    % Find files matching the pattern in the modality directory
    file_list = dir(fullfile(modality_dir, pattern));
    
    % Check if files matching the pattern were found
    if isempty(file_list)
        fprintf('No files matching pattern "%s" were found in directory %s\n', pattern, modality_dir);
    else
        fprintf('Found %d files matching pattern "%s" in directory %s\n', numel(file_list), pattern, modality_dir);
    end
    
    % Collect full file paths
    file_paths = fullfile({file_list.folder}, {file_list.name});
    
    % If no files found, return empty
    if isempty(file_paths)
        file_paths = {};
    end
end
