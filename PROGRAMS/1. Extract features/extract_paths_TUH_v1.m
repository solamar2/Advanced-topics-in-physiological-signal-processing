function [path, tests, date, montage_config, clinician_report_contents] = extract_paths_TUH_v1(patient_id, session,main_path)
%% Header
% extract_paths_TUH_v1 - Extracts .edf file paths, number of individual tests
% within a single recording session, recording date, montage type, and
% clinician report of the TUH recordings into the MATLAB environment.
%
% This function takes the input of patient_id and session to locate the
% corresponding data in the TUH dataset. It iterates through the
% directories, extracts relevant information, and stores the results
% in output variables.
%
% Inputs:
% - patient_id: ID of the patient for whom the data is being extracted.
% - session: Session identifier of the recording session.
%
% Outputs:
% - path: A cell array containing the paths of .edf files.
% - tests: A cell array containing the number of tests during the session.
% - date: A string containing the recording date.
% - montage_config: A string indicating the montage type used in the
% recordings.
% - clinician_report_contents: A string specifying the path of the clinician report.
%
% Author: Hamza Imtiaz
% Date: 2023-06-21

%% Initialize output variables
path = [];
tests = [];
date = [];

% Specify the parent directory of the TUH EEG dataset
tuh_dataset_dir = strcat(main_path,'\edf');


% Specify the subdirectories within the parent directory
% ar = average reference; le = linked ear
sub_dir = ["\01_tcp_ar", "\02_tcp_le", "\03_tcp_ar_a", "\04_tcp_le_a"];

% Specify the destination folder for the clinician report
destination_folder = 'D:\Users\Alon Friedman\Documents\Hamza\Analysis\Updated Temple Analysis\Clinician Reports';

% Set this to 1 if you'd like the clinician reports to be copied to
% the destination folder
copy_variable = 0;

%% Enter subdirectories
for z = 1:numel(sub_dir)
    found = 0; % Variable to track if patient_id is found

    % Format patient_id to be a string with 8 digits by adding leading zeros
    % based on the number of digits in patient_id
    num_of_digits = numel(num2str(patient_id));

    switch num_of_digits
        case 1
            patient_id_zeroes = append('0000000', string(patient_id));
        case 2
            patient_id_zeroes = append('000000', string(patient_id));
        case 3
            patient_id_zeroes = append('00000', string(patient_id));
        case 4
            patient_id_zeroes = append('0000', string(patient_id));
        case 5
            patient_id_zeroes = append('000', string(patient_id));
        case 6
            patient_id_zeroes = append('00', string(patient_id));
        case 7
            patient_id_zeroes = append('0', string(patient_id));
        case 8
            patient_id_zeroes = string(patient_id);
        otherwise
            error('Error: Unable to process patient ID')
    end

    % Extract the relevant portion from the patient_id_zeroes
    patient_id_first = patient_id_zeroes{1}(4:6);

    % Create the new path by joining the directory components
    new_path = join([tuh_dataset_dir, sub_dir(z), filesep, patient_id_first, filesep, patient_id_zeroes], "");

    % Read subdirectories and look for session (input parameter)
    % example: s005. Read subdirectories first 4 characters and match.
    % The remaining characters (6:end) are date in yyyy_mm_dd format.

    % Loop through each subfolder
    S = dir(new_path);
    N = setdiff({S([S.isdir]).name},{'.','..'}); % list of subfolders of S.

    if isempty(S)
        continue
    end

    % Loop through each subfolder
    for i = 1:numel(N)
        % Check if the subfolder matches the specified session
        if strcmp(session, N{i}(1:4))==1

            % Get the full path of the current subfolder
            final_dir = (fullfile(new_path,N{i}));

            % Retrieve the directory listing of EDF files in the subfolder

            tests_struc = dir(fullfile(final_dir,'*.edf'));

            %% Extract path, tests, and date
            % Loop through each EDF file
            tests = 0;
            for j = 1:numel(tests_struc)

                % Append the file name to the tests string
                tests = tests + 1;

                % Prepend the full path of the file to the path string
                path = [fullfile(final_dir, tests_struc(j).name) path];

                % Extract the date from the final_dir and append it to the date string
                pattern = '\d{4}_\d{2}_\d{2}';
                match = regexp(final_dir, pattern, 'match');
                date = string(match);
            end

            %% Find clinician reports and save as variable
            % Find the clinician report file
            clinician_report_files = dir(fullfile(final_dir, '*.txt'));

            % Check if any report files were found
            if ~isempty(clinician_report_files)
                % Extract the path of the first report file
                clinician_report_path = fullfile(clinician_report_files(1).folder, clinician_report_files(1).name);

                % Read the contents of the report file into a MATLAB variable
                clinician_report_contents = fileread(clinician_report_path);

                % Remove leading and trailing whitespace
                clinician_report_contents = strtrim(clinician_report_contents);

                % Uncomment and modify this section to remove special
                % characters
                %                 % Define the special characters to remove
                %                 specialChars = '[!@#$%^&*(),.?":{}|<>]';
                %
                %                 % Remove the special characters from the text
                %                 cleaned_text = regexprep(clinician_report_contents, specialChars, '');

                % Save the clinician report contents to a variable in consolidatedData
                consolidatedData.clinician_report_contents{i} = clinician_report_contents;
            else
                % Handle the case when no clinician report file is found
                consolidatedData.clinician_report_contents{i} = '';
            end

            % Save clinician reports to the destination folder if required
            %             if copy_variable == 1
            %                 copyfile(clinician_report,destination_folder);
            %             end

            % Set the montage configuration
            montage_config = sub_dir(z);

            % Set the found flag to 1
            found = 1;
        end
    end

    % Exit the loop once the .edf and clinician files are found
    if found == 1
        break
    end
end

%% For debugging
% Check if tests is empty
if isempty(tests)
    tests = 'No tests available';
end

% Check if path is empty
if isempty(path)
    path = 'No path available';
end

% Check if date is empty
if isempty(date)
    date = 'No date available';
end

% Check if montage_config exists
if ~exist('montage_config', 'var') || isempty(montage_config)
    montage_config = 'No montage config available';
end

% Check if clinician_report variable exists
if ~exist('clinician_report_contents', 'var') || isempty(clinician_report_contents)
    clinician_report_contents = 'No clinician report available';
end

end