function consolidatedData = extract_info_TUH(excelFilePath,main_path)
%% Header
% extract_info_TUH - This function reads data stored in the
% _AAREADME_v01.xlsx excel file provided by the TUH corpus. It consolidates
% the data into a single sheet and provides information on data distribution.
% This function also extracts signal paths that can be input into further
% analysis. 
%
% Syntax: -extract_info_TUH(excelFilePath)
%
% Inputs:
%   - excelFilePath: Path to the excel file
%
% Outputs:
%   - consolidatedData: Table with the following columns:
%         consolidatedData.PatientId: TUH patient id
%         consolidatedData.Session: recording session
%         consolidatedData.EEGType: EEG type (EMU/ICU/Outpatient/Inpatient)
%         consolidatedData.EEGSubtype: subtype of EEG (ex. BURN/NICU/OR)
%         consolidatedData.LTM_Routine: LTM or routine recording
%         consolidatedData.path: signal path to .edf file
%         consolidatedData.tests: path to all tests within a single session
%         consolidatedData.date: recording date in YYYY_MM_DD format
%         consolidatedData.montage_config: either ar, le, ar_a, le_a, where
%         ar = average reference and le = linked ear
%         consolidatedData.clinician_report: clinician report for the .edf
%         file which may include demographic information, medical history,
%         neurologist impression and more. clinician reports are
%         unstructured and have varying information
%
% Example:
%   consolidatedData = extract_info_TUH(excelFilePath)
%
% Author: Hamza Imtiaz
% Date: 2023-06-21

%% Load and consolidate the excel file provided by TUH

% Specify the path to the Excel file (for debugging)
% excelFilePath = 'F:\tuh_eeg\v1.1.0\_AAREADME_v01.xlsx';

% Read the sheet names from the Excel file
sheetNames = sheetnames(excelFilePath);

% Initialize an empty table to store the consolidated data
consolidatedData = table();

% Iterate over each sheet
for sheetIndex = 1:numel(sheetNames)
    % Read the data from the current sheet
    sheetData = readtable(excelFilePath, 'Sheet', sheetNames{sheetIndex});

    % Append the current sheet's data to the consolidated table
    consolidatedData = [consolidatedData; sheetData];
end

%% Extract paths, number of tests, recording date, montage configuration
% and clinician report

% Loop through each row of consolidatedData
for i = 1:height(consolidatedData)
    % Get the patient_id and session for the current row
    patient_id = consolidatedData.PatientId(i);
    session = consolidatedData.Session{i};

    % Call the extract_paths_TUH function
    [path, tests, date, montage_config, clinician_report] = extract_paths_TUH_v1(patient_id, session,main_path);

    % Update the corresponding columns in consolidatedData
    consolidatedData.path{i} = path;
    consolidatedData.tests{i} = tests;
    consolidatedData.date{i} = date;
    consolidatedData.montage_config{i} = montage_config;
    consolidatedData.clinician_report{i} = clinician_report;
end

%% Save results to a .csv file

% % Specify the filename for the CSV file
% filename = 'TUH_Recordings_Metadata.csv';
% 
% % Save the table as a CSV file
% writetable(consolidatedData, filename, 'Delimiter', ',');
% 
% disp('consolidatedData saved as CSV successfully.');

%% ========================================================================
% Data description: The remaining code below is used to generate pie charts
% showing the distribution of the EEG types and subtypes within the TUH
% cohort. 
%% ========================================================================

% # of unique patients
% # of unique sessions per patient
% Counting unique patients and sessions
totalUniquePatients = numel(unique(consolidatedData.PatientId));
totalUniqueSessions = numel(unique(consolidatedData.Session));

% Counting LTM and routine recordings
numLTMRecordings = sum(consolidatedData.LTM_Routine == "LTM");
numRoutineRecordings = sum(consolidatedData.LTM_Routine == "Routine");

% Displaying the counts
disp(['Total unique patients: ', num2str(totalUniquePatients)]);
disp(['Total unique sessions: ', num2str(height(consolidatedData))]);
disp(['Number of LTM recordings: ', num2str(numLTMRecordings)]);
disp(['Number of routine recordings: ', num2str(numRoutineRecordings)]);

%% Pie charts of EEG types and subtypes
% Preprocess column
% Convert the 'EEGType' column to lowercase for case-insensitive matching
consolidatedData.EEGType = lower(consolidatedData.EEGType);

% Combine "outpatient" and "Outpatient" as the same category
consolidatedData.EEGType = replace(consolidatedData.EEGType, 'outpatient', 'outpatient');

% Combine "iCU" and "ICU" as the same category
consolidatedData.EEGType = replace(consolidatedData.EEGType, 'icu', 'icu');

% Remove "**> Cannot find info for" labels if encountered
consolidatedData = consolidatedData(~strcmp(consolidatedData.EEGType, '**> cannot find info for'), :);

% Count the occurrence of each EEG type
eegTypeCounts = countcats(categorical(consolidatedData.EEGType));

% Get the labels and counts for the pie chart
labels = categories(categorical(consolidatedData.EEGType));
counts = eegTypeCounts;

% Add count to labels
labelsWithCount = strcat(labels, ' (', string(counts), ')');

% Create the pie chart
pie(counts, labelsWithCount);

% Add a title to the pie chart
title('Distribution of EEG Types');

%% EMU subtypes
% Filter the data for the "emu" type
% Filter the data for the "emu" type
close all;
emuData = consolidatedData(consolidatedData.EEGType == "emu", :);

% Count the occurrence of each EEG subtype within the "emu" type
subtypeCounts = countcats(categorical(emuData.EEGSubtype));

%{
% Get the labels and counts for the pie chart
subtypeLabels = strcat(categories(categorical(emuData.EEGSubtype)), " (", string(subtypeCounts), ")");

% Create the pie chart
pie(subtypeCounts, subtypeLabels);
title('Distribution of EEG Subtypes for EMU');
%}

%% ICU subtypes
% Filter the data for the "ICU" type
close all;
icuData = consolidatedData(consolidatedData.EEGType == "icu", :);

% Preprocess the EEG subtypes
icuData.EEGSubtype = replace(icuData.EEGSubtype, 'CCU', 'CICU');
icuData.EEGSubtype = replace(icuData.EEGSubtype, {'RiCU', 'RICU'}, 'RiCU');

% Count the occurrence of each EEG subtype within the "ICU" type
subtypeCounts = countcats(categorical(icuData.EEGSubtype));
%{
% Get the labels and counts for the pie chart
subtypeLabels = strcat(categories(categorical(icuData.EEGSubtype)), " (", string(subtypeCounts), ")");

% Create the pie chart
pie(subtypeCounts, subtypeLabels);
title('Distribution of EEG Subtypes for ICU');
%}

%% Inpatient
% Filter the data for the "inpatient" type
close all;
inpatientData = consolidatedData(consolidatedData.EEGType == "inpatient", :);

% Preprocess the EEG subtypes if needed

% Count the occurrence of each EEG subtype within the "inpatient" type
subtypeCounts = countcats(categorical(inpatientData.EEGSubtype));

%{
% Get the labels and counts for the pie chart
subtypeLabels = strcat(categories(categorical(inpatientData.EEGSubtype)), " (", string(subtypeCounts), ")");

% Create the pie chart
figure;
pie(subtypeCounts, subtypeLabels);
title('Distribution of EEG Subtypes for Inpatient');
%}

%% Outpatient
% Filter the data for the "outpatient" type
close all;
outpatientData = consolidatedData(consolidatedData.EEGType == "outpatient", :);

% Preprocess the EEG subtypes if needed

% Count the occurrence of each EEG subtype within the "outpatient" type
subtypeCounts = countcats(categorical(outpatientData.EEGSubtype));

% Get the labels and counts for the pie chart
subtypeLabels = strcat(categories(categorical(outpatientData.EEGSubtype)), " (", string(subtypeCounts), ")");

% Create the pie chart
figure;
pie(subtypeCounts, subtypeLabels);
title('Distribution of EEG Subtypes for Outpatient');

%% Close
close all;

end
%%
