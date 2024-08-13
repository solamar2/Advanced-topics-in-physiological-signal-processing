function csv_output = export_PSWE(eegDataset)
%% Header
% export_PSWE - This function imports edf files into the MATLAB environment,
% picks channels and assigns channel locations
%
% Inputs:
%   - eegDataset: contains PSWE processed data
%
% Outputs:
%
% Author: Hamza Imtiaz
% Date: 2023-06-21

%% Patient identifier and recording information
patientID = {eegDataset.PatientId};
session = {eegDataset.Session};
recordingDate = {eegDataset.date};
numberOfTests = {eegDataset.tests};

% Extract the recording length information
recordingLength = {eegDataset.hdr.recording_length};

% Create a table with the extracted information
tableData_id = table(patientID', session', recordingDate', numberOfTests', recordingLength,...
    'VariableNames',{'Patient', 'Session', 'Recording_date', 'Number_of_tests',...
    'Recording_length'});

%% PSD
%{
delta_PSD = mean(eegDataset.PSD.BW.delta);
theta_PSD = mean(eegDataset.PSD.BW.theta);
alpha_PSD = mean(eegDataset.PSD.BW.alpha);
beta_PSD = mean(eegDataset.PSD.BW.beta);
gamma_PSD = mean(eegDataset.PSD.BW.gamma);

% Calculate the sum of all frequency bands
total_PSD = delta_PSD + theta_PSD + alpha_PSD + beta_PSD + gamma_PSD;

% Normalize the PSD values
delta_normalized = delta_PSD / total_PSD;
theta_normalized = theta_PSD / total_PSD;
alpha_normalized = alpha_PSD / total_PSD;
beta_normalized = beta_PSD / total_PSD;
gamma_normalized = gamma_PSD / total_PSD;

% Create a table with the normalized PSD values
tableData_PSD = table(delta_normalized, theta_normalized, alpha_normalized, beta_normalized,...
    gamma_normalized, ...
    'VariableNames', {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma'});
%}

%% Global PSWE
% Extract the relevant information from eegDataset
persistentSlowing = {eegDataset.PSWE.persistent_slowing}';
PSWEPerMin = {eegDataset.PSWE.PSWE_per_min}';
timeInPSWE = {eegDataset.PSWE.time_in_events}';
meanDuration = {eegDataset.PSWE.meanDur}';
medianDuration = {eegDataset.PSWE.medianDur}';

% Create the table
tableData_PSWE_global = table(persistentSlowing, PSWEPerMin, timeInPSWE, meanDuration, medianDuration, ...
    'VariableNames', {'Persistent_slowing', 'PSWE_per_min', 'Time_in_PSWE', 'Mean_duration', 'Median_duration'});
%{


%% Per channel PSWE
% Define the channel names
channelNames = {eegDataset.hdr.chanlocs.labels};

% --- PSWE per minute ---
% Initialize the table with channel names as variable names
tableData_PSWE_per_chan = table();

% Add columns for num_of_events_per_chan
for i = 1:numel(channelNames)
    channel = channelNames{i};
    channelValues = eegDataset.PSWE.num_of_events_per_chan(i);
    columnName = strcat(channel, '_PSWE_per_chan');
    tableData_PSWE_per_chan.(columnName) = channelValues;
end

% Add columns for time_in_events_per_chan
for i = 1:numel(channelNames)
    channel = channelNames{i};
    channelValues = eegDataset.PSWE.time_in_events_per_chan(i);
    columnName = strcat(channel, '_time_in_events');
    tableData_PSWE_per_chan.(columnName) = channelValues;
end

% Add columns for medianDur_per_chan
for i = 1:numel(channelNames)
    channel = channelNames{i};
    channelValues = eegDataset.PSWE.medianDur_per_chan(i);
    columnName = strcat(channel, '_median_duration');
    tableData_PSWE_per_chan.(columnName) = channelValues;
end

%% Clinician report
clinician_report = {eegDataset.clinician_report}';

% Create the table
tableData_PSWE_clinician_report = table(clinician_report, ...
    'VariableNames', {'Clinician_report'});
%}

%% Save and export to csv
% Concatenate the tables horizontally
csv_output = [tableData_id, tableData_PSWE_global];

%% Master
% disp('hZ')

end