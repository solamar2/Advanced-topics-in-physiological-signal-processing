function eegDataset_import = preprocessing_partial(eegDataset,main_path)
%% Initialize variables
% Input paths for channel location layout
lookup_filepath_partial= '\0. Find global PSWEs\standard-10-5-cap385.elp';
lookup_filepath = strcat(main_path,lookup_filepath_partial);
channelloc_filepath_patial = '\0. Find global PSWEs\Standard-10-20-Cap19_custom_edf.ced';
channelloc_filepath = strcat(main_path,channelloc_filepath_patial);


% Preprocessing status
eegDataset.preprocess_check = "Unprocessed";

%% ====================================================================
% 1) Import module
% =====================================================================

% Specify the file paths of the .edf files in an array
filePath = char(eegDataset.path{1});

% Initialize the merged EEG object
EEG = [];

% Loop through the file paths
for j = 1:size(filePath,3)
    % Load the current file
    currentEEG = pop_biosig(filePath(:,:,j));

    % Check the duration of the current file
    durationMinutes = currentEEG.pnts / currentEEG.srate / 60; % Calculate duration in minutes

    % Only merge if the duration is at least 10 minutes
    if durationMinutes >= 10
        % Merge the current file with the existing EEG object
        if isempty(EEG)
            EEG = currentEEG;
        else
            EEG = pop_mergeset(EEG, currentEEG);
        end
    else
        %                 fprintf('Skipping file "%s" as it is less than 10 minutes long.\n', filePath(:,:,j));
    end
end

% Perform any necessary post-processing on the merged EEG object
EEG = eeg_checkset(EEG);

% If all tests are less than 10 minutes long,
if isempty(EEG)
    eegDataset_import = handlePreprocessFailure('Import', eegDataset);
    return
end

%% Select channels to include in the analysis
% List of electrodes from the EEG variable
electrode_labels = {EEG.chanlocs.labels}.';
electrode_labels_preremoval = electrode_labels;

% Remove unwanted substrings using case-insensitive regular expressions
patternsToRemove = {'^EEG', 'REF', 'LE', '-', 'AVE', ' '};
for i = 1:numel(patternsToRemove)
    pattern = patternsToRemove{i};
    electrode_labels = regexprep(electrode_labels, pattern, '', 'ignorecase');
end

% List of electrodes to select
selected_electrodes = {'Fp1', 'Fp2', 'F3', 'F4', 'C3', 'C4', 'P3', 'P4',...
    'O1', 'O2', 'F7', 'F8', 'T3', 'T4', 'T5', 'T6', 'Fz', 'Cz', 'Pz'};

selected_electrodes = char(selected_electrodes);  % Convert selected_electrodes to char array

% Convert electrode_labels and selected_electrodes to lowercase
electrode_labels = upper(char(string(electrode_labels)));
selected_electrodes = upper(selected_electrodes);

% Select specific electrodes from electrode_labels
[idx, ~] = ismember(electrode_labels, selected_electrodes, 'rows');
electrode_removed = electrode_labels(idx == 0, :); % for debugging
electrode_labels = electrode_labels(idx == 1, :);

% Select data corresponding to selected electrodes
try
    data_selected = EEG.data(idx == 1, :);
catch
    % Handle case when path is empty
    data_selected = [];
end

% Remove omitted electrodes from EEG.chanlocs
EEG.chanlocs(idx==0,:) = [];

% Update EEG properties
EEG.nbchan = length(electrode_labels);
EEG.data = data_selected;

% Assign processed electrode labels to EEG.chanlocs
for jj = 1:length(electrode_labels)
    EEG.chanlocs(jj).labels = electrode_labels(jj,:);
end

% Check the number of selected channels
if EEG.nbchan ~= 19
    eegDataset_import = handlePreprocessFailure('Channels', eegDataset);
    eegDataset_import.hdr.chanlocs = electrode_labels_preremoval';
    return;  % Exit the function or loop if desired
end

%% Import channel layout
% Load channel layout file
EEG = pop_chanedit(EEG, 'lookup', lookup_filepath, 'load', {channelloc_filepath, 'filetype', 'autodetect'});
EEG = eeg_checkset(EEG);

%% ====================================================================
% 2) Preprocessing module
% =====================================================================
% DC removal:
numElectrodes = size(EEG.data, 1); 
for i = 1:numElectrodes
    % Subtract the mean of each electrode's signal
    EEG.data(i, :) = EEG.data(i, :) - mean(EEG.data(i, :));
end

% Filter 1-45 Hz (EEGLAB default - Hamming windowed sinc FIR filter)
% See Widmann et. al (J Neurosci Methods. 2015)
try
    EEG = pop_eegfiltnew(EEG, 1, 45);

catch
    eegDataset_import = handlePreprocessFailure('Filter', eegDataset);
    return;  % Exit the function or loop if desired
end


%% Save results
% Convert to a structure
eegDataset_import = table2struct(eegDataset);

% Add recording length
eegDataset_import.hdr.recording_length = EEG.xmax / 60; % In minutes

% If preprocessing failed, remove subject from analysis
if ~strcmp(eegDataset.preprocess_check, "Unprocessed")
    eegDataset_import = handlePreprocessFailure('Other', eegDataset);
elseif strcmp(eegDataset.preprocess_check, "Unprocessed")
    eegDataset_import.preprocess_check = "Preprocessed";
    eegDataset_import.hdr.chanlocs = EEG.chanlocs;
    eegDataset_import.hdr.nbchan= EEG.nbchan;
    eegDataset_import.hdr.srate= EEG.srate;
    eegDataset_import.data = EEG.data;
end

end

% Subfunction to handle preprocess failure
function eegDataset_import = handlePreprocessFailure(failureMessage, eegDataset)

% Initialize output variable
eegDataset_import = table2struct(eegDataset);

% Update failure message
eegDataset_import.preprocess_check = failureMessage;

eegDataset_import.data = [];
eegDataset_import.hdr.chanlocs = [];
eegDataset_import.hdr.nbchan = [];
eegDataset_import.hdr.srate = [];
end