function PSWE_Main_TUH
%% Header
% PSWE_Main_TUH - This function performs a PSWE analysis on the TUH Corpus
%
% Inputs:
%   - Path to the excel document provided with the TUH corpus, which 
%       contains patient number and session along with some meta data
%
% Outputs:
%   - PSWE analysis results in an excel file
%
% Author: Hamza Imtiaz
% Date: 2023-06-21

%% Read me
% need to chenge path : 
% main_path - the main foled with all the codes and data

clear; clc; close all;
%% A) Extract information from TUH excel sheet
main_path='C:\Users\solam\OneDrive\Desktop\Master\Semester A\Advanced topics in physiological signal processing\final project';

excel_path = '\0. Find global PSWEs\relevant_data.xlsx';
excelFilePath = strcat(main_path,excel_path);

%for debug:
% excelFilePath = 'C:\Users\solam\OneDrive\Desktop\Thesis\TUH data\161 epilepsy v 192 non epilepsy - verified\check.xlsx';

consolidatedData = extract_info_TUH(excelFilePath,main_path);
filename = 'PSWE_TUH';

% Clear
clear excelFilePath
%% B) Select outpatient recordings
outpatientData = consolidatedData(strcmp(consolidatedData.EEGType, 'outpatient'), :);

outpatientData_ar = table();
count = 0;

for i = 1:height(outpatientData)
    if strcmp(outpatientData.montage_config{i}, '\01_tcp_ar')
        count = count + 1;
        outpatientData_ar(count, :) = outpatientData(i, :);
    end
end

%% For Sol data:
outpatientData_ar=outpatientData; 
% in case all the edf files are in the same folder without subfolders:
edf_path = strcat(main_path,'\edf');

for i = 1:size(outpatientData_ar,1)
    current_path=strcat(edf_path,'\',outpatientData_ar.edf_file_names(i));
    outpatientData_ar.path(i)= current_path;
end

clearvars -except outpatientData_ar main_path

%% Processing and feature extraction in batches

% Initialize variables
eegDataset_failed = [];     % Initialize the structure for failed datasets
eegDataset_preprocessed = [];  % Initialize the structure for preprocessed datasets
csv_output = []; % Variable for storing the output in a table

% Start EEGLAB
eeglab nogui


for yf = 1:height(outpatientData_ar)
    % 1)Import and preprocess
    %eegDataset_preprocessed_temp = import_edf(outpatientData_ar(yf, :),flag_preprocessing,flagICA);
    %eegDataset_preprocessed_temp = preprocessing_sol(outpatientData_ar(yf, :));

    eegDataset_preprocessed_temp = preprocessing_partial(outpatientData_ar(yf, :),main_path);

    % 2)Feature extraction
    eegDataset_final_temp = PSWE_v2(eegDataset_preprocessed_temp);

    % Remove data after feature extraction to conserve memory
    eegDataset_final_temp = rmfield(eegDataset_final_temp, 'data');

    % Check if processing is successful
    if strcmp(eegDataset_final_temp.preprocess_check, 'Preprocessed')
        if isempty(eegDataset_preprocessed)
            eegDataset_preprocessed = eegDataset_final_temp;
        else
            eegDataset_preprocessed(end + 1) = eegDataset_final_temp;
        end
    else
        if isempty(eegDataset_failed)
            eegDataset_failed = eegDataset_final_temp;
        else
            eegDataset_failed(end + 1) = eegDataset_final_temp;
        end
    end

    % Clear temporary variables
    clear eegDataset_preprocessed_temp eegDataset_final_temp;

end

%% 3) Export results

% Initialize the csv_output table
csv_output = [];

% Loop through each row in eegDataset_preprocessed
for i = 1:length(eegDataset_preprocessed)
    % Get the current row of eegDataset_preprocessed
    currentRow = eegDataset_preprocessed(i);
    
    % Call the export_PSWE function for the current row
    output = export_PSWE(currentRow);
    
    % Append the output to csv_output
    csv_output = [csv_output; output];
end

%% Save file
% Generate the CSV file with the specified filename and .csv extension
filename = 'PSWE_TUH';
% filename='check';

filenameWithExtension = strcat(filename, '.xlsx');
writetable(csv_output, filenameWithExtension);


%% 4. Sol: Save properties of global PSWE: for future studies
%{
presistently_slow=[];
no_pswe_no_presistently_slow=[];

for i=1:length(eegDataset_preprocessed)
    eegDataset = eegDataset_preprocessed(i);
    
    if eegDataset.PSWE.flag_PSWE==0
        % : 
        start_ind_GPSWE=eegDataset.PSWE.start_ind_global_pswe;
        end_ind_GPSWE=eegDataset.PSWE.end_ind_global_pswe;
        num_of_chan_GPSWE=cell2mat(eegDataset.PSWE.num_of_chan)';
        involved_chan_GPSWE=eegDataset.PSWE.involved_chan'; involved_chan_GPSWE= cell2mat(involved_chan_GPSWE);
        
        for t=1:length(start_ind_GPSWE)
             if(start_ind_GPSWE(t)==0)
                 start_ind_GPSWE(t)=[];
                 end_ind_GPSWE(t)=[];
            end
        end
   
        % in order to convert involved chan to table:
        % adding 0 because they aren't in the same size
        involved_mat=zeros(length(num_of_chan_GPSWE),max(num_of_chan_GPSWE));
        k=1;
        for t=1:length(num_of_chan_GPSWE)
            involved_mat(t,1:num_of_chan_GPSWE(t))=involved_chan_GPSWE(k:k+num_of_chan_GPSWE(t)-1);
            k=k+num_of_chan_GPSWE(t);
        end

        %saving global pswe data into xlsx file:
        GPSWE_Prop=table(start_ind_GPSWE,end_ind_GPSWE,num_of_chan_GPSWE,involved_mat,'VariableNames', {'Start_ind','End_ind','Num_of_channels','Involved_channels'});
        filenameWithExtension = strcat(eegDataset.edf_file_names  , '.xlsx');
        writetable(GPSWE_Prop, filenameWithExtension);
        
    elseif eegDataset.PSWE.flag_PSWE==1
         % if signal is persistently slow - all the variables are NaN 
         output = export_PSWE(eegDataset);
         % Append the output 
         presistently_slow= [presistently_slow; output];
         
    elseif eegDataset.PSWE.flag_PSWE==2
        % if signal has no PSWE and is not persistently slow - all the
        % variables are 0
         output = export_PSWE(eegDataset);
         % Append the output to csv_output
         no_pswe_no_presistently_slow= [no_pswe_no_presistently_slow; output];
    end
end

filename1 = 'presistently_slow.xlsx';
writetable(presistently_slow, filename1);

filename2 = 'no_pswe_no_presistently.xlsx';
writetable(no_pswe_no_presistently_slow, filename2);
%}
end