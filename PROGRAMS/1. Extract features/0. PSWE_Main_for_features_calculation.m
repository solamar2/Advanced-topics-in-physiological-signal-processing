function PSWE_Main_for_features_calculation
% After detecting PSWE and Global PSWE stored in xlsx files (using Hamza
% codes - PSWE_Main_TUH )
% this code uses the stat and end time of each global PSWE and calculate
% the features. 
% features calculation done using PSWE_fetures function

% Author: Sol Amara
% 25.01.2024

%% Read me
% need to chenge path : 
% main_path - the main foled with all the codes and data

clear; clc; close all;
main_path='C:\Users\solam\OneDrive\Desktop\Master\Semester A\Advanced topics in physiological signal processing\final project';

%% extract PSWE characteristics table from Verifed table:
%{
featre_table_filepath_partial = '\‏‏PSWE_TUH_characteristic.xlsx';
featre_table_filepath = strcat(main_path,featre_table_filepath_partial);
Verified_data_filepath_partial='\161 epilepsy v 192 non epilepsy - verified\161 nonepilepsy vs 188 epilepsy vs 715 NLP -V2 for Sol.xlsx';
Verified_data_filepath = strcat(main_path,Verified_data_filepath_partial);

featre_table=readtable(featre_table_filepath); featre_table = sortrows(featre_table, 'Patient');
Verified_data=readtable(Verified_data_filepath); Verified_data = sortrows(Verified_data, 'Patient_ID');

for i=1:length(featre_table.Patient)
    ind_verified=find(Verified_data.Patient_ID==featre_table.Patient(i));
    
   featre_table.Epilepsy_verified(i)=Verified_data.Epilepsy_verified(ind_verified);
   featre_table.Age(i)=Verified_data.Age(ind_verified);
   featre_table.Sex(i)=Verified_data.Sex(ind_verified);
end

writetable(featre_table,featre_table_filepath);
%}

%% extract PSWE channels characteristics:
%{
featre_table_filepath_partial = '\‏‏PSWE_TUH_characteristic.xlsx';
featre_table_filepath = strcat(main_path,featre_table_filepath_partial);

featre_table=readtable(featre_table_filepath); 

global_pswe_per_patient_filepath='C:\Users\solam\OneDrive\Desktop\Master\Semester A\Advanced topics in physiological signal processing\final project\Results\per patient';
% Get a list of all files in the folder
files = dir(fullfile(global_pswe_per_patient_filepath, '*.xlsx'));
% Extract and display the names of XLSX files
xlsxFileNames = {files.name};
xlsxFileNames_patient=regexprep(xlsxFileNames, '^0+|_.*', ''); 
xlsxFileNames_patient = str2double(xlsxFileNames_patient);


for i=1:length(featre_table.Patient)
    if featre_table.Persistent_slowing(i)==0 % if there is PSWE and it isn't presisitenly slow
    ind=find(xlsxFileNames_patient==featre_table.Patient(i));
    filenameWithExtension = strcat(global_pswe_per_patient_filepath,'\',xlsxFileNames(ind));
    PSWE_properties=readtable(strjoin(filenameWithExtension));

    featre_table.num_of_PSWE(i)=size(PSWE_properties,1);
    featre_table.mean_num_of_chan(i)=mean(PSWE_properties.Num_of_channels);
    featre_table.min_num_of_chan(i)=min(PSWE_properties.Num_of_channels);
    featre_table.max_num_of_chan(i)=max(PSWE_properties.Num_of_channels);
    end
end

writetable(featre_table,featre_table_filepath);

%}

%% extract EEG data:
% A) Extract information from TUH excel sheet
excelFilePath_partial = '\0. Find global PSWEs\relevant_data.xlsx';
excelFilePath = strcat(main_path,excelFilePath_partial);

%for debug:
% excelFilePath = 'C:\Users\solam\OneDrive\Desktop\Thesis\TUH data\161 epilepsy v 192 non epilepsy - verified\check.xlsx';

consolidatedData = extract_info_TUH(excelFilePath,main_path);
filename = 'PSWE_TUH';
% Clear
clear excelFilePath
% B) Select outpatient recordings 
outpatientData = consolidatedData(strcmp(consolidatedData.EEGType, 'outpatient'), :);
% all the patients are AR in our data. For other data set:
%{ 
outpatientData_ar = table();
count = 0;
for i = 1:height(outpatientData)
    if strcmp(outpatientData.montage_config{i}, '\01_tcp_ar')
        count = count + 1;
        outpatientData_ar(count, :) = outpatientData(i, :);
    end
end
outpatientData=outpatientData_ar; 
%}

% in case all the edf files are in the same folder without subfolders:
edf_path = strcat(main_path,'\edf');
for i = 1:size(outpatientData,1)
    current_path=strcat(edf_path,'\',outpatientData.edf_file_names(i));
    outpatientData.path(i)= current_path;
end
clearvars -except outpatientData main_path

%% extract Global PSWE characteristics table data:
global_pswe_per_patient_filepath_partial='\Results\per patient';
global_pswe_per_patient_filepath = strcat(main_path,global_pswe_per_patient_filepath_partial);

% Get a list of all files in the folder
files = dir(fullfile(global_pswe_per_patient_filepath, '*.xlsx'));
% Extract and display the names of XLSX files
xlsxFileNames = {files.name};
clear files; 

%% extract PSWE features table per patient:
feature_table_filepath_partial = '\‏‏PSWE_TUH.xlsx';
feature_table_filepath = strcat(main_path,feature_table_filepath_partial);

feature_table_per_patient=readtable(feature_table_filepath); 


%% Loop for each patient:
% Processing and feature extraction in batches

% Start EEGLAB
eeglab nogui


for i = 1:height(outpatientData) % loop for each patient
    % 1)Import and preprocess
    eegDataset=preprocessing_partial(outpatientData(i, :),main_path);
    Fs=eegDataset.hdr.srate;

    % 2) finding the ind in the xlsx file name:
    ind_xlsxfilenames=find_ind1(eegDataset,xlsxFileNames);

    % 3) extract only PSWE signal:
    if length(ind_xlsxfilenames)>0
    filenameWithExtension = strcat(global_pswe_per_patient_filepath,'\',xlsxFileNames(ind_xlsxfilenames));
    filenameWithExtension = strjoin(filenameWithExtension, ', ');
    Global_PSWE_properties=readtable(filenameWithExtension);
    % Global_PSWE_properties is a tables that contain start and end time for
    % each global PSWE and which electrode involved in the specific event

    % loop for each event:
    for t=1:length(Global_PSWE_properties.Start_ind) % t indicate the row in Global_PSWE_properties
% means that t= specific event

        %Global_PSWE_properties contain the start ind and end ind of PSWE
        %in sec
        start_index=Global_PSWE_properties.Start_ind(t)*Fs; % sec * Fs
        end_index=(Global_PSWE_properties.End_ind(t)+1)*Fs-1; % sec * Fs
        %extract the data from the specific electrode that involve in PSWE:
        for s=1:Global_PSWE_properties.Num_of_channels(t)
            electrode_num=table2array(Global_PSWE_properties(t,s+3));
            PSWE_data(s,:)=eegDataset.data(electrode_num,start_index:end_index);
            %PSWE_data stored only the PSWE spesific signals in the involved electrode.
        end
        
       % 4) for each signal: calculate the features
       x(t,:)=PSWE_features(PSWE_data,t,Fs,outpatientData.PatientId(i),main_path);
       clear PSWE_data;
    end

    % 5) for each patient:
    new_features_per_patient(i,1:30)=varfun(@mean, x); %mean on all electrodes
    new_features_per_patient.Patient(i)=outpatientData.PatientId(i);
    end
         
end

new_features_path_partial = '\new_features.xlsx';
new_features_path = strcat(main_path,new_features_path_partial);

writetable(new_features_per_patient,new_features_path);

%{
feature_table_filepath_partial = '\‏‏PSWE_TUH_characteristic.xlsx';
feature_table_filepath = strcat(main_path,feature_table_filepath_partial);
feature_table_per_patient=readtable(feature_table_filepath); 

new_features_path_partial = '\new_features.xlsx';
new_features_per_patient = strcat(main_path,new_features_path_partial);
new_features_per_patient=readtable(new_features_per_patient_path); 

joinedTable = outerjoin(feature_table_per_patient,new_features_per_patient, 'Keys', 'Patient', 'Type', 'left');

rows_to_modify = joinedTable.Persistent_slowing == 2; % Find rows where 'Persistent_slowing' column equals 2
joinedTable{rows_to_modify, 8:end} = 0; % Change columns 8 to the end to 0 for the selected rows

joinedTable_path_partial = '\joinedTable.xlsx';
joinedTable_path = strcat(main_path,joinedTable_path_partial);

writetable(joinedTable,joinedTable_path);
%}



end