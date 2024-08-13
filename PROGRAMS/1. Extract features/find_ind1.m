function ind_xlsxfilenames=find_ind1(eegDataset,xlsxFileNames)
% this function compare the patient ID, session and time and find the
% specific row in feature table and xlsxfilenames that incidcate the same
% patient as in eegdataset.
% used in : 'PSWE_Main_for_features_calculation' function. 

 %% 2.2: xlsx file name:
    [~, filename, ~] = fileparts(eegDataset.edf_file_names);
    filename=strcat(filename,'.edf.xlsx');
    ind_xlsxfilenames=find(strcmp(xlsxFileNames,filename));

end