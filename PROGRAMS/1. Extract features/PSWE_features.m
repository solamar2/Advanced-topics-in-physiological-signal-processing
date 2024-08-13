function  features=PSWE_features(PSWE_data,t,Fs,patientID,main_path)
% After detecting global PSWE, this function calculate the characteristics
% for each PSWE in each electrode that involve in the event.
% PSWE start, end and involved electroed stored in xlsx files in folder
% called 'per patient'
% The characteristics for each patient stored in 'TUH_PSWE_characteristic'

% Input: 
% Pswe_data: only the signal data in the part of spesific PSWE (t). each row is diffrent electrode that detcted this event.
% Global_PSWE_properties : table for each patient with the time of start,
% end and elctrode that involved
% t: the spesific event in this patient


%%
filename = '\Results\features per event\';

filenameWithExtension = strcat(main_path,filename,num2str(patientID),'_event_number_',num2str(t),'_per_segment', '.xlsx');

PSWE_features_per_electrode(PSWE_data,Fs,filenameWithExtension);
num_of_electrode=size(PSWE_data,1);

for i=1:num_of_electrode
electrode_features_per_segment = readtable(filenameWithExtension, 'Sheet', i);
electrode_features(i,:)=varfun(@mean, electrode_features_per_segment); %mean on all the segments
std_MPF(i)=std(electrode_features_per_segment.MPFperseg);
meadian_MPF(i)=median(electrode_features_per_segment.MPFperseg);
high_MPF=max(electrode_features_per_segment.MPFperseg);
end

features=varfun(@mean, electrode_features); %mean on all electrodes
features.stdMPF=mean(std_MPF);
features.meadianMPF=mean(meadian_MPF);
features.highMPF=mean(high_MPF);

filenameWithExtension2 = strcat(main_path,filename,num2str(patientID),'_event_number_',num2str(t),'.xlsx');
writetable(features, filenameWithExtension2,'Sheet', 'mean on electrodes');
writetable(electrode_features, filenameWithExtension2,'Sheet', 'data per electrode');

end