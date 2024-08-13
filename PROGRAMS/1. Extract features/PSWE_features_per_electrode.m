function  PSWE_features_per_electrode(PSWE_data,Fs,filenameWithExtension)
% calculate features per electrode in 2 sec segment with 1 sec overlap

% Input: only the signal data in the part of spesific PSWE. each row is
% diffrent electrode that detcted this event.

% OUTPUT: xlsx file in 'features per event' folder for every single event
% seperatly. sheet represent spsific electrode results.

% Author: Sol Amara
% 25.01.2024

%% LOOP FOR EVERY ELECTRODE SEPERATLY IN ORDER TO CALCULATE THE FEATURES:
features_names = {'Energy', 'Maximalampllitude', 'Amplitudeintegral', 'STD', 'Skewness', 'Kutrosis', 'Curvelength', 'Harmonicity', 'Zerocrossing', 'fractledimnsion', 'Peak', 'meanNE', 'STDNE', 'KutrosisNE', 'SkewnessNE', 'STDFFT', 'KutrosisFFT', 'SkewnessFFT', 'meadianFFT', 'SumFFT', 'delta', 'theta', 'alpha', 'beta', 'lowgamma', 'highgamma','MPFperseg'};

lwindow = Fs*2;
num_of_electrode=size(PSWE_data,1);
for i=1:num_of_electrode %for every electrode:
    bdata = buffer (PSWE_data(i,:),lwindow,lwindow*0.5);%buffering the data
  
    %Nvision 26 features:
    x= Nvision_parameters_5_vs_26 (bdata,Fs);
    k=size(x,2)+1;
    
    
    % Initialize medFreq array
    numSegments = size(bdata, 2);medFreq = zeros(1, numSegments);

    % Compute median power frequency for each buffer segment
    for j = 1:numSegments
        % Load buffer segment
        y = bdata(:, j);
    
        % Compute median power frequency using built-in function
        [medFreq(j), ~] = medfreq(y, Fs);
    end

    x(:,k)=medFreq'; 
    
     %% Save each window features to xlsx:
    x = array2table(x, 'VariableNames', features_names);
    writetable(x, filenameWithExtension,'Sheet', i);
    % sheet number is for each electrode.
    clear bdata;
end




























end