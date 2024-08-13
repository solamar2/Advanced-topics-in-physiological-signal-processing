function PSWE_Main_for_PSP
% Main for Final project on "Advanced topics in physiological signal processing"

% This function read the ‏‏'PSWE_TUH' xlsx that contain
% the PSWE features per patients and patients characterictics such as
% epilepsy\non-epilepsy , age, sex etc.

% Author: Sol Amara
% 25.01.2024

%% Read me
% need to change path : 
% main_path - the main foled with all the codes and data
% c in line 53

clear; clc; close all;
%% A) Extract information from TUH excel sheet
main_path='C:\Users\solam\OneDrive\Desktop\Master\Semester A\Advanced topics in physiological signal processing\final project';

Main_Table_path_partial = '\DATABASE\PSWE_TUH';
Main_Table_path = strcat(main_path,Main_Table_path_partial);

Main_Table=readtable(Main_Table_path);

%% Distribution of the data
%{
epillepsy=find(Main_Table.Epilepsy_verified==1);
nonepillepsy=find(Main_Table.Epilepsy_verified==0);

epi=find(Main_Table.Sex(epillepsy)==1);
epi=mean(Main_Table.Age(nonepillepsy));
%}


%% Dividing into train, validation and test tables:
% 80% for tarin and 20% for test
% 20% out of the train is for validation
% meaning: 0.64% train, %0.16% validation, 20% test

[Train_data,Val_data,Test_data]=separate_data(0.64,0.16,Main_Table);
%{
% check ditributions:
epillepsy=find(Train_data.Epilepsy_verified==1);
nonepillepsy=find(Train_data.Epilepsy_verified==0);

epi=find(Train_data.Sex(epillepsy)==1);
epi=std(Train_data.Age(epillepsy));

clearvars -except Train_data Test_data Val_data
%}

%% Feature selection:
c=8; % First col of the features
xTrain=table2array(Train_data(:,c:end));
yTrain=Train_data.Epilepsy_verified;
xVal=table2array(Val_data(:,c:end));
yVal=Val_data.Epilepsy_verified;

training_data=[];
other_features=xTrain;

for m=1:10 % m = number of features
   models_data{m}=find_optimal_feature(training_data,other_features,xTrain,yTrain,xVal,yVal);
    training_data=models_data{1, m}.training_data;
    other_features=models_data{1, m}.other_features;
    J(m)=models_data{1,m}.auc;
end

figure
scatter(1:10,J, 'filled');
ylabel('J - optimal creiterion')
xlabel(' m - number of features')

m_optimal=find(J==max(J)); m_optimal=m_optimal(1)
training_data_optimal=models_data{1,m_optimal}.training_data;
[~, cols_selected_features] = ismember(training_data_optimal', xTrain', 'rows'); 

mdl = models_data{1,m_optimal}.model;

[~, scores] = predict(mdl, training_data_optimal);
[X_ROC, Y_ROC, ~, AUC] = perfcurve(yTrain, scores(:, 2), 1);

% Plot ROC curve
figure;
plot(X_ROC, Y_ROC);
xlabel('False Positive Rate');
ylabel('True Positive Rate');
title('ROC Curve');
grid on;

%% checking performance on the Test data
clearvars -except mdl Test_data cols_selected_features c

%load('Test_data.mat')
xTest=table2array(Test_data(:,c-1+cols_selected_features));
yTest=Test_data.Epilepsy_verified;

[~, scores] = predict(mdl, xTest); % Predicted probabilities
[X_ROC, Y_ROC, T, AUC] = perfcurve(yTest, scores(:, 2), 1);

% Calculate sensitivity and specificity for each threshold
sensitivity = Y_ROC;
specificity = 1 - X_ROC;
% Calculate Youden's J statistic
J = sensitivity + specificity - 1;
% Find the index of the optimal point (maximal J)
[optimal_J, optimal_idx] = max(J);
% Display the optimal point
optimal_threshold = T(optimal_idx);
optimal_sensitivity = sensitivity(optimal_idx);
optimal_specificity = specificity(optimal_idx);


y_pred=scores(:, 2)>= optimal_threshold;
% Accuracy
accuracy = sum(y_pred == yTest) / numel(yTest);

indepi=find(yTest==1);sum(y_pred(indepi))
indnonepi=find(yTest==0);31-sum(y_pred(indnonepi))

%plot results:
figure;
plot(X_ROC, Y_ROC);
hold on
scatter(1-optimal_specificity,optimal_sensitivity,'filled');
xlabel('False Positive Rate');
ylabel('True Positive Rate');
title('ROC Curve');
grid on;

%% ODD Ratio : 
coefficients = mdl.Coefficients.Estimate;
feature_index = find(abs(coefficients)==max(abs(coefficients))); 
feature_coefficient = coefficients(feature_index)


%% Replicate the article model:
clear; clc; close all;
main_path='C:\Users\solam\OneDrive\Desktop\Master\Semester A\Advanced topics in physiological signal processing\final project';

Main_Table_path_partial ='\DATABASE\PSWE_TUH';
Main_Table_path = strcat(main_path,Main_Table_path_partial);

Main_Table=readtable(Main_Table_path);

[Train_data,~,Test_data]=separate_data(0.8,0,Main_Table);

% in this reaserch they used only the 4 features:
xTrain=table2array(Train_data(:,8:11));
yTrain=Train_data.Epilepsy_verified;
xTest=table2array(Test_data(:,8:11));
yTest=Test_data.Epilepsy_verified;

% training the model:
mdl = fitglm(xTrain, yTrain, 'Distribution', 'binomial', 'Link', 'logit');

[~,scores] = predict(mdl, xTest); % Predicted probabilities
[X_ROC, Y_ROC, T, AUC] = perfcurve(yTest, scores(:, 2), 1);


% Calculate sensitivity and specificity for each threshold
sensitivity = Y_ROC;
specificity = 1 - X_ROC;
% Calculate Youden's J statistic
J = sensitivity + specificity - 1;
% Find the index of the optimal point (maximal J)
[optimal_J, optimal_idx] = max(J);
% Display the optimal point
optimal_threshold = T(optimal_idx);
optimal_sensitivity = sensitivity(optimal_idx);
optimal_specificity = specificity(optimal_idx);


y_pred= scores(:, 2)>= optimal_threshold;
% Accuracy
accuracy = sum(y_pred == yTest) / numel(yTest);

indepi=find(yTest==1);sum(y_pred(indepi))
indnonepi=find(yTest==0);31-sum(y_pred(indnonepi))

%plot results:
figure;
plot(X_ROC, Y_ROC);
hold on
scatter(1-optimal_specificity,optimal_sensitivity,'filled');
xlabel('False Positive Rate');
ylabel('True Positive Rate');
title('ROC Curve');
grid on;

% ODD ration: 
coefficients = mdl.Coefficients.Estimate;
feature_index =  find(abs(coefficients)==max(abs(coefficients))); 
feature_coefficient = coefficients(feature_index);



