function [Train_data,Val_data,Test_data]=separate_data(precent_train,precent_validation,Main_Table)
% this function seprate the data into train, validation and test.

epillepsy=find(Main_Table.Epilepsy_verified==1);
nonepillepsy=find(Main_Table.Epilepsy_verified==0);

%train:
numtrain_epi=round(precent_train*length(epillepsy)); 
numtrain_nonepi=round(precent_train*length(nonepillepsy));
%validation:
numval_epi=round(precent_validation*length(epillepsy));
numval_nonepi=round(precent_validation*length(nonepillepsy)); 
% test:
numtest_epi=length(epillepsy)-numtrain_epi-numval_epi;
numtest_nonepi=length(nonepillepsy)-numtrain_nonepi-numval_nonepi;

random_IDs = randperm(length(epillepsy)); 
selected_indices_train=random_IDs(1:numtrain_epi); selected_indices_val=random_IDs(numtrain_epi+1:numtrain_epi+numval_epi); selected_indices_test=random_IDs(numtrain_epi+numval_epi+1:end);
Train_epi = epillepsy(selected_indices_train); Val_epi = epillepsy(selected_indices_val); Test_epi = epillepsy(selected_indices_test);

random_IDs = randperm(length(nonepillepsy)); 
selected_indices_train=random_IDs(1:numtrain_nonepi); selected_indices_val=random_IDs(numtrain_nonepi+1:numtrain_nonepi+numval_nonepi); selected_indices_test=random_IDs(numtrain_nonepi+numval_nonepi+1:end);
Train_nonepi = nonepillepsy(selected_indices_train); Val_nonepi = nonepillepsy(selected_indices_val); Test_nonepi = nonepillepsy(selected_indices_test);

Train_data=Main_Table(Train_epi,:);
Train_data(numtrain_epi+1:numtrain_epi+numtrain_nonepi,:)=Main_Table(Train_nonepi,:);
Train_data = sortrows(Train_data, 'Patient');

Val_data=Main_Table(Val_epi,:);
Val_data(numval_epi+1:numval_epi+numval_nonepi,:)=Main_Table(Val_nonepi,:);
Val_data = sortrows(Val_data, 'Patient');

Test_data=Main_Table(Test_epi,:);
Test_data(numtest_epi+1:numtest_epi+numtest_nonepi,:)=Main_Table(Test_nonepi,:);
Test_data = sortrows(Test_data, 'Patient');

end