function best_feature=find_optimal_feature(training_data,other_features,xTrain,yTrain,xVal,yVal)
% this function is a part of the FSS algoritm, 
% training_data is the already chosen features, and this function check who
% is the next features that will chosen

d= size(other_features,2); % intial number of features
auc=zeros(1,d);
models = cell(1, d); % num_iterations is the number of iterations

 for i=1:d
     check=[training_data, other_features(:,i)];
     [~,col]= ismember( check',xTrain','rows');
     current_val_data=xVal(:,col);
     [auc(i),models{i}] = AUC_of_Logistic_regression(check, yTrain,current_val_data,yVal);
 end

 ind=find(auc==max(auc));
 training_data=[training_data,other_features(:,ind)];
 other_features(:,ind)=[];

 best_feature.training_data=training_data;
 best_feature.other_features=other_features;
 best_feature.auc=auc(ind);
 best_feature.model=models{ind};

