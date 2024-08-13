function [auc,mdl]= AUC_of_Logistic_regression(xTrain, yTrain,xVal,yVal)
% This function get vector of features and vector of predict value (1=
% epilepsy anf 0 non-epilepsy) and training a logitic regression model
% After training, the output of function is AUC of this model 

    % Train logistic regression model
    mdl = fitglm(xTrain, yTrain, 'Distribution', 'binomial', 'Link', 'logit');

    % Predict probabilities on the validation data
    [~, scores] = predict(mdl, xVal);

    % Compute AUC on the training data
    [~, ~, ~, auc] = perfcurve(yVal, scores(:, 2), 1);
end




