function [trainedClassifier, validationAccuracy, f1Score, perStageAccuracy] = trainClassifier3(trainingData, responseData)

% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19'};
predictors = inputTable(:, predictorNames);
response = responseData;
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false];

% Train a classifier
% This code specifies all the classifier options and trains the classifier.
classificationKNN = fitcknn(...
    predictors, ...
    response, ...
    'Distance', 'Euclidean', ...
    'Exponent', [], ...
    'NumNeighbors', 1, ...
    'DistanceWeight', 'Equal', ...
    'Standardize', true, ...
    'ClassNames', categorical({'REM'; 'None'; 'N3'; 'N2'; 'N1'; 'Wake'}, {'REM' 'None' 'N3' 'N2' 'N1' 'Wake'}));

% Create the result struct with predict function
predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
knnPredictFcn = @(x) predict(classificationKNN, x);
trainedClassifier.predictFcn = @(x) knnPredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
trainedClassifier.ClassificationKNN = classificationKNN;
trainedClassifier.About = 'This struct is a trained model exported from Classification Learner R2022b.';
trainedClassifier.HowToPredict = sprintf('To make predictions on a new predictor column matrix, X, use: \n  yfit = c.predictFcn(X) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nX must contain exactly 19 columns because this model was trained using 19 predictors. \nX must contain only predictor columns in exactly the same order and format as your training \ndata. Do not include the response column or any columns you did not import into the app. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19'};
predictors = inputTable(:, predictorNames);
response = responseData;
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false];

validationPredictions = predict(trainedClassifier.ClassificationKNN, trainingData);
validationAccuracy = sum(validationPredictions == responseData)/numel(responseData);

% Now, calculate the F1 score

% Create confusion matrix
cm = confusionmat(responseData, validationPredictions);

% Calculate precision, recall, and F1 score for each class
precision = diag(cm) ./ sum(cm, 2); % Precision per class
recall = diag(cm) ./ sum(cm, 1)';   % Recall per class

% Avoid NaN values by setting them to zero (for classes with no predictions)
precision(isnan(precision)) = 0;
recall(isnan(recall)) = 0;

% Calculate F1 score per class
f1 = 2 * (precision .* recall) ./ (precision + recall);

% The F1 score for the entire model is the average of F1 scores per class
f1Score = mean(f1);

% Calculate per-stage accuracy (accuracy for each class)
perStageAccuracy = diag(cm) ./ sum(cm, 2); % Accuracy per class

% Replace NaN values (if any class has no instances in the validation set) with 0
perStageAccuracy(isnan(perStageAccuracy)) = 0;
end
