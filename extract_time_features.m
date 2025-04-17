function time_features = extract_time_features(signal)
    % signal is a cell array, where each cell contains a matrix
    % each row of the matrix represents data from one epoch
    % the returned time_features will also be a cell array
    
    num_patients = length(signal);
    time_features = cell(num_patients, 1);
    
    for patient = 1:num_patients
        patient_data = signal{patient};
        num_epochs = size(patient_data, 1);
        
        % Initialize feature matrix for current patient's epochs
        patient_features = zeros(num_epochs, 8);
        
        % Extract features for each epoch
        for epoch = 1:num_epochs
            epoch_signal = patient_data(epoch, :);
            
            % Calculate time domain features
            mean_val = mean(epoch_signal);
            var_val = var(epoch_signal, 0);
            skew_val = skewness(epoch_signal, 0);
            kurt_val = kurtosis(epoch_signal, 0);
            
            % Calculate zero-crossing rate (ZCR)
            zcr = sum(diff(sign(epoch_signal)) ~= 0);
            
            % Hjorth parameters:
            activity = var_val;
            mobility = std(diff(epoch_signal)) / std(epoch_signal);
            
            % Calculate complexity
            complexity = std(diff(diff(epoch_signal))) / std(diff(epoch_signal)) / mobility;
            
            % Combine all features into one row
            patient_features(epoch, :) = [mean_val, var_val, skew_val, kurt_val, zcr, activity, mobility, complexity];
        end
        
        % Store features for current patient
        time_features{patient} = patient_features;
    end
end
% function time_features = extract_time_features(signal)
%     % Each row of 'signal' is an individual signal.
%     
%     % Compute the mean for each signal (along columns)
%     mean_val = mean(signal, 2);
%     
%     % Compute the variance for each signal (along columns)
%     var_val = var(signal, 0, 2);
%     
%     % Compute the skewness for each signal (along columns)
%     skew_val = skewness(signal, 0, 2);
%     
%     % Compute the kurtosis for each signal (along columns)
%     kurt_val = kurtosis(signal, 0, 2);
%     
%     % Compute the zero-crossing rate (ZCR) for each signal:
%     % diff computes differences along the columns, and then sum counts the number of zero crossings.
%     zcr = sum(diff(sign(signal), 1, 2) ~= 0, 2);
%     
%     % Hjorth Parameters:
%     % Activity is simply the variance.
%     activity = var_val;
%     
%     % Mobility: ratio of the standard deviation of the first derivative to the standard deviation of the signal.
%     mobility = std(diff(signal, 1, 2), 0, 2) ./ std(signal, 0, 2);
%     
%     % Complexity: ratio of the standard deviation of the second derivative to that of the first derivative,
%     % normalized by the mobility.
%     complexity = std(diff(diff(signal, 1, 2), 1, 2), 0, 2) ./ std(diff(signal, 1, 2), 0, 2) ./ mobility;
%     
%     % Combine all features into one matrix:
%     % Each row corresponds to a signal, and each column is a computed feature.
%     time_features = [mean_val, var_val, skew_val, kurt_val, zcr, activity, mobility, complexity];
% end
