function time_features = extract_time_features(signal)
    % Each row of 'signal' is an individual signal.
    
    % Compute the mean for each signal (along columns)
    mean_val = mean(signal, 2);
    
    % Compute the variance for each signal (along columns)
    var_val = var(signal, 0, 2);
    
    % Compute the skewness for each signal (along columns)
    skew_val = skewness(signal, 0, 2);
    
    % Compute the kurtosis for each signal (along columns)
    kurt_val = kurtosis(signal, 0, 2);
    
    % Compute the zero-crossing rate (ZCR) for each signal:
    % diff computes differences along the columns, and then sum counts the number of zero crossings.
    zcr = sum(diff(sign(signal), 1, 2) ~= 0, 2);
    
    % Hjorth Parameters:
    % Activity is simply the variance.
    activity = var_val;
    
    % Mobility: ratio of the standard deviation of the first derivative to the standard deviation of the signal.
    mobility = std(diff(signal, 1, 2), 0, 2) ./ std(signal, 0, 2);
    
    % Complexity: ratio of the standard deviation of the second derivative to that of the first derivative,
    % normalized by the mobility.
    complexity = std(diff(diff(signal, 1, 2), 1, 2), 0, 2) ./ std(diff(signal, 1, 2), 0, 2) ./ mobility;
    
    % Combine all features into one matrix:
    % Each row corresponds to a signal, and each column is a computed feature.
    time_features = [mean_val, var_val, skew_val, kurt_val, zcr, activity, mobility, complexity];
end
