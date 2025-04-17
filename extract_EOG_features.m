function EOG_features = extract_EOG_features(signal)
    % Extract EOG features from signal data
    % Input:
    %   signal - a cell array where each cell contains a matrix
    %            each row of the matrix represents data from one epoch
    % Output:
    %   EOG_features - cell array containing extracted features for each patient
    
    num_patients = length(signal);
    EOG_features = cell(num_patients, 1);
    
    % Parameters
    epochLength = 30; % in seconds
    Fs = 50; % sampling frequency
    
    % Define frequency bands for spectral analysis
    slowBand = [0.3 0.8]; % Slow eye movement band
    remBand = [0.8 5];    % REM band
    
    for patient = 1:num_patients
        patient_data = signal{patient};
        num_epochs = size(patient_data, 1);
        
        % Pre-allocate matrix: rows = epochs, columns = features
        % Order of features:
        % 1. blinkCount
        % 2. movementDensity
        % 3. blinkRate
        % 4. avgBlinkAmplitude
        % 5. avgBlinkWidth
        % 6. slowPower
        % 7. remPower
        % 8. ratio_RemToSlow
        features = zeros(num_epochs, 8);
        
        for epoch = 1:num_epochs
            % Get epoch signal
            epoch_signal = patient_data(epoch, :);
            
            % Preprocessing
            epoch_signal = detrend(epoch_signal); % Remove linear trend
            
            % Apply bandpass filter (0.1-20 Hz) to focus on relevant EOG frequencies
            [b, a] = butter(4, [0.1 20]/(Fs/2), 'bandpass');
            filtered_signal = filtfilt(b, a, epoch_signal);
            
            % Compute adaptive threshold for blink detection
            % Typically blinks have higher amplitude than baseline
            threshold = mean(filtered_signal) + 2.5 * std(filtered_signal);
            
            % 1-5: Blink detection and characteristics
            [blinkAmplitudes, blinkLocations, blinkWidths] = findpeaks(filtered_signal, ...
                'MinPeakHeight', threshold, ...
                'MinPeakDistance', round(0.2*Fs), ... % Minimum 200ms between blinks
                'WidthReference', 'halfheight');
            
            if ~isempty(blinkAmplitudes)
                % Convert widths from samples to seconds
                blinkWidths_sec = blinkWidths / Fs;
                
                % 1. Blink count
                features(epoch, 1) = length(blinkAmplitudes);
                
                % 3. Blink rate (blinks per minute)
                features(epoch, 3) = features(epoch, 1) / (epochLength / 60);
                
                % 4. Average blink amplitude
                features(epoch, 4) = mean(blinkAmplitudes);
                
                % 5. Average blink width (in seconds)
                features(epoch, 5) = mean(blinkWidths_sec);
            else
                % No blinks detected
                features(epoch, 1) = 0;
                features(epoch, 3) = 0;
                features(epoch, 4) = 0;
                features(epoch, 5) = 0;
            end
            
            % 2. Movement density (improved)
            % Use a smoothed derivative and threshold to reduce noise sensitivity
            diff_signal = diff(filtered_signal);
            smooth_diff = movmean(abs(diff_signal), round(0.05*Fs)); % 50ms window
            movement_threshold = 0.2 * std(smooth_diff);
            features(epoch, 2) = sum(smooth_diff > movement_threshold) / length(smooth_diff);
            
            % 6-8: Spectral analysis
            % Using Welch's method with 4-second windows and 50% overlap
            window_length = 4 * Fs;
            overlap = window_length / 2;
            [pxx, f] = pwelch(filtered_signal, window_length, overlap, [], Fs);
            
            % 6. Slow eye movement power
            features(epoch, 6) = bandpower(pxx, f, slowBand, 'psd');
            
            % 7. REM power
            features(epoch, 7) = bandpower(pxx, f, remBand, 'psd');
            
            % 8. Ratio of REM power to slow power
            if features(epoch, 6) > 0
                features(epoch, 8) = features(epoch, 7) / features(epoch, 6);
            else
                features(epoch, 8) = 0;
            end
        end
        
        EOG_features{patient} = features;
    end
    
    % Add feature names to the output for better readability
    feature_names = {'blinkCount', 'movementDensity', 'blinkRate', ...
                    'avgBlinkAmplitude', 'avgBlinkWidth', 'slowPower', ...
                    'remPower', 'ratio_RemToSlow'};
    
%     % Create a structure for each patient with named fields
%     named_features = cell(num_patients, 1);
%     for patient = 1:num_patients
%         patient_struct = struct();
%         for feat = 1:length(feature_names)
%             patient_struct.(feature_names{feat}) = EOG_features{patient}(:, feat);
%         end
%         named_features{patient} = patient_struct;
%     end
%     
%     EOG_features = named_features;
end