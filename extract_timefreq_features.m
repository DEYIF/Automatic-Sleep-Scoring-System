function features_timefreq = extract_timefreq_features(EEG, Fs)
    % EEG is a cell array where each cell contains a matrix
    % each row of the matrix represents data from one epoch
    % Fs is the sampling frequency in Hz
    
    num_patients = length(EEG);
    features_timefreq = cell(num_patients, 1);
    
    for patient = 1:num_patients
        patient_data = EEG{patient};
        num_epochs = size(patient_data, 1);
        num_samples = size(patient_data, 2);
        
        % Initialize feature matrix for current patient
        patient_features = zeros(num_epochs, 2); % One column for Wavelet energy, another for Spectral Entropy
        
        for epoch = 1:num_epochs
            signal = patient_data(epoch, :);
            
            % --- Wavelet Coefficients ---
            [cA, cD] = dwt(signal, 'db4'); % Decomposition with Daubechies 4 wavelet
            wavelet_energy = sum(abs(cD).^2); % Energy of detailed coefficients
            
            % --- Spectral Entropy ---
            psd = abs(fft(signal)).^2; % Compute Power Spectral Density (PSD)
            psd = psd(1:floor(num_samples/2)); % Use only the positive half
            psd = psd / sum(psd); % Normalize to get probabilities
            spectral_entropy = -sum(psd .* log2(psd + eps)); % Compute spectral entropy
            
            % Store features
            patient_features(epoch, :) = [wavelet_energy, spectral_entropy];
        end
        
        % Store features for current patient
        features_timefreq{patient} = patient_features;
    end
end
% function features_timefreq = extract_timefreq_features(EEG, Fs)
%     % Number of epochs and samples per epoch
%     [num_epochs, num_samples] = size(EEG);
%     
%     % Initialize feature matrix
%     features_timefreq = zeros(num_epochs, 2); % One column for Wavelet energy, another for Spectral Entropy
%     
%     for i = 1:num_epochs
%         signal = EEG(i, :);
%         
%         % --- Wavelet Coefficients ---
%         [cA, cD] = dwt(signal, 'db4');  % Decomposition with Daubechies 4 wavelet
%         wavelet_energy = sum(abs(cD).^2); % Energy of detailed coefficients
%         
%         % --- Spectral Entropy ---
%         psd = abs(fft(signal)).^2;  % Compute Power Spectral Density (PSD)
%         psd = psd(1:floor(num_samples/2)); % Use only the positive half
%         psd = psd / sum(psd); % Normalize to get probabilities
%         spectral_entropy = -sum(psd .* log2(psd + eps)); % Compute spectral entropy
%         
%         % Store features
%         features_timefreq(i, :) = [wavelet_energy, spectral_entropy];
%     end
% end
