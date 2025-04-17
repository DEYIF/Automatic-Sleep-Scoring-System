function [EOG_L_filt, EOG_R_filt] = preprocess_eog(EOG_L, EOG_R, Fs)
    % EOG_L and EOG_R are cell arrays where each cell contains a matrix
    % each row of the matrix represents data from one epoch
    % Fs is the sampling frequency
    
    num_patients = length(EOG_L);
    EOG_L_filt = cell(num_patients, 1);
    EOG_R_filt = cell(num_patients, 1);
    
    % Design low-pass Butterworth filter (cutoff at 15 Hz)
    [b, a] = butter(4, 15 / (Fs / 2), 'low');
    
    for patient = 1:num_patients
        % Get current patient's data
        patient_EOG_L = EOG_L{patient};
        patient_EOG_R = EOG_R{patient};
        
        % Get dimensions
        num_epochs = size(patient_EOG_L, 1);
        
        % Initialize filtered data matrices for current patient
        patient_EOG_L_filt = zeros(size(patient_EOG_L));
        patient_EOG_R_filt = zeros(size(patient_EOG_R));
        
        % Process each epoch for this patient
        for epoch = 1:num_epochs
            % Extract current epoch
            epoch_EOG_L = patient_EOG_L(epoch, :);
            epoch_EOG_R = patient_EOG_R(epoch, :);
            
            % Remove baseline (mean subtraction)
            epoch_EOG_L = epoch_EOG_L - mean(epoch_EOG_L);
            epoch_EOG_R = epoch_EOG_R - mean(epoch_EOG_R);
            
            % Apply filter
            patient_EOG_L_filt(epoch, :) = filtfilt(b, a, epoch_EOG_L);
            patient_EOG_R_filt(epoch, :) = filtfilt(b, a, epoch_EOG_R);
        end
        
        % Store processed data for current patient
        EOG_L_filt{patient} = patient_EOG_L_filt;
        EOG_R_filt{patient} = patient_EOG_R_filt;
    end
end
% function [EOG_L_filt, EOG_R_filt] = preprocess_eog(EOG_L, EOG_R, Fs)
%     % Remove baseline (mean subtraction)
%     EOG_L = EOG_L - mean(EOG_L, 2);
%     EOG_R = EOG_R - mean(EOG_R, 2);
% 
%     % Low-pass Butterworth filter (cutoff at 15 Hz)
%     [b, a] = butter(4, 15 / (Fs / 2), 'low');
% 
%     % Filter each epoch (row-wise)
%     EOG_L_filt = filtfilt(b, a, EOG_L')';
%     EOG_R_filt = filtfilt(b, a, EOG_R')';
% end
% 