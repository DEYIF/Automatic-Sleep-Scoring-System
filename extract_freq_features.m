function [freq_features, allPxx, freqs] = extract_freq_features(signal, Fs)
    % Each row represents an individual signal.
    [numSignals, N] = size(signal);
    
    % Create a common frequency axis (considering only positive frequencies)
    freqs = linspace(0, Fs/2, round(N/2+1)); 

    % Preallocate output:
    % Total features: 5 band powers + 1 spectral edge frequency + 5 relative band powers = 11 features
    freq_features = zeros(numSignals, 11);
    %allPxx = cell(numSignals, 1); % Store power spectral density for each signal
    allPxx = zeros(numSignals,round(N/2+1));

    % Define EEG frequency bands: Delta, Theta, Alpha, Beta, Gamma
    bands = [0.5 4; 4 8; 8 12; 12 30; 30 40];

    for s = 1:numSignals
        sig = signal(s, :);
        fft_signal = fft(sig);
        % Compute the normalized power spectral density (PSD) and take only the positive frequencies
        Pxx = abs(fft_signal(1:round(N/2+1))).^2 / N;
        allPxx(s,:) = Pxx;
        
        % Compute power in each frequency band
        power_bands = zeros(1, size(bands,1));
        for i = 1:size(bands,1)
            idx = find(freqs >= bands(i,1) & freqs < bands(i,2));
            if ~isempty(idx)
                idx = round(idx);  % Ensure indices are integers
                idx = max(1, min(idx, length(Pxx))); % Keep indices within valid range
                power_bands(i) = trapz(freqs(idx), Pxx(idx)); 
            else
                power_bands(i) = NaN;
            end
        end

        % Compute Spectral Edge Frequency (SEF at 80%)
        total_power = sum(Pxx);
        cum_power = cumsum(Pxx) / total_power;
        spec_edge_idx = find(cum_power >= 0.8, 1, 'first');
        if isempty(spec_edge_idx) || spec_edge_idx > length(freqs)
            spectral_edge = NaN;
        else
            spec_edge_idx = int32(round(spec_edge_idx)); % Ensure it is an integer
            spec_edge_idx = max(1, min(spec_edge_idx, length(freqs)));
            spectral_edge = freqs(spec_edge_idx);
        end

        % Compute relative power for each frequency band
        total_band_power = sum(power_bands, 'omitnan');
        relative_power = power_bands / total_band_power;
        
        % Combine feature vector: [band powers, spectral edge frequency, relative band powers]
        freq_features(s, :) = [power_bands, spectral_edge, relative_power];
    end
end
