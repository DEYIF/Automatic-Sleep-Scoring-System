function EOG_features=extract_EOG_features(signal)    

    % Blink Detection 
    %In an EOG blinks appear as large amplitude spikes, so they will be
    %detected using a threshold 
    threshold=80; %choosen based on the EOG plot
    blinkDetection= sum(max(signal, [], 2) > threshold);  

    % Movement density 
    %Total amount of activity (movement or change) in the signal
    movementDensity= sum(abs(diff(signal, 1, 2)), 2);

    % Blink rate and characteristics 
    [blinkAmplitude, blinkLocation, blinkWidth]= findpeaks(signal(1,:), 'MinPeakHeight', threshold); 
    %findpeaks() detects peaks in a single epoch
    %Rate and characteristics calculated:
    epochLength=30;
    Fs=200;
    blinkRate= length(blinkLocation) / (epochLength / 60); %number of blinks per minute
    avgBlinkAmplitude= mean(blinkAmplitude); %average blink amplitude in microvolts
    avgBlinkWidth= mean (blinkWidth/Fs); %duration in seconds of a blink

    % Sleep Eye movement detection ->  Slow/rapid eye movement patterns
    %We know that rapid movements correspond with high frequencies and slow
    %movements correspond with low frequencies. Slow eye movements have a
    %band frequency of 0.3-0.8 Hz, meanwhile, REM have a band frequency of
    %0.8-5Hz. 
    % Use Welch's method to compute Power Spectral Density (PSD) for an epoch
    [pxx_L,f]= pwelch(signal(1,:), [], [], [], Fs);
    % Calculate power in the slow movement band (0.3–0.8 Hz)
    slowPower= bandpower(pxx_L, f, [0.3 0.8], 'psd');
    % Calculate power in the REM band (0.8–5 Hz)
    remPower= bandpower(pxx_L, f, [0.8 5], 'psd');
    %Ratio used to understand the dominant type of mevement:
    ratio_RemToSlow= remPower/slowPower;
    
    
    %gather all the extracted features and assign it to EOG_features
    EOG_features.blinkDetection      = blinkDetection;
    EOG_features.movementDensity     = movementDensity;
    EOG_features.blinkRate           = blinkRate;
    EOG_features.avgBlinkAmplitude   = avgBlinkAmplitude;
    EOG_features.avgBlinkWidth       = avgBlinkWidth;
    EOG_features.slowPower           = slowPower;
    EOG_features.remPower            = remPower;
    EOG_features.ratio_RemToSlow     = ratio_RemToSlow;


end 