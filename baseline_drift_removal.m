function filtered = baseline_drift_removal(origin, Fs, Fc)
    if nargin < 3
        Fc = 0.5;  % default Fc
    end
    [b,a] = butter(4,Fc/(Fs/2),'high'); % 4-order butter hp filter
    filtered = filtfilt(b,a,origin);
    t = linspace(0,length(origin)/Fs,length(origin));

    figure;
    subplot(3,1,1)
    plot(t, origin)
    xlim([0,t(end)])
    title('Original Signal')
    
    subplot(3,1,2)
    plot(t, filtered)
    xlim([0,t(end)])
    title('Filtered Signal')
    
    subplot(3,1,3)
    plot(t, origin, t, filtered)
    xlim([0,t(end)])
    legend('original', 'filtered')
    title('Comparison')
end