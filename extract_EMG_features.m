function features_cell = extract_EMG_features(signal_cell, env_cell, Fs)
% extract_EMG_features  对多名受试者、多个 epoch 的 EMG 信号提取特征
%
% 输入：
%   signal_cell - cell 数组，{p} 是第 p 位受试者的 raw 信号矩阵 [epochs × samples]
%   env_cell    - cell 数组，{p} 是第 p 位受试者的包络矩阵     [epochs × samples]
%   Fs          - 采样频率（Hz）
%
% 输出：
%   features_cell - cell 数组，{p} 是 struct，字段为：
%       .RMS             [epochs×1] – 均方根幅值
%       .MeanFrequency   [epochs×1] – 平均频率
%       .BurstRate       [epochs×1] – 每分钟爆发次数
%       .MuscleTone      [epochs×1] – 平均肌张力（包络均值）

    num_patients = numel(signal_cell);
    features_cell = cell(num_patients, 1);
    

    for p = 1:num_patients
        % 强制双精度，防止 pwelch 等函数报错
        sig_mat = double(signal_cell{p});  
        env_mat = double(env_cell{p});
        [num_epochs, num_samples] = size(sig_mat);

        features = zeros(num_epochs, 4);
        % 预分配
        RMS   = zeros(num_epochs,1);
        MF    = zeros(num_epochs,1);
        BR    = zeros(num_epochs,1);
        MT    = zeros(num_epochs,1);

        for e = 1:num_epochs
            x   = sig_mat(e, :);
            env = env_mat(e, :);

            % ——1. RMS Amp——
            features(e, 1) = rms(x);

            % ——2. Mean Frequency——
            [pxx, f] = pwelch(x, [], [], [], Fs);
            features(e, 2) = sum(f .* pxx) / sum(pxx);

            % ——3. Burst Rate——
            thr = 0.4 * max(env);
            mask = env > thr;
            dm = diff([0, mask, 0]);
            starts = find(dm == 1);
            ends   = find(dm == -1) - 1;
            durations = (ends - starts) / Fs;
            valid = durations > 0.05;           % 过滤 <50 ms 的伪爆发
            n_bursts = sum(valid);
            epoch_min = num_samples / (Fs * 60);
            features(e, 3) = n_bursts / epoch_min;

            % ——4. Muscle Tone——
            MT(e) = mean(env);
            features(e, 4) = mean(env);

        end

        features_cell{p} = features;
    end
end

% function EMG_features=extract_EMG_features(signal, EMG_envelope, Fs)   
% %RMS amplitude
% %It captures muscle activity intensity.
% % Wakefulness → High RMS (high muscle tone and movement).
% % NREM sleep → Lower RMS (progressively decreasing through N1, N2, and N3).
% % REM sleep → Very low RMS (due to muscle atonia, meaning the muscles are paralyzed to prevent acting out dreams).
%     rms_value = rms(signal);
% 
% 
% %Frequency content
% %Checks how fast the muscle activity is oscillating.
% % Wakefulness and light sleep (N1) → high-frequency components (due to voluntary movements, small twitches).
% % N2, N3 and REM → low-frequency components.
%     [pxx, f] = pwelch(signal, [], [], [], Fs);
%     mean_frequency = sum(f .* pxx) / sum(pxx);
% 
% 
% %Burst detection 
% %Identifies episodes of short, voluntary muscle contractions.
% % Wakefulness → frequent bursts (due to muscle adjustments, posture shifts).
% % NREM sleep → bursts decrease (few movements).
% % REM sleep → bursts should be almost absent due to muscle atonia.
%     burst_threshold = 0.4 * max(EMG_envelope); % cutoff to detect real bursts without detecting every small fluctuation
%     burst_mask = EMG_envelope > burst_threshold;
%     burst_mask = burst_mask(:);
%     burst_diff = diff([0; burst_mask; 0]);
%     burst_starts = find(burst_diff == 1);
%     burst_ends = find(burst_diff == -1) - 1;
% 
%     min_burst_duration = 0.05; % Only keep bursts that are longer than 0.05 seconds (50 ms), because very short spikes are probably noise
%     valid_bursts = (burst_ends - burst_starts) / Fs > min_burst_duration;
%     n_bursts = sum(valid_bursts);
% 
%     epoch_duration_min = length(signal) / (Fs * 60); % Duration in minutes
%     burst_rate = n_bursts / epoch_duration_min;       % Bursts per minute
% 
% 
% %Muscle tone analysis
% %Muscle tone is the background activation of the muscle when the person is not moving.
% % Wakefulness → High muscle tone.
% % N1/N2 (light sleep) → Tone decreases.
% % N3 (deep sleep) → Very low tone.
% % REM sleep → Muscle tone disappears (almost zero).
%     muscle_tone = mean(EMG_envelope);
% 
% 
% %Output 
%     EMG_features.RMS = rms_value;
%     EMG_features.MeanFrequency = mean_frequency;
%     EMG_features.BurstRate = burst_rate;
%     EMG_features.MuscleTone = muscle_tone;
% 
% 
% end
% 
