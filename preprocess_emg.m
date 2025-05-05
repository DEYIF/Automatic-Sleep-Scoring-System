function [emg_clean_cell, emg_envelope_cell] = preprocess_emg(emg_raw_cell, Fs)
% preprocess_emg  对多名受试者的 EMG 数据按 epoch 进行预处理
% 
% 输入：
%   emg_raw_cell   - cell 数组，每个 cell 为 一个 matrix，[num_epochs × num_samples]
%   Fs             - 采样频率
%
% 输出：
%   emg_clean_cell    - cell 数组，每个 cell 为 去伪迹、平滑后的包络 [num_epochs × num_samples]
%   emg_envelope_cell - cell 数组，每个 cell 为 原始包络         [num_epochs × num_samples]

    num_patients = numel(emg_raw_cell);
    emg_clean_cell    = cell(num_patients, 1);
    emg_envelope_cell = cell(num_patients, 1);

    % 设计带通滤波器：20–450 Hz（如果 Fs<=900，上限取 Fs/2-1）
    if Fs <= 900
        high_cutoff = Fs/2 - 1;
    else
        high_cutoff = 450;
    end
    [b, a] = butter(4, [20, high_cutoff] / (Fs/2), 'bandpass');

    for p = 1:num_patients
        % 强制转为 double，确保 filtfilt 能运行
        data_raw = double(emg_raw_cell{p});  
        [num_epochs, num_samples] = size(data_raw);

        env_mat   = zeros(num_epochs, num_samples);
        clean_mat = zeros(num_epochs, num_samples);

        for e = 1:num_epochs
            x = data_raw(e, :);

            % 1) 零相位带通滤波
            x_filt = filtfilt(b, a, x);

            % 2) Hilbert 包络提取
            env = abs(hilbert(x_filt));

            % 3) 异常伪迹检测与插值
            z = (env - mean(env)) / std(env);
            bad = abs(z) > 3;
            env(bad) = NaN;
            if any(bad)
                t = 1:num_samples;
                env(bad) = interp1(t(~bad), env(~bad), t(bad), 'linear', 'extrap');
            end

            env_mat(e, :)   = env;
            clean_mat(e, :) = env;  
        end

        emg_envelope_cell{p} = env_mat;
        emg_clean_cell{p}    = clean_mat;
    end
end


% function [emg_clean, emg_envelope] = preprocess_emg(emg_raw, Fs)
%     % ---------------------------------------------------------------
%     % Preprocess EMG signal:
%     % 1. Bandpass filtering (20–450 Hz or adjusted to Fs)
%     % 2. Envelope extraction using Hilbert transform
%     % 3. Artifact removal with z-score + linear interpolation
%     % ---------------------------------------------------------------
% 
%     % Adjust high cutoff if sampling rate is too low
%     if Fs <= 900
%         high_cutoff = Fs / 2 - 1;
%     else
%         high_cutoff = 450;
%     end
% 
%     % 1. Bandpass filtering
%     [b, a] = butter(4, [20 high_cutoff] / (Fs / 2), 'bandpass');
%     emg_filtered = filtfilt(b, a, emg_raw')';
% 
%     % 2. Envelope extraction (Hilbert transform)
%     emg_envelope = abs(hilbert(emg_filtered')');
% 
%     % 3. Artifact removal + interpolation
%     emg_clean = emg_envelope;
%     for i = 1:size(emg_envelope, 1)
%         % Z-score threshold
%         z = (emg_envelope(i,:) - mean(emg_envelope(i,:))) / std(emg_envelope(i,:));
%         outliers = abs(z) > 3;
%         emg_clean(i, outliers) = NaN;
% 
%         % Interpolation over NaNs
%         x = 1:length(emg_clean(i,:));
%         nan_idx = isnan(emg_clean(i,:));
%         if any(nan_idx)
%             emg_clean(i, nan_idx) = interp1(x(~nan_idx), emg_clean(i, ~nan_idx), x(nan_idx), 'linear', 'extrap');
%         end
%     end
% end
% 
