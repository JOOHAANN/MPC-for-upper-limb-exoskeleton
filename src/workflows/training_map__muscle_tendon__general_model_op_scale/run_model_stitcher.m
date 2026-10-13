function final_signal = run_model_stitcher(tasks, combination_raw, baseline_file)
    % RUN_MODEL_STITCHER 根据输入的任务列表，自动计算MAP、插值并缝合到时间轴上
    % Inputs:
    %   tasks: 结构体数组，包含P, Start, Speed, Weight等信息
    %   total_len: 最终信号的总长度
    %   baseline_file: 传递给MAP_model的文件名

    % 2. 循环处理每一个片段
    for i = 1:length(tasks)
        task = tasks(i);
        fprintf('Processing Segment %d: Speed=%d, Weight=%d...\n', i, task.speed, task.weight);

        % --- A. 运行 MAP 模型 ---
        % 注意：这里调用你原本的 MAP_model
        [~, ~, MEP_out] = MAP_model(task.speed, task.weight, 'SM', baseline_file);

        % --- B. 计算原始激活度 (简化版公式) ---
        % 假设 MEP 是 4 列
        raw_val = (MEP_out(:, 1) + MEP_out(:, 2)) * 0.25 + ...
                  (MEP_out(:, 3) + MEP_out(:, 4)) * 0.25;

        % --- C. 全局归一化 ---
        % 使用辅助函数归一化
        processed_val = normalize_max_min(raw_val, task.norm_max, task.ref_min);

        % --- D. 局部二次归一化 (50:150) ---
        if task.local_norm
            % 提取 50:150 区间 (如果数据不够长要做保护)
            roi = 50:150;
            segment_part = processed_val(roi);

            % 这里原本的逻辑比较复杂，需要参考数据的 Min
            % 我们假设 task.ref_min 或者 task.local_idx 对应的 Reference 值
            % 为了简化，这里演示逻辑，你可能需要传入具体的 ref_segment_min
            % 假设局部归一化的下限就是全局下限
            processed_val(roi) = normalize_max_min(segment_part, max(segment_part), task.ref_min);
        end

        % --- E. 插值 (Interpolation) ---
        % 从原始长度 (假设200) 拉伸到 task.P_length
        original_len = length(processed_val); % 通常是 200
        k = 1:original_len;
        t_target = linspace(1, original_len, task.P_length);

        interpolated_segment = interp1(k', processed_val, t_target');

        % --- F. 缝合 (Stitching) ---
        % 计算在总画布上的结束位置
        end_pos = task.start_idx + task.P_length - 1;

        % 越界检查
        if end_pos > total_len
            warning('Segment %d 长度超出总长度，将被截断！', i);
            valid_len = total_len - task.start_idx + 1;
            combination_raw(task.start_idx:end, 2) = interpolated_segment(1:valid_len);
        else
            % 直接填入数组，无需计算 padding
            combination_raw(task.start_idx:end_pos, 2) = interpolated_segment;
        end
    end
    final_signal = combination_raw(:, 2);
end

function out = normalize_max_min(in, target_max, target_min)
    data_min = min(in);
    data_max = max(in);
    % 防止除以0
    if data_max == data_min
        out = in;
    else
        out = ((in - data_min) ./ (data_max - data_min)) * (target_max - target_min) + target_min;
    end
end
