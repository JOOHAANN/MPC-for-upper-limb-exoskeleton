function final_signal = run_model_stitcher(tasks, Combination_raw, baseline_file)
    % RUN_MODEL_STITCHER 根据输入的任务列表，自动计算MAP、插值并缝合到时间轴上
    % Inputs:
    %   tasks: 结构体数组，包含P, Start, Speed, Weight, norm_max, ref_min 等信息
    %   Combination_raw: 初始画布
    %   baseline_file: 传递给MAP_model的文件名

    % 1. 初始化一个巨大的全零画布
    total_len = length(Combination_raw);
    final_signal = Combination_raw;

    % 2. 循环处理每一个片段
    for i = 1:length(tasks)
        task = tasks(i);
        fprintf('Processing Segment %d: Speed=%d, Weight=%d...\n', i, task.speed, task.weight);

        % --- A. 运行 MAP 模型 ---
        [~, ~, MEP_out] = MAP_model(task.speed, task.weight, 'SM', baseline_file);

        % --- B. 计算原始激活度 (简化版公式) ---
        raw_val = (MEP_out(:, 1) + MEP_out(:, 2)) * 0.25 + ...
                  (MEP_out(:, 3) + MEP_out(:, 4)) * 0.25;

        % --- C. 找到原信号的关键特征点 (峰顶 & 谷底) ---
        P_old = max(raw_val);

        TF_min = islocalmin(raw_val);
        local_mins = raw_val(TF_min);
        if isempty(local_mins)
            % 如果找不到局部低谷，为了程序不报错，临时用最小值代替
            warning('Segment %d 未检测到局部低谷，使用全局最小值代替', i);
            V_old = min(raw_val);
        else
            V_old = min(local_mins);
        end

        % --- D. 执行分段数值映射 (Piecewise Mapping) ---
        % 核心修改：不再使用线性公式 (y=ax+b)，而是使用映射表
        % 目的：锁死0点，防止负数出现，同时拉伸峰值和谷值

        % 1. 定义原信号的锚点 [基线0, 原谷值, 原峰值]
        min_raw = min(raw_val);
        map_old = [min_raw, V_old, P_old];

        % 2. 定义目标信号的锚点 [基线0, 目标谷值, 目标峰值]
        map_new = [0, task.ref_min, task.norm_max];

        % 3. 处理原信号可能存在的微小负数噪声 (防止 interp1 报错或映射失效)
        min_raw = min(raw_val);
        if min_raw < 0
            % 如果原信号有负值，把负值强行映射到 0
            map_old = [min_raw, map_old];
            map_new = [0, map_new];
        end

        % 4. 排序和去重 (interp1 的硬性要求：X轴必须单调递增)
        [map_old, sortIdx] = sort(map_old);
        map_new = map_new(sortIdx);
        [map_old, uniqueIdx] = unique(map_old);
        map_new = map_new(uniqueIdx);

        % 5. 应用映射 (使用 'pchip' 保持曲线光滑连续)
        % 最后的 0 表示：如果遇到范围以外的值，默认给 0
        raw_val_rescaled = interp1(map_old, map_new, raw_val, 'pchip', 0);

        % --- E. 插值 (Interpolation - 长度拉伸) ---
        % 从原始长度拉伸到 task.P_length
        original_len = length(raw_val_rescaled);
        k = 1:original_len;
        t_target = linspace(1, original_len, task.P_length);

        interpolated_segment = interp1(k', raw_val_rescaled, t_target');

        % --- F. 缝合 (Stitching) ---
        end_pos = task.start_idx + task.P_length - 1;

        if end_pos > total_len
            warning('Segment %d 长度超出总长度，将被截断！', i);
            valid_len = total_len - task.start_idx + 1;
            final_signal(task.start_idx:end) = interpolated_segment(1:valid_len);
        else
            final_signal(task.start_idx:end_pos) = interpolated_segment;
        end
    end
end
