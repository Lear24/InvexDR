function data_winsor = util_winsorize(data, lower_pct, upper_pct, dim)
% 功能: 对数据进行缩尾处理（Winsorization）
% 输入:
%   data      - 输入数据（向量或矩阵）
%   lower_pct - 下分位数阈值（默认 1%）
%   upper_pct - 上分位数阈值（默认 99%）
%   dim       - 处理维度（1=按列处理，2=按行处理，默认按列）
% 输出:
%   data_winsor - 缩尾后的数据

% 参数检查与默认值
if nargin < 2, lower_pct = 1; end
if nargin < 3, upper_pct = 99; end
if nargin < 4, dim = 1; end % 默认按列处理

% 计算分位数阈值
lower_bound = prctile(data, lower_pct, dim);
upper_bound = prctile(data, upper_pct, dim);

% 替换极端值
data_winsor = data;
if dim == 1 % 按列处理
    for i = 1:size(data, 2)
        col_data = data(:, i);
        col_data(col_data < lower_bound(i)) = lower_bound(i);
        col_data(col_data > upper_bound(i)) = upper_bound(i);
        data_winsor(:, i) = col_data;
    end
else % 按行处理
    for i = 1:size(data, 1)
        row_data = data(i, :);
        row_data(row_data < lower_bound(i)) = lower_bound(i);
        row_data(row_data > upper_bound(i)) = upper_bound(i);
        data_winsor(i, :) = row_data;
    end
end
end