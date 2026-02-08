function [Y_cell, Y_all] = generate_response(m, n, p, d,real_coef, X_cell, e_cell, varargin)
% GENERATE_RESPONSE 生成响应变量
% 输入:
%   m         - 客户端数量
%   n         - 每个客户端的样本量 (m×1 向量)
%   p         - 样本维度
%   real_coef - 各客户端真实系数 (m×1 元胞数组)
%   X_cell    - 样本数据 (m×1 元胞数组, 每个元胞 n(i)×p 矩阵)
%   e_cell    - 噪声数据 (m×1 元胞数组, 每个元胞 n(i)×1 向量)
%   varargin  - 模式参数:
%               模式1: {'MeanFunction', 'sin'} 或 {'MeanFunction', 'exp'}
%               模式2: {'NewtonRaphson', 'vector'} 或 {'NewtonRaphson', 'matrix'}
%
% 输出:
%   Y_cell - m×1 元胞数组，每个元胞包含 n(i)×1 响应向量
%   Y_all  - 所有样本拼接的 (sum(n)×1) 响应向量

% 解析可变输入参数
main_mode = '';
sub_mode = '';

if length(varargin) >= 2
    main_mode = varargin{1};
    sub_mode = varargin{2};
end

% 初始化输出
Y_cell = cell(m,1);
Y_all = zeros(sum(n),1);
ni = 1;  % 样本索引指针

for i = 1:m
    % 提取当前客户端数据
    nn = n(i);
    e = cell2mat(e_cell(i,1));
    X = cell2mat(X_cell(i,1));
    beta = cell2mat(real_coef(i)); %(p-by-d)

    % 计算响应变量 (调用辅助函数)
    Y = compute_Y(X, beta, e, n(i),p,d,main_mode, sub_mode);

    % 存储结果
    Y_cell(i,1) = mat2cell(Y,size(Y,1),size(Y,2));
    Y_all(ni:ni+nn-1) = Y;
    ni = ni + nn;

end
end

function Y = compute_Y(X, beta, e, n, p,d, main_mode, sub_mode)
% COMPUTE_Y 根据指定模式计算响应变量
% 输入:
%   X        - n×p 样本矩阵
%   beta     - p×1 系数向量
%   e        - n×1 噪声向量
%   main_mode - 主模式字符串 ('MeanFunction' 或 'NewtonRaphson')
%   sub_mode  - 子模式字符串 ('sin','exp','vector','matrix')

% 默认为 NewtonRaphson matrix 模式
if isempty(main_mode) && isempty(sub_mode)
    main_mode = 'NewtonRaphson';
    sub_mode = 'vector';
end

% 执行对应模式的计算
switch main_mode
    case 'MeanFunction'
        switch sub_mode
            case 'linear'
                Y = sum(X' * beta,2) + e;
            case 'sin'
                Y = sum(sin(X' * beta),2) + e;    % 正弦函数模式
            case 'exp'
                Y = sum(1*exp(X' * beta),2) + e;  % 指数函数模式
            case 'fun1'
                Y = fun1(X,e,beta,n);  % 指数函数模式
            case 'fun2'
                Y = fun2(X,e,beta,n);  % 指数函数模式
            case 'ex1'
                Y = ex1(X,e,beta,n);  % 指数函数模式
            case 'ex2'
                Y = ex2(X,e,beta,n);  % 指数函数模式
            otherwise
                error('无效子模式: %s。使用''lieanr'', ''sin'' 或 ''exp''', sub_mode);
        end

    case 'NewtonRaphson'
        switch sub_mode
            case 'vector'
                Y = X' * beta(:) + e;       % 向量计算模式
            case 'matrix'
                Y = zeros(n,1);
                for i = 1:n
                    Y(i) = trace(X(:,(i-1)*d+1:i*d)'* beta) + e(i);
                end

            otherwise
                error('无效子模式: %s。使用 ''vector'' 或 ''matrix''', sub_mode);
        end

    otherwise
        error('无效主模式: %s。使用 ''MeanFunction'' 或 ''NewtonRaphson''', main_mode);
end
end


function Y = ex1(X,e,beta,n)

xb = X'* beta;
xb1 = xb(:,1);
xb2 = xb(:,2);
% m = sin(2*xb1).*exp(2+xb2);
% m = sin(2*xb1)+ exp(xb2-1);
m = sin(2*xb1).*exp(xb2);
Y = m + e;

end

function Y = ex2(X,e,beta,n)

xb = X'* beta;
xb1 = xb(:,1);
xb2 = xb(:,2);
m = 3*xb1./(1+(1+xb2).^2);
% m = xb1./(0.5+(1.5+xb2).^2);
Y = m + e;

end




function Y = fun1(X,e,beta,n)

xb = X'* beta;
xb1 = xb(:,1);
% xb2 = xb(:,2);
m = 3*xb1./(1+(1+xb1).^2);
% m = xb1./(0.5+(1.5+xb2).^2);
Y = m + e;
end

function Y = fun2(X,e,beta,n)

xb = X'* beta;
xb1 = xb(:,1);
xb2 = xb(:,2);
m = 5*xb1./(1+(1+xb2).^2);
% m = xb1./(0.5+(1.5+xb2).^2);
Y = m + e;

end
