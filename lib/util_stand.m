function [X_cell]  = util_stand(X,m)

X_cell = cell(m,1);Y_cell = cell(m,1);

for i = 1:m
    x = X{i,1}'; 
    x = stand(x);
    X_cell{i,1} = x';
end




end 


function z = stand(x)
    if ismatrix(x) && size(x,2) > 1
        n = size(x,1);
        p = size(x,2);
        xb = mean(x,1); % 计算每列的均值，相当于R的apply(x,2,mean)
        xb_mat = repmat(xb, n, 1); % 复制均值向量以匹配x的维度
        x1 = x - xb_mat; % 减去均值
        sigma = (x1' * x1) / (n-1); % 计算协方差矩阵，相当于R的t(x1)%*%x1/(n-1)
        [eve, eva_mat] = eig(sigma); % 特征分解，eva_mat是特征值矩阵
        eva = diag(eva_mat); % 提取特征值向量
        % 计算Sigma^(-1/2)，注意MATLAB的eig返回单位特征向量，无需额外归一化
        sigmamrt = eve * diag(1./sqrt(eva)) * eve'; 
        z = (sigmamrt * x1')'; % 标准化并转置为n×p矩阵，相当于R的t(z)
    else
        % 处理向量输入
        z = (x - mean(x)) / std(x);
    end
end