function [hX,hY,Ew,Ewx] = KerHatXY(X,Y,beta,w,m,n,p,d,varargin)
% Input
%   M:   n-by-1
%   M1:  n-by-d
%   X: cell with p-by-n
%   Y: cell with n-by-1
% Output:
%   hX: n-by-pd in cell(m,1)
%   hY: n-by-1  in cell(m,1)

if max(size(n)) ==1; n = n*ones(m,1);end
options = varargin{1};
if isempty(varargin), mode = "nope"; else, mode = options{1}; end
% if ismatrix(X), mode = "pool"; end 


switch mode
    case "mave"
        Ew = []; Ewx = [];
        h = 0.7;
        [hX,hY] = HatXYonMAVE (X,Y,beta,m,n,d,h);
    case "nope"
        [M,M1]=KerMeanFuncion(X,Y,beta,-1,'reduction');
        hX = cell(m,1);hY = cell(m,1);
        Ew = KerMeanInvVar(X,Y,beta,w,m,n);
        Ewx = KerMeanXInvVar(X,Y,beta,w,m,n,p,d);
        for i=1:m
            hXi = zeros(n(i),p*d);    hYi = zeros(n(i),1);
            Xi = X{i};Yi = Y{i};
            M1i = M1{i}; Mi = M{i};
            bi = beta{i};
            Ewxi =  Ewx{i};    Ewi =  Ew{i};

            for j = 1:n(i)
                % temp = (Xi(:,j) - 0*Ewxi(j,:)'/Ewi(j))*M1i(j,:); % p-by-d
                % temp = (Xi(:,j) - Ewxi(j,:)'/Ewi(j))*M1i(j,:); % p-by-d
                temp = (Xi(:,j) - Ewxi(j,:)')*M1i(j,:); % p-by-d
                hXi(j,:) = temp(:)'; % 1-by-pd
                hYi(j) = hXi(j,:)*bi(:) + Yi(j) - Mi(j); % scalar
            end
            hX{i} = hXi';
            hY{i} = hYi;
        end
    case "pool"
        Ew = []; Ewx = [];
        h = 0.7;
        % 1. based on MAVE's Kernel; 2. hX,hY: matrix.
        [hX,hY] = poolHatXY (X,Y,beta,m,n,d,h); 
end







end

%% pooled Hat X Y  based on MAVE

function [hX,hY] = poolHatXY (X,Y,beta,m,n,d,h)
% 1. based on MAVE's Kernel; 
% 2. hX,hY: pd-by-n and n-by-1 matrix.
% 
b = beta{1};
x = X'; y  = Y; % n-by-p; n-by-1;
n = size(x,1);    p = size(x,2);
ab = raibi(x, y, b, d, h);
m = ab(1,:); % 1-by-n
m1 = ab(2:end,:); % d-by-n
hx = zeros(n,p*d);    hy = zeros(n,1);
ExRed = zeros(n,p); % control the estimated equation
ExRed = KerMeanRed(x,b,h); % n-by-p;
for j = 1:n
    temp = (x(j,:) - ExRed(j,:))'*m1(:,j)'; % p-by-d
    hx(j,:) = temp(:)'; % pd-by-1;
    hy(j) = hx(j,:)*b(:) + y(j) - m(j);
end

hX = hx'; % pd-by-n
hY = hy;  % n-by-1
end


%% Hat X Y  based on MAVE

function [hX,hY] = HatXYonMAVE (X,Y,beta,m,n,d,h)
% beta: p-by-d; x: n-by-p; X{i}: p-by-n;

hX = cell(m,1);hY = cell(m,1);

for i=1:m
    x = X{i}'; y=Y{i};b = beta{i};
    n = size(x,1);
    p = size(x,2);
    ab = raibi(x, y, b, d, h);
    m = ab(1,:); % 1-by-n
    m1 = ab(2:end,:); % d-by-n
    hx = zeros(n,p*d);    hy = zeros(n,1);
    ExRed = zeros(n,p); % control the estimated equation 
    ExRed = KerMeanRed(x,b,h); % n-by-p;
    for j = 1:n
        temp = (x(j,:) - ExRed(j,:))'*m1(:,j)'; % p-by-d
        hx(j,:) = temp(:)'; % pd-by-1;
        hy(j) = hx(j,:)*b(:) + y(j) - m(j);
    end

    hX{i} = hx'; % pd-by-n
    hY{i} = hy;  % n-by-1
end

end

function ExRed = KerMeanRed(x,beta,h)

n = size(x,1);
p = size(x,2);
ExRed = zeros(n,p);
for i = 1:n 
    
    w = kern(x * beta, h, i); % n-by-1;
    temp = x.*w; % n-by-p;
    ExRed(i,:) = mean(temp); % 1-by-p; 
end

end

%%  Hat X Y  based on KerReg



function Ew = KerMeanInvVar(X,Y,beta,w,m,n)
Ew = cell(m,1);
Ker = KerWeight(X,Y,beta,m,'reduction','bandwidth',[],'kernel','normal');
for i=1:m
    Ewi=zeros(n(i),1);
    wi =  w{i};
    Kerw = Ker{i};
    for j = 1:n(i)
        Kerw(j,j) = 0;
        Ewi(j,1) = sum(Kerw(:,j).* wi(j))/sum(Kerw(:,j));
    end
    Ew{i} = Ewi;
end

end

function Ewx = KerMeanXInvVar(X,Y,beta,w,m,n,p,d)
Ewx = cell(m,1);
Ker = KerWeight(X,Y,beta,m,'reduction','bandwidth',[],'kernel','normal');
for i=1:m
    Ewxi=zeros(n(i),p);
    x = X{i};
    wi =  w{i};
    Kerw = Ker{i};
    for j = 1:n(i)
        Kerw(j,j) = 0;
        temp = sum(Kerw(:,j).*wi.*x')/sum(Kerw(:,j));
        Ewxi(j,:) = temp';
    end
    Ewx{i} = Ewxi;
end
end


%%  mave function

% 1.   standardize x
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

% 2.   matrix power
function ai = matpower(a, alpha)
small = 1e-8; % 小阈值，用于忽略小特征值
p1 = size(a,1);
[eve, eva_mat] = eig(a);
eva = diag(eva_mat); % 提取特征值向量
% R代码中归一化特征向量，但MATLAB的eig已返回单位向量，故跳过归一化步骤
index = find(abs(eva) > small); % 找到非零特征值索引
evai = eva;
evai(index) = eva(index).^alpha; % 对特征值取幂
ai = eve * diag(evai) * eve'; % 重构矩阵
end

% 3.   compute kernel
function w1 = kern(x, h, i)
if ismatrix(x) && size(x,2) > 1
    r = size(x,2); % 预测变量维度
    xi = x(i,:); % 第i个观测点
    del = x - xi; % 差值矩阵，相当于R的t(t(x)-x[i,])
    ndel = sum(del.^2, 2); % 欧氏距离平方，相当于R的diag(del%*%t(del))
else
    r = 1;
    del = x - x(i); % 向量情况
    ndel = del.^2;
end
% 计算高斯权重
w = (1/sqrt(2*pi))^r * (1/h^r) * exp(-0.5 * ndel / h^2);
w = w(:); % 确保为列向量
w1 = w / sum(w); % 归一化权重
end

% 4. weighted least squares coefficients
function b = wls(x, y, w)
% 加权设计矩阵：相当于R的x*w（元素乘，w广播）
x_weighted = x .* w;
xtwx = x_weighted' * x; % 相当于R的t(x*w)%*%x
xtwy = x_weighted' * y; % 相当于R的t(x*w)%*%y
b = matpower(xtwx, -1) * xtwy; % 使用matpower求逆
end

% 5. 计算ai bi，给定beta
function ab = aibi(x, y, beta, h)
n = size(x, 1);
p = size(x, 2);
ab = [];

for i = 1:n
    % 构建设计矩阵：第一列为1，其余为(x_j - x_i)*beta
    x_diff = x - x(i,:);  % 相当于R的t(t(x) - x[i,])
    del = [ones(n,1), x_diff * beta];  % 添加截距项

    % 计算核权重，带宽调整与R一致
    w = kern(x, h * n^(-1/(p+4)), i);

    % 加权最小二乘计算ai和bi
    abi = wls(del, y, w);
    ab = [ab, abi];  % 相当于R的cbind
end

end

% 6. 计算beta，给定ai和bi
function beta_mat = bb(x, y, h, ab)
n = size(x, 1);
p = size(x, 2);
d = size(ab, 1) - 1;  % beta的维度

ai = ab(1, :);  % 提取ai（第一行）
bi = ab(2:d+1, :);  % 提取bi（剩余行）

w = [];
yma = [];
bxmx = [];

for i = 1:n
    % 计算y - ai[i]
    yma = [yma; y - ai(i)];

    % 计算核权重
    w = [w; kern(x, h * n^(-1/(p+4)), i)];

    % 计算(x_j - x_i)与bi[,i]的Kronecker积
    x_diff = x - x(i,:);  % 相当于R的t(t(x)-x[i,])
    kronecker_prod = kron(bi(:,i)', x_diff);  % 相当于R的kronecker
    bxmx = [bxmx; kronecker_prod];
end

% 加权最小二乘求解beta
beta_vec = wls(bxmx, yma, w);
beta_mat = reshape(beta_vec, p, d);  % 重塑为p×d矩阵
end

% 7. 计算目标函数
function objec = obj(x0, y, ab, beta, h)
n = size(x0, 1);
p = size(x0, 2);
x = stand(x0);  % 标准化x
d = size(ab, 1) - 1;

ai = ab(1, :);
bi = ab(2:d+1, :);

objec = 0;

for i = 1:n
    % 计算残差：y - ai[i] - (x - x[i,])*beta*bi[,i]
    x_diff = x - x(i,:);
    residual = y - ai(i) - sum(x_diff * beta * bi(:,i), 2);

    % 计算加权平方残差
    kernel_weight = kern(x, h * n^(-1/(p+4)), i);
    weighted_residual = (residual.^2) .* kernel_weight;

    objec = objec + sum(weighted_residual);
end
end


% 8. MAVE函数
function eta = mave(x, y, d, h, beta0, nmave)
z = stand(x);
v = matpower(cov(x), 0.5) * beta0;

for i = 1:nmave
    ab = aibi(z, y, v, h);
    v = bb(z, y, h, ab);
    % 可选：打印目标函数值（与R代码注释部分对应）
    % fprintf('Objective: %f\n', obj(z, y, ab, v, h));
end

eta = matpower(cov(x), -0.5) * v;
end

% 9. OPG函数
function v3 = opg(x, y, d, h)
n = size(x, 1);
p = size(x, 2);
z = stand(x);
b = [];

for i = 1:n
    % 构建设计矩阵：第一列为1，其余为(z_j - z_i)
    z_diff = z - z(i,:);  % 相当于R的t(t(z) - z[i,])
    del = [ones(n,1), z_diff];  % 相当于R的cbind(1, ...)

    % 计算核权重，带宽调整与R一致
    w = kern(z, h * n^(-1/(p+4)), i);

    % 加权最小二乘，取第2到p+1个系数（去掉截距项）
    bi_full = wls(del, y, w);
    bi = bi_full(2:p+1);  % 相当于R的[2:(p+1)]
    b = [b, bi];  % 水平连接，相当于R的cbind
end

% 计算b*b'的特征向量，取前d个主方向
[eve, ~] = eigs(b * b');
v1 = eve(:, 1:d);  % 取前d个特征向量

% 调整协方差矩阵的影响
v3 = matpower(cov(x), -0.5) * v1;
end

% 10. RMAVE函数中的raibi函数
function ab = raibi(x, y, beta, d, h)
n = size(x, 1);
p = size(x, 2);
ab = [];

for i = 1:n
    % 构建设计矩阵：第一列为1，其余为(x_j - x_i)*beta
    x_diff = x - x(i,:);
    del = [ones(n,1), x_diff * beta];  % 相当于R的cbind(1, t(t(x)-x[i,])%*%beta)

    % 计算核权重，在降维空间x*beta上计算
    w = kern(x * beta, h, i);  % 相当于R的kern(x%*%beta, h, i)

    % 加权最小二乘计算系数
    abi = wls(del, y, w);
    ab = [ab, abi];  % 水平连接
end
end

% refined beta
function beta_mat = rbb(x, y, h, ab, v)
n = size(x, 1);
p = size(x, 2);
d = size(ab, 1) - 1;

ai = ab(1, :);
bi = ab(2:d+1, :);

w = [];
yma = [];
bxmx = [];

for i = 1:n
    yma = [yma; y - ai(i)];
    % 在降维空间x*v上计算核权重，带宽调整使用d而不是p
    w = [w; kern(x * v, h * n^(-1/(d+4)), i)];

    x_diff = x - x(i,:);
    kronecker_prod = kron(bi(:,i)', x_diff);
    bxmx = [bxmx; kronecker_prod];
end

beta_vec = wls(bxmx, yma, w);
beta_mat = reshape(beta_vec, p, d);
end

% refined MAVE
function eta = rmave(x, y, d, h, beta0, nmave)
z = stand(x);
v = matpower(cov(x), 0.5) * beta0;
u = v;

for i = 1:nmave
    ab = raibi(z, y, u, d, h);
    v = rbb(z, y, h, ab, u);
    u = v;
    % 可选：打印目标函数值（与R代码注释部分对应）
    % fprintf('Objective: %f\n', obj(z, y, ab, v, h));
end

eta = matpower(cov(x), -0.5) * v;
end

% sample score function
function res = sco(x, y, d, h, beta)
p = size(x, 2);
n = size(x, 1);
z = stand(x);
v = matpower(cov(x), 0.5) * beta;
ab = raibi(z, y, v, h);
ai = ab(1, :);
bi = ab(2:d+1, :);

res = zeros(d, p);

for i = 1:n
    res = res + bi(:, i) * z(i, :) * (y(i) - ai(i));
end

res = res' / n;
end