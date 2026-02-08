function [real_coef, real_theta_para]=generate_RandomCoef(m,p,d,theta_1,theta_2,varargin)
% Input:
% m - node number
% p - variable dimension
% d - reduction dimension
% theta - similarity of angle
%
% varargin:
% 1. mode
%    rot
%    set
% 2. set: 1 2 3
% 3. rot: binomial normal sparse normal
% Output:
% real_coef - real central mean subspace  as p-by-b matrix
% real_theta_para the real cosine value between the corresponding eigenvector of central mean subspace on each node.
%
% 使用说明:
% 1. p d 很接近 随机生成 binomial 分布可能是不满秩的
%
%


mode = 'binomial'; %默认为向量计算形式


if ~isempty(varargin)                                                      % Check if there are any optional parameters to input
    input_mode = varargin{1};                                              % Extract the first optional parameter
    if ischar(input_mode)                                                  % Verify whether the input is a string
        mode = lower(input_mode);                                          % Convert to lowercase to avoid case sensitivity issues
    end
end


switch mode
    case 'coef1'
        d = 2;
        beta_1 = zeros(p,d);
        beta_1(1:4,1) = [1,0,1,1]; beta_1(1:4,2) = [0,1,-1,1];
        fprintf('In this setting beta on node 1 is an %d-by-%d matrix \n',p,d)
    case 'coef2'
        d = 2;
        beta_1 = zeros(p,d);
        beta_1(3:6,1) = [1,0,1,1]; beta_1(3:6,2) = [0,1,-1,1]; %2-1
        % beta_1(3:6,1) = [1,0,0,1]; beta_1(3:6,2) = [0,1,0,-1]; %2-2
        % beta_1(3:6,1) = [1,0,0,2]; beta_1(3:6,2) = [0,2,0,-1]; %2-5
        % beta_1(3:6,1) = [2,-2,0,2]; beta_1(3:6,2) = [0,2,2,-1]; %2-6
        % beta_1(3:6,1) = [1,0,0,1]; beta_1(3:6,2) = [0,0,1,-1]; %2-4
        % beta_1(3:6,1) = [1,0,-1,1]; beta_1(3:6,2) = [0,1,1,1]; %2-3
        fprintf('In this setting beta on node 1 is an %d-by-%d matrix \n',p,d)
    case 'binomial'
        beta_1=binornd(1,0.5,p,d); % d*1 vector
        beta_1 = util_gram_schmidt(beta_1);
    case 'normal'
        beta_1=randn(1,0.5,p,d); % d*1 vector
        beta_1 = util_gram_schmidt(beta_1);
    case 'setting1'
        d = 1;
        beta_1 = zeros(p,d);
        beta_1(1:4,1) = [1,1,-1,1];
        fprintf('In this setting beta on node 1 is an %d-by-%d matrix \n',p,d)
    case 'setting2'
        d = 2;
        beta_1 = zeros(p,d);
        beta_1(1:4,1) = [1,0,1,1];beta_1(1:4,2) = [0,1,-1,1];
        fprintf('In this setting beta on node 1 is an %d-by-%d matrix \n',p,d)
    case 'setting3'
        d = 2;
        beta_1 = zeros(p,d);
        beta_1(end-3:end,1) = [1,0,1,1];beta_1(end-3:end,2) = [0,1,-1,1];
        fprintf('In this setting beta on node 1 is an %d-by-%d matrix \n',p,d)
end

if theta_1 == theta_2 && theta_2 == 0
    real_coef=cell(m,1);
    for i=1:m
        real_coef{i} = beta_1;
        real_theta_para = ones(m,d);
    end
    return
end


real_coef=cell(m,1); % 我们的真实参数矩阵是以 m*d 的矩阵储存的
real_coef(1)=mat2cell(beta_1,p,d);
I=eye(p);
real_theta_para = zeros(m,d);
real_theta_para(1,:) = ones(1,d);
beta_1 = util_gram_schmidt(beta_1);


for i = 2:m
    beta_i = zeros(p,d);
    % error_contol=1;
    rot_axis = binornd(1,0.5,p,1); %旋转轴 - n 维关于旋转轴的旋转不是唯一的 我们可以有不同的确定旋转平面的方式
    tau= theta_1 + (theta_2 - theta_1)*rand(1,1);
    tau = tau*(2*binornd(1,0.5,1,1)-1);
    for j = 1:d
        error_contol=1;
        A=[beta_1(:,j),rot_axis,I];
        [~,J]=rref(A);
        Q = util_gram_schmidt(A(:,J));
        beta_i(:,j)= util_rotmnd(Q(:,3:p),tau)*beta_1(:,j);

        while error_contol==1
            if abs(beta_1(:,j)'* beta_i(:,j)/norm( beta_i(:,j),2)/norm(beta_1(:,j),2)-cos(tau))<10^(-3)
                fprintf('the cosine value of %d -th eigenvector between node 1 and %d: %1.4f Pi \n',j,i,tau/pi)
                real_coef(i)=mat2cell(beta_i,p,d);
                % real_theta_para(i,j) = tau/pi;
                real_theta_para(i,j) = cos(tau);
                error_contol=0;
            else
                fprintf('node %d, the %d-th eigenvector : Regenetate\n',i,i)
                rot_axis = binornd(1,0.8,p,1);
                tau= theta_1 + (theta_2 - theta_1)*rand(1,1);
                A=[beta_1(:,j),rot_axis,I];
                [~,J]=rref(A);
                Q = util_gram_schmidt(A(:,J));
                beta_i(:,j)= util_rotmnd(Q(:,3:p),tau)*beta_1(:,d);
            end

        end
    end






end

end