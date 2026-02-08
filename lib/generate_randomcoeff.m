function [real_coef,real_coef_cell ,real_theta_para]=generate_randomcoeff(m,p,d,theta_1,theta_2,varargin) 
% [real_beta, real_theta_para]=generate_randomcoeff(m,d,theta_1,theta_2,varargin)
% Input:
% m - client number
% p - original dimension  
% d - reduced dimension  
% theta_1 - similarity
% theta_2 -similarity
% varargin:
% mode binomial: 0-1
% mode set: given coeff
% Output:
% real_coef: 


omega_1=binornd(1,0.5,p,1); % d*1 vector
real_coef=zeros(m,p); % 我们的真实参数矩阵是以 m*d 的矩阵储存的
real_coef(1,:)=omega_1';
I=eye(p);
real_theta_para = zeros(m,1);



for i = 2:m
    error_contol=1;
    rot_axis = binornd(1,0.5,p,1); %旋转轴 - n 维关于旋转轴的旋转不是唯一的 我们可以有不同的确定旋转平面的方式
    tau= theta_1 + (theta_2 - theta_1)*rand(1,1);
    tau = tau*(2*binornd(1,0.5,1,1)-1);
    %     disp(cos(tau));
    %     T=eye(d);
    %     T(2,2)=cos(tau);T(2,3)=-sin(tau);
    %     T(3,2)=sin(tau);T(3,3)=cos(tau);
    %     A=[rot_axis,I];
    A=[omega_1,rot_axis,I];
    %     A=[omega_1,I]; % 这样也可以产生旋转了tau角度的向量,但是总是绕着同一个轴,无法引入更多的随机性. 至少应该随机取d-2个平面作为 旋转轴
    [~,J]=rref(A); % J 中为A的极大线性无关组的坐标信息
    %     disp(rank(A)==d);
    %     disp(J(1)==1);
    Q = util_gram_schmidt(A(:,J)); % 正交化
    real_coef(i,:)= util_rotmnd(Q(:,3:p),tau)*omega_1;
    
    while error_contol==1
        if (real_coef(1,:)*real_coef(i,:)'/norm(real_coef(1,:),2)/norm(real_coef(i,:),2)-cos(tau))<10^(-4)
            fprintf('设备 %d 与原向量角度为 %1.4f Pi \n',i,tau/pi)
            real_theta_para(i) = tau/pi;
            error_contol=0;
        else
            fprintf('设备%d 无法确定旋转角度, 建议重新生成\n',i)
            rot_axis = binornd(1,0.5,p,1);
            tau= theta_1 + (theta_2 - theta_1)*rand(1,1);
            A=[omega_1,rot_axis,I];
            [~,J]=rref(A); % J 中为A的极大线性无关组的坐标信息
            Q = gram_schmidt(A(:,J)); % 正交化
            real_coef(i,:)= util_rotmnd(Q(:,3:p),tau)*omega_1;
        end
        
        
    end 
end


real_coef_cell = mat2cell(real_coef,repelem(1,m));





end