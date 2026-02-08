function [beta_estimated,objective_value] = ...
        opt_poolDHDR_orth(beta_initial,X,Y,m,n,p,d,T_max,K_max,lambda,varargin)
% varargin: 'original', 'cosine' or 'invex'
% - Input:
%   - beta_initial: p-by-d matrix, collected as m-by-1 cell.
%   - X: p*d-by-n matrix, collected as m-by-1 cell. 尽量用 'vector'
%   - Y: n-by-1 matrix, collected as m-by-1 cell. 
%   - m,n,p,d: the number of nodes, samples, original dimension and
%              reduction dimensions
%   - T_max: communication round; 
%   - K_max: local update rounds.
% - Output:
%   - beta_estimated: p-by-d matrix, collected as m-by-1 cell.
%   - objective_value: the objective in each communication round

mode = varargin{1};
modeL = mode{1};modeP = mode{2};modegrad = mode{3};

objective_value = zeros(T_max+1,1);
beta_estimated=beta_initial;
objGiNew = zeros(m,1);

% [objective_value(1),objGi] = CostEachNode(beta_initial,X,Y,m,n,lambda,mode);

[temp ,~]= obj_DHDR(beta_initial{1},beta_initial,m,sum(n),Y,X,lambda,mode);
objGi = repmat(temp,m,1);
objective_value(1) = sum(objGi);


fprintf('初始函数值为%.4f \n',objective_value(1));

step = 0.005;
aclr = 0.9;

coef_v = cell(m,1);coef_intial = cell(m,1);coef_update = cell(m,1);


for T = 1:T_max
    cost_temp=zeros(m,1);
    beta_initial = beta_estimated;
    
    for mi=1
        
        Xi=X;
        Yi=Y;
        n_i = sum(n);
        alpha = cell2mat(beta_initial(mi));        
        beta_cell = beta_initial; %在每一次通信内部 beta_cell 不变
        % optimalization part on local client------------------------------
        for K = 1:K_max
            
            if K>1; alpha = coef_update{mi}; end
            
            beta_cell(mi) = mat2cell(alpha,p,d);
            [~ ,grad]= obj_DHDR(alpha,beta_cell,m,n_i,Yi,Xi,lambda,mode);
            Walpha = grad*alpha' - alpha*grad';
            Wupdate = (eye(p)+step*Walpha/2)^(-1)*(eye(p)-step*Walpha/2);
            alphanew = Wupdate*alpha;
            coef_update(mi) = mat2cell(alphanew,p,d);
            
        end
        % beta_cell(mi) = mat2cell(alphanew,p,d);
        [objGiNew(mi) ,~]= obj_DHDR(alphanew,beta_cell,m,n_i,Yi,Xi,lambda,mode);
        % feadback the optimalization result on local client---------------
        beta_estimated(mi) = mat2cell(alphanew,p,d);
        
    end
    
    for mi=2:m
        objGiNew(mi) = objGiNew(1);
        beta_estimated{mi} = beta_estimated{1};
        coef_update{mi} = coef_update{1};
    end




    % At the server, Broadcast the value of the objective function---------
    if  mod(T,200)==0
        fprintf('在轮次%d时, 目标函数值cost为 %f \n',T,sum(objGiNew));
    end
    % At the server, determine whether to exit the loop--------------------
    if abs(sum(objGiNew)-sum(objGi))<10^(-5)
        fprintf('在轮次%d时, 目标函数值cost为 %f \n',T,sum(objGiNew));
        fprintf('算法收敛,总共通讯了%d轮\n',T);
        objective_value(T+1:end) = sum(objGiNew);% 算法收敛, 记录最终的目标函数值
%         objective_value(T+20:end) = [];
        break;
    elseif sum(objGiNew)-sum(objGi)>0.03 && T >= 2
        fprintf('损失函数增加,算法收敛,总共通讯了%d轮\n',T);
        fprintf('在轮次%d时, 目标函数值cost为 %f \n',T,sum(objGi)); % 保留更新之前的结果
        beta_estimated = beta_initial;
        objective_value(T+1:end) =  sum(objGi);
        % objective_value(T+20:end) = [];
        break;
    else
        % continue loop-----------
        objGi=objGiNew;
        objective_value(T+1) = sum(objGi); % record objective function value
    end
    % The conmunication round --------------------------------------------- 
    
end




end

function [objG,objGi] = CostEachNode(beta,X,Y,m,n,lambda,varargin)
mode = varargin{1};
objGi = zeros(m,1);

for i = 1:m
    Xi=cell2mat(X(i));
    Yi=cell2mat(Y(i));
    alpha = cell2mat(beta(i));
    objGi(i) =  obj_DHDR(alpha,beta,m,n(i),Yi,Xi,lambda,mode);
end

objG = sum(objGi);

end
