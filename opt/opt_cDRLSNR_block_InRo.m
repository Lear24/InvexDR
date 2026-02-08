function [beta_estimated,trsimi,objall,results,stop] = ...
    opt_cDRLSNR_block_InRo(beta_initial,X,Y,m,n,p,d,T_max,K_max,lambda,N_max,real_beta,varargin)
% stop: 0 - max iter; 1 - convergence; (ignore the increase case )
% initialization robustness

if ~isempty(varargin{1}); optionsOPT = varargin{1}; else; optionsOPT = 'nope'; end
if ~isempty(varargin{2}); optionsInvVar = varargin{2}; else; optionsInvVar = 'nope'; end
if ~isempty(varargin{3}); optionsInRo = varargin{3}; else; optionsInvVar = 'nope'; end


if any ([strcmp(optionsOPT,'nope'), strcmp(optionsInvVar,'nope'), strcmp(optionsInvVar,'nope')])
    error('实验设置输入不足')
end

% 实验设置部分: 选取需要进行实验的结果
tRound = optionsInRo.tRound;
IniNum = optionsInRo.IniNum;
IniOption = optionsInRo.IniOption;


trsimi = zeros(N_max+1,m);
objall = zeros(N_max+1,T_max+1);
N=1; stop =0;
beta_initial = BlockCoef(beta_initial,m,d);
ts = util_TraceSimilarity(beta_initial,real_beta,m,d);
trsimi(N,:) = ts';

while (N<=N_max) && (stop ==0)
    stop =0;
    fprintf('通信轮次%d时,',N);
    N = N + 1;


    InvVar = KerInvVar(X,Y,beta_initial,m,n,optionsInvVar);
    [hX,hY,~,~] = KerHatXY(X,Y,beta_initial,InvVar,m,n,p,d);
    [beta_estimated,objective_value] = ...
        opt_cDHDR_block(beta_initial,hX,hY,m,n,p,d,T_max,K_max,lambda,optionsOPT);


    if N-1 == tRound
        results = test_InitialRobust(hX,hY,real_beta,IniNum,IniOption,...
            m,n,p,d,T_max,K_max,lambda,optionsOPT); % 在第tRound-th round 进行初值稳定性的实验
    end



    % 记录机制
    beta_initial = beta_estimated;
    ts = util_TraceSimilarity(beta_estimated, real_beta,m,d);
    trsimi(N,:) = ts';
    objall(N,:) = objective_value';
    fprintf('----------------------------------------------------------\n\n');
    % 退出机制
    if (abs(objall(N-1,1) - objall(N,1) ) < 0.0001)
        fprintf('在通信轮次%d时, 目标函数值cost为 %f \n',N,objall(N,1));
        fprintf('算法收敛,总共通讯了%d轮---------------------------------\n\n\n',N);
        stop = 1;
        trsimi(N+1:end,:) = [];
        objall(N+1:end,:) = [];
        break;
    elseif  (objall(N,1) - objall(N-1,1) >0.5) && N>3 %若目标函数值增加
        fprintf('在通信轮次%d时, 损失函数增加,目标函数值cost为 %f \n',N,objall(N-1,1));
        fprintf('算法收敛,总共通讯了%d轮---------------------------------\n\n\n',N-1);
        stop = 2;
        trsimi(N+1:end,:) = []; %还是保留增加了的 看看情况
        objall(N+1:end,:) = [];
        break;
    end



end

end

function beta = BlockCoef(beta,m,d)

for i = 1:m
    b = beta{i};
    b(1:d,1:d)=eye(d);
    beta{i} = b;
end

end

function ts = util_TraceSimilarity(beta, real_beta,m,d)

ts = zeros(m,1);
for mi = 1: m
    A = util_projection(beta{mi}); B = util_projection(real_beta{mi});
    ts(mi) = trace(A*B)/d;
end
end



function  results = test_InitialRobust(hX,hY,real_beta,IniNum,IniOption,m,n,p,d,T_max,K_max,lambda,optionsOPT)


beta = cell(IniNum,1);
objvalue = zeros(T_max+1,IniNum);
trsimi = zeros(m,IniNum);

for i = 1:IniNum
    beta_initial=set_initialization(m,p,d,IniOption{1},IniOption{2});
    beta_initial = BlockCoef(beta_initial,m,d);
    [beta_estimated,objective_value] = ...
        opt_cDHDR_block(beta_initial,hX,hY,m,n,p,d,T_max,K_max,lambda,optionsOPT);
    ts = util_TraceSimilarity(beta_estimated, real_beta,m,d);

    beta{i} = beta_estimated;
    objvalue(:,i) = objective_value;
    trsimi(:,i)= ts;
end

results.EstimatedCoef   = beta;
results.ObjectiveValue  = objvalue;
results.TraceSimilarity  = trsimi;
end