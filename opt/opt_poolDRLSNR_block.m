function [beta_estimated,trsimi,objall,fnorm,stop] = ...
    opt_poolDRLSNR_block(beta_initial,X,Y,X_all_matrix,Y_all_matrix,m,n,p,d,T_max,K_max,lambda,N_max,real_beta,varargin)
% stop: 0 - max iter; 1 - convergence; (ignore the increase case )
optionsOPT = varargin{1};
optionsInvVar = varargin{2};
optionsKer = varargin{3};
trsimi = zeros(N_max+1,m);
objall = zeros(N_max+1,T_max+1);
N=1; stop =0;
beta_initial = BlockCoef(beta_initial,m,d);

for i = 2:m
    beta_initial{i} = beta_initial{1}; % estimate one parameter with all samples.
end 

ts = util_TraceSimilarity(beta_initial,real_beta,m,d);
f = util_Fnorm(beta_initial, real_beta,m,d);
trsimi(N,:) = ts';
fnorm(N,:) = f;

mm = 1; nn = sum(n);
while (N<=N_max) && (stop ==0)
    stop =0;

    fprintf('通信轮次%d时,',N);

    
    N = N + 1;


    InvVar = KerInvVar(X_all_matrix,Y_all_matrix,beta_initial,mm,nn,optionsInvVar); % mm = 1; nn = N (sum n); 
    [hX,hY,~,~] = KerHatXY(X_all_matrix,Y_all_matrix,beta_initial,InvVar,m,n,p,d,optionsKer); % pd-by-nn; nn-by-1;

    [beta_estimated,objective_value] = ...
        opt_poolDHDR_block(beta_initial,hX,hY,m,n,p,d,T_max,K_max,lambda,optionsOPT);
    % 记录机制
    beta_initial = beta_estimated;
    ts = util_TraceSimilarity(beta_estimated, real_beta,m,d);
    f = util_Fnorm(beta_estimated, real_beta,m,d);
    trsimi(N,:) = ts';
    fnorm(N,:) = f;
    objall(N,:) = objective_value';
    fprintf('----------------------------------------------------------\n\n');
    % 退出机制
    if (abs(objall(N-1,1) - objall(N,1) ) < 0.0001)
        fprintf('在通信轮次%d时, 目标函数值cost为 %f \n',N,objall(N,1));
        fprintf('算法收敛,总共通讯了%d轮---------------------------------\n\n\n',N);
        stop = 1;
        trsimi(N+1:end,:) = [];
        objall(N+1:end,:) = [];
        fnorm(N+1:end,:) = [];
        break;
    elseif  (objall(N,1) - objall(N-1,1) >0.01) && N>3 %若目标函数值增加
        fprintf('在通信轮次%d时, 损失函数增加,目标函数值cost为 %f \n',N,objall(N-1,1));
        fprintf('算法收敛,总共通讯了%d轮---------------------------------\n\n\n',N-1);
        stop = 2;
        trsimi(N+1,:) = trsimi(N-1,:);
        objall(N+1,:) = objall(N-1,:);
        fnorm(N+1,:) = fnorm(N-1,:);
        trsimi(N+2:end,:) = []; %还是保留增加了的 看看情况
        objall(N+2:end,:) = [];
        fnorm(N+2:end,:) = [];
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


function fnorm = util_Fnorm(beta, real_beta,m,d)

fnorm = zeros(m,1);
for mi = 1: m
    A = util_projection(beta{mi}); B = util_projection(real_beta{mi});
    fnorm(mi) = norm(A-B,'fro');
end
end