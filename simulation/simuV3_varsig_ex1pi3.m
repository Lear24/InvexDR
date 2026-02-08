%% The main experiment with sig varies: ex1-pi3
%% Setting parameters
clear all;
clc
rng('default')
rng(1)
% simulation parameter setting


ex = 'ex1';coef = 'coef1';setting = 'set1'; p=16; d=2; m=10;
% ex = 'ex2';coef = 'coef1';setting = 'set2'; p=10;d=2;  m=5;
nn = 300*ones(m,1);

% 生成真实样本 X
switch setting
    case 'set1'
        % uniform distribution with [-2,2]
        para = [-2 2];optionsGenerateX = {'uniform',para};
    case 'set2'
        % normal distribution with Sigma = (0.5^(k-l))
        para = 0.5; optionsGenerateX = {'normal',para};
end

% 真实参数生成
theta_1=0;theta_2=pi/3;
[real_beta,real_beta_similarity]= ...
    generate_RandomCoef(m,p,d,theta_1,theta_2,coef);

%% Conduct Part

%variation part
repetition = 100;num_var = 3;
% var is the variable varying in the simulation
variable = [0.8,1,1.2];
%storage part
beta_estimated = cell(12,num_var,repetition);
trsimi = cell(12,num_var,repetition);
objall = cell(12,num_var,repetition);
fnorm = cell(12,num_var,repetition);

beta_initial=set_initialization(m,p,d,'binomial','orth');
options = {'vector', 'project', 'grad'};
optionsInvVar = {'Constant',1};

% T_max= 10; K_max = 1;N_max=3;
T_max= 1000; K_max = 1;N_max=2000;

lamALL = 1;
% lamALL = 1.2;
% lamALL = 1.6;

%% 
% load('SimuV2_varsig_ex1pi3_rep100.mat')
for rep = 1:repetition
    for var = 1:num_var

        rng(rep);
        % 生成噪声
        sig= variable(var);  sigma= sig*ones(m,1);
        e = generate_noise(m,nn,sigma);
        % 生成样本
        [X,X_all_matrix] = generation_sample(m,nn,p,'distribution',...
            optionsGenerateX{1},'parameter',optionsGenerateX{2});
        % 生成真实样本Y
        [Y,Y_all_matrix] = ...
            generate_response(m,nn,p,d,real_beta,X,e,'MeanFunction',ex);

        
        %============================  MAVE ====================================
        fprintf("========================== rep = %d var = %d ======================= \n", rep, variable(var))
        % 1. MAVE
        hmave = 0.7; nmave = 10;
        [~, beta_estimated{1,var,rep},~,~,~,~,trsimi{1,var,rep},fnorm{1,var,rep}] = ...
            opt_drmave(X, Y, real_beta, m, d, hmave, nmave);

        %============================  NR-B/NR-O ===============================
        fprintf("========================== rep = %d var = %d ======================= \n", rep, variable(var))
        optionsKer = {'mave'};
        % 2. NR-B
        lam = 0;
        T_max= 1000; K_max = 1;N_max=1000;
        [beta_estimated{2,var,rep},trsimi{2,var,rep},objall{2,var,rep},fnorm{2,var,rep}] = ...
            opt_cDRLSNR_block(beta_estimated{1,var,rep},X,Y,m,nn,p,d,T_max,K_max,lam,...
            N_max,real_beta,options,optionsInvVar,optionsKer);
        % 3. NR-O
        lam = 0;
        T_max= 1000; K_max = 1;N_max=2000;
        [beta_estimated{3,var,rep},trsimi{3,var,rep},objall{3,var,rep},fnorm{3,var,rep}]  = ...
            opt_cDRLSNR_orth(beta_estimated{1,var,rep},X,Y,m,nn,p,d,T_max,K_max,lam,...
            N_max,real_beta,options,optionsInvVar,optionsKer);

        %============================  InvexDR ====================================
        optionsKer = {'mave'};
        fprintf("========================== rep = %d var = %d ======================= \n", rep, variable(var))
        % 4. mNRO-IDR
        lam = lamALL;
        T_max= 1000; K_max = 1;N_max=2000;
        [beta_estimated{4,var,rep},trsimi{4,var,rep},objall{4,var,rep},fnorm{4,var,rep}]  = ...
            opt_DRLSNR(beta_estimated{3,var,rep},X,Y,m,nn,p,d,T_max,K_max,lam,...
            N_max,real_beta,options,optionsInvVar,optionsKer);

        % 5. m-IDR
        lam = lamALL;
        T_max= 1000; K_max = 1;N_max=2000;
        [beta_estimated{5,var,rep},trsimi{5,var,rep},objall{5,var,rep},fnorm{5,var,rep}]  = ...
            opt_cDRLSNR_orth(beta_estimated{1,var,rep},X,Y,m,nn,p,d,T_max,K_max,lam,...
            N_max,real_beta,options,optionsInvVar,optionsKer);

        % 6. mNRO-IDR
        % lam = lamALL;
        % T_max= 1000; K_max = 1;N_max=2500;
        % [beta_estimated{6,var,rep},trsimi{6,var,rep},objall{6,var,rep},fnorm{6,var,rep}]  = ...
        %     opt_cDRLSNR_orth(beta_estimated{3,var,rep},X,Y,m,nn,p,d,T_max,K_max,lam,...
        %     N_max,real_beta,options,optionsInvVar,optionsKer);


        %============================  InvexDR-B/O ====================================
        optionsKer = {'mave'};
        fprintf("========================== rep = %d var = %d ======================= \n", rep, variable(var))
        % 7. mNR
        lam = 0;
        T_max= 1000; K_max = 1;N_max=1000;
        [beta_estimated{7,var,rep},trsimi{7,var,rep},objall{7,var,rep},fnorm{7,var,rep}]  = ...
            opt_DRLSNR(beta_estimated{1,var,rep},X,Y,m,nn,p,d,T_max,K_max,lam,...
            N_max,real_beta,options,optionsInvVar,optionsKer);

        % 8. InvexDR-B
        lam = lamALL;
        T_max= 1000; K_max = 1;N_max=2000;
        [beta_estimated{8,var,rep},trsimi{8,var,rep},objall{8,var,rep},fnorm{8,var,rep}]  = ...
            opt_cDRLSNR_block(beta_estimated{7,var,rep},X,Y,m,nn,p,d,T_max,K_max,lam,...
            N_max,real_beta,options,optionsInvVar,optionsKer);
        % 9. InvexDR-O
        lam = lamALL;
        T_max= 1000; K_max = 1;N_max=2000;
        [beta_estimated{9,var,rep},trsimi{9,var,rep},objall{9,var,rep},fnorm{9,var,rep}]  = ...
            opt_cDRLSNR_orth(beta_estimated{7,var,rep},X,Y,m,nn,p,d,T_max,K_max,lam,...
            N_max,real_beta,options,optionsInvVar,optionsKer);


        %============================  pool MAVE  NR-B/O ====================================

        fprintf("========================== rep = %d var = %d ======================= \n", rep, variable(var))
        optionsKer = {'mave'};
        lam = 0;T_max= 1000; K_max = 1;N_max=2000;
        [beta_estimated{10,var,rep},trsimi{10,var,rep},objall{10,var,rep},fnorm{10,var,rep}]  = ...
            opt_cDRLSNR_orth(beta_estimated{7,var,rep},X,Y,m,nn,p,d,T_max,K_max,lam,...
            N_max,real_beta,options,optionsInvVar,optionsKer);

        optionsKer = {'pool'};
        lam = 0;T_max= 1000; K_max = 1;N_max=2000;
        [beta_estimated{11,var,rep},trsimi{11,var,rep},objall{11,var,rep},fnorm{11,var,rep}]  = ...
            opt_poolDRLSNR_block(beta_estimated{1,var,rep},X,Y,X_all_matrix,Y_all_matrix, ...
            m,nn,p,d,T_max,K_max,lam,...
            N_max,real_beta,options,optionsInvVar,optionsKer);

        lam = 0;T_max= 1000; K_max = 1;N_max=2000;
        [beta_estimated{12,var,rep},trsimi{12,var,rep},objall{12,var,rep},fnorm{12,var,rep}]  = ...
            opt_poolDRLSNR_orth(beta_estimated{1,var,rep},X,Y,X_all_matrix,Y_all_matrix, ...
            m,nn,p,d,T_max,K_max,lam,...
            N_max,real_beta,options,optionsInvVar,optionsKer);


    

    end

    if rep ==5, save('SimuV2_varsig_ex1pi3_rep5.mat'), end
    if rep ==10, save('SimuV2_varsig_ex1pi3_rep10.mat'), end

    save('SimuV2_varsig_ex1pi3_rep100.mat')
end

%%
% clear
clc
% load('SimuV2_varn_ex1pi3_rep100.mat')

% 计算tr
repetition = rep; ind_invex = 4;
EvalTr = zeros(12,num_var);EvalF = zeros(12,num_var);
distri_Tr= zeros(12,m);  distri_F = zeros(12,m);
for method = 1:12
    if method == 6, continue; end
    for var = 1:num_var
        tempTr = zeros(1,m);tempF = zeros(1,m);
        for rep = 1:repetition
            resultTr = (trsimi{method,var,rep});
            resultF = (fnorm{method,var,rep});
            r_Tr = resultTr(end,:);r_F = resultF(end,:);
            tempTr = tempTr + r_Tr(end,:);
            tempF = tempF + r_F(end,:);
        end
        if var == 3
            distri_Tr(method,:) = tempTr/repetition; 
            distri_F(method,:)= tempF/repetition;
            dTr = distri_Tr([ind_invex,1,2,3,11,12],:);
            dF  = distri_F([ind_invex,1,2,3,11,12],:);
        end
        EvalTr(method,var) = sum(tempTr)/m/repetition;
        EvalF(method,var) = sum(tempF)/m/repetition;
    end
end



%%
clc
fprintf(' & InvexDR-1 & f  & %5.4f & %5.4f & %5.4f  \\\\\n', EvalF(4,:));
fprintf(' &           & tr & %5.4f & %5.4f & %5.4f  \\\\\n', EvalTr(4,:));

fprintf(' & InvexDR-2 & f  & %5.4f & %5.4f & %5.4f  \\\\\n', EvalF(5,:));
fprintf(' &           & tr & %5.4f & %5.4f & %5.4f  \\\\\n', EvalTr(5,:));

fprintf(' & InvexDR-3 & f  & %5.4f & %5.4f & %5.4f  \\\\\n', EvalF(6,:));
fprintf(' &           & tr & %5.4f & %5.4f & %5.4f  \\\\\n', EvalTr(6,:));


fprintf(' & MAVE      & f  & %5.4f & %5.4f & %5.4f  \\\\\n', EvalF(1,:));
fprintf(' &           & tr & %5.4f & %5.4f & %5.4f  \\\\\n', EvalTr(1,:));

fprintf(' & NR-B      & f  & %5.4f & %5.4f & %5.4f  \\\\\n', EvalF(2,:));
fprintf(' &           & tr & %5.4f & %5.4f & %5.4f  \\\\\n', EvalTr(2,:));
fprintf(' & NR-O      & f  & %5.4f & %5.4f & %5.4f  \\\\\n', EvalF(3,:));
fprintf(' &           & tr & %5.4f & %5.4f & %5.4f  \\\\\n', EvalTr(3,:));

fprintf(' & InvexDR-B & f  & %5.4f & %5.4f & %5.4f  \\\\\n', EvalF(8,:));
fprintf(' &           & tr & %5.4f & %5.4f & %5.4f  \\\\\n', EvalTr(8,:));
fprintf(' & InvexDR-O & f  & %5.4f & %5.4f & %5.4f  \\\\\n', EvalF(9,:));
fprintf(' &           & tr & %5.4f & %5.4f & %5.4f  \\\\\n', EvalTr(9,:));

fprintf(' & NR-O      & f  & %5.4f & %5.4f & %5.4f  \\\\\n', EvalF(10,:));
fprintf(' &           & tr & %5.4f & %5.4f & %5.4f  \\\\\n', EvalTr(9,:));

fprintf(' & pNR-B     & f  & %5.4f & %5.4f & %5.4f  \\\\\n', EvalF(11,:));
fprintf(' &           & tr & %5.4f & %5.4f & %5.4f  \\\\\n', EvalTr(11,:));
fprintf(' & pNR-O     & f  & %5.4f & %5.4f & %5.4f  \\\\\n', EvalF(12,:));
fprintf(' &           & tr & %5.4f & %5.4f & %5.4f  \\\\\n', EvalTr(12,:));