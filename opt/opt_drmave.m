function [bb2,bb2_o,bb1,bb1_o, bb0, bb0_o,...
            trsimi,fnorm]= opt_drmave(X, Y, real_beta,m, d, h, nmave,varargin)
% X,Y are m-by-1 cells with p-by-n and n-by-1 matrix, respectively

bb2 = cell(m,1); bb2_o = cell(m,1);
bb1 = cell(m,1); bb1_o = cell(m,1);
bb0 = cell(m,1); bb0_o = cell(m,1);
% if isempty(varargin)
%     mode = "rmave";
% else
%     mode = varargin{1};
% end


for  i = 1:m

    x = X{i,1}'; % n-by-p matrix
    y = Y{i,1};

    disp(i);
    [beta2,beta1,beta0] =  opt_rmave(x, y, d, h, nmave); % x should be a n-by-p matrix


    % switch mode
    %     case "rmave"
    %         beta =  opt_rmave(x, y, d, h, nmave); % x should be a n-by-p matrix
    %     case "mave"
    %         [~,beta,~] =  opt_rmave(x, y, d, h, nmave); % x should be a n-by-p matrix
    % end

    if ~isreal(beta2);   beta2 = real(beta2);   end
    if ~isreal(beta1);   beta2 = real(beta2);   end
    if ~isreal(beta0);   beta2 = real(beta2);   end

    
    bb2{i,1} = beta2; bb2_o{i,1} = orth(beta2);
    bb1{i,1} = beta1; bb1_o{i,1} = orth(beta1);
    bb0{i,1} = beta0; bb0_o{i,1} = orth(beta0);

end

trsimi = util_TraceSimilarity(bb2,real_beta,m,d);
fnorm = util_Fnorm(bb2, real_beta,m,d);

trsimi = trsimi';
fnorm = fnorm';
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