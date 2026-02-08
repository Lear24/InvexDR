function [X_cell,X_all] = generation_sample(m,n,p,varargin)
% [X_cell,X_all] = generation_sample(m,n,p,varargin)
% Input:
% X = generation_sample(m,n,d,nu)
% m - number of client
% n - sample size (m×1)
% d - dimension of samples
% varargin = 'normal', 't' or 'uniform'
% Output:
% X_cell: m-by-1 cell with cell(m_i) is an p-by-n matrix
% X_all: total nomber of samples across all of nodes


% okkernels = {'normal' 'uniform'};
okargs = {'distribution'  'parameter'};
defaults = {'normal'     0.3   };
[mode,para] = internal.stats.parseArgs(okargs, defaults, varargin{:});


switch mode
    case 'normal'
        [X_cell,X_all] = normal(m,n,p,para);
    case 'uniform'
        [X_cell,X_all] = uniform(m,n,p,para);
end

end

function [X_cell,X_all] = uniform(m,n,p,para)

X_cell = cell(m,1);
X_all = zeros(sum(n),p);
nn = 1;
for i =1:m
    
    X = unifrnd(para(1), para(2), n(i),p);
    X_cell(i,1) = mat2cell(X',p,n(i));
    X_all(nn:nn+n(i)-1,:) = X;
    nn = nn+ n(i);
end
X_all = X_all';

end



function [X_cell,X_all] = normal(m,n,p,para)

X_cell = cell(m,1);
Sigma_X = zeros(p,p);
for i = 1:p
    for j = 1:p
        Sigma_X(i,j) = para^abs(i-j);
    end
end

X_all = zeros(sum(n),p);
nn = 1;
for i =1:m
    
    X = mvnrnd(zeros(p,1),Sigma_X,n(i));
    X_cell(i,1) = mat2cell(X',p,n(i));
    X_all(nn:nn+n(i)-1,:) = X;
    nn = nn+ n(i);
end
X_all = X_all';

end
