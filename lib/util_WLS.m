function [a,b,beta] = util_WLS(x,y,w,varargin)
% x: p-by-n
% y: n-by-1
% w: n-by-1
[p,n] = size(x);
mode = varargin{1};
if strcmp(mode,'intercept')
    X = [x; ones(1, n) ];
end

W = diag(w);

beta = (X * W * X')^(-1) * (X * W * y);

a = beta(p+1);
b = beta(1:p);
end 