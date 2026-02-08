function [Reg, KerWgt, Hat, Yhat] = KerReg(X,Y,X0,varargin)

% To estimate the conditional mean function E(Y | X). We evaluate at X0.

% INPUT: 
% X: an n x p matrix
% Y: an n x d matrix
% X0: an n0 x p matrix

% OUTPUT: 
% Reg: an n0 x d matrix


% COPYRIGHT
% LIPING ZHU @ NOV 5 2015

% EXAMPLE
% clear all, clc, close all
% n = 100;
% p = 2;
% d = 3;
% X = normrnd(0,1,n,p);
% Y = X(:,1).^2 * ones(1,d) + normrnd(0,1,n,d);
% X0 = X;
% [X0(:,1),IND] = sort(X(:,1));
% X0(:,2:end) = X0(IND,2:end);
% Y0 = X0(:,1).^2;
% [Reg, KerWgt, Hat, Yfit] = KerReg(X,Y,X0);
% figure(1),clf
% plot(X0(:,1),Reg(:,1),'r.'), hold on
% plot(X0(:,1),Reg(:,2),'b.'), hold on
% plot(X0(:,1),Reg(:,3),'g.'), hold on
% plot(X0(:,1), Y0,'-')


[n, p] = size(X);
[m, d] = size(Y);
if n ~= m
    error('Y must be univariate and has the same length as X')
end

[n0, p0] = size(X0);
if n0 == 1 && p0 == 1
    X0 = X0 * ones(1,p);
end
[n0, p0] = size(X0);
if n0 == p0
    if n0 ~= p
        error('Size of X0 does not match X')
    end
else
    if n0 == p
        X0 = X0';
    end
end
n0 = size(X0,1);

% Process additional name/value pair arguments
okargs = {'bandwidth'  'kernel'};
defaults = {[ ]        'epanechinikov'};
[h,kernel] = internal.stats.parseArgs(okargs, defaults, varargin{:});


okkernels = {'normal' 'epanechinikov' 'box' 'triangle'};
if isempty(kernel)
   kernel = okkernels{1};
elseif ~(isa(kernel,'function_handle') || isa(kernel,'inline'))
   if ~ischar(kernel)
      error('Smoothing kernel must be a function.');
   end
   knum = strcmpi(kernel, okkernels);
   if (length(knum) == 1)
      kernel = okkernels{knum};
   end
end

% Default window parameter is optimal for normal distribution
if (isempty(h))
   sig = mad(X,1,1) / 0.6745;
   h = sig * (4/(3*n))^(1/(p+4));
end

h = h(:); l = length(h);   
if l == 1
    h = h*ones(p,1);
end
if length(h) ~= p
    error('Length of h must equal to the number of predictors')
end
KerUni = zeros(n0,n,p);
for dim = 1:p
    Xt = X(:,dim);
    Xt0 = X0(:,dim);
    U = (Xt0 * ones(1,n) - ones(n0,1) * Xt')./h(dim);
    KerUni(:,:,dim) = feval(kernel, U);
end
KerWgt = prod(KerUni,3);
KerSum =  sum(KerWgt,2);
Hat = KerWgt./(KerSum * ones(1,n));                                        % Smoothing matrix, denoted by S(lambda) in some nonparametric context: Yhat = S Y.
Reg = (KerWgt * Y)./(KerSum * ones(1,d));
Yhat = Hat * Y;
% EOF


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                           PLUG-IN FUNCTIONS                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = normal(z)
%NORMAL Normal density kernel.
%f = normpdf(z);
f = exp(-0.5 * z .^2) ./ sqrt(2*pi);


function f = epanechinikov(z)
%EPANECHINIKOV Epanechinikov's asymptotically optimal kernel.
% a = sqrt(5);
% z = max(-a, min(z,a));
% f = .75 * (1 - .2*z.^2) / a;
f = .75 * (1 - z.^2) .* (abs(z) <= 1);


function f = box(z)
%BOX    Box-shaped kernel
a = sqrt(3);
f = (abs(z) <= a) ./ (2 * a);


function f = triangle(z)
%TRIANGLE Triangular kernel.
a = sqrt(6);
z = abs(z);
f = (z<=a) .* (1 - z/a) / a; 




