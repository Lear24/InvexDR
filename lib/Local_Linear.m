function [ffit, fdev] = Local_Linear(y,x,u,varargin)

%   y:  n times 1 response variable
%   x:  n times 1 design points
%   u:  n times 1 grid points


%   ffit = Local_Linear(y,x,u) computes the estimate of regression function
%   evaluated at u. The default estimate is based on a normal kernel
%   function using a window parameter (bandwidth) that is a function of the
%   number of points in univariate x.

%   [ffit, fdev] = Local_Linear(...) also returns the derivative function.

%
%   [ffit, fdev] = Local_Linear(...,'PARAM1',val1,'PARAM2',val2,...) specifies
%   parameter name/value pairs to control the regression estimation.  Valid
%   parameters are the following:
%
%      Parameter    Value
%      'kernel'     The type of kernel smoother to use, chosen from among
%                   'normal' (default), 'box', 'triangle', and 'epanechinikov'.
%      'width'          The bandwidth of the kernel smoothing window.  The default
%                   is optimal for estimating normal densities, but you
%                   may want to choose a smaller value to reveal features
%                   such as multiple modes.
%
%
%   Example:
%         x = randn(30,1); y = exp(x(:,1));
%         f = Local_Linear(y,x,x, 'kernel','epanechinikov','width',1);
%         figure(3),clf,
%         plot(x(:,1),f,'go');
%         hold on,
%         plot(x(:,1),y,'b.');

u = u(:);
x = x(:);
y = y(:);


% % % Process additional name/value pair arguments
% % okargs = {'width'  'kernel'};
% % defaults = {[]  'normal'};
% % [emsg,h,m,kernel] = statgetargs(okargs, defaults, varargin{:});
% % error(emsg);

okargs =   {'width'    'kernel'};
defaults = {[]         'normal'};
[h,kernel] = internal.stats.parseArgs(okargs, defaults, varargin{:});

% [eid,emsg,h,kernel] = ...
%     internal.stats.getargs(okargs, defaults, varargin{:});
% if ~isempty(eid)
%     error(sprintf('stats:locallinear:%s',eid),emsg);
% end

% h=[];
% kernel = [];

% Default window parameter is optimal for normal distribution
if (isempty(h))
   med = median(x);
   n = length(x);
   sig = median(abs(x-med)) / 0.6745;
   h = sig * (4/(3*n))^(1/5);
end


okkernels = {'normal' 'epanechinikov' 'box' 'triangle'};
if isempty(kernel)
   kernel = okkernels{1};
elseif ~(isa(kernel,'function_handle') || isa(kernel,'inline'))
   if ~ischar(kernel)
      error('Smoothing kernel must be a function.');
   end
   knum = strcmpi(kernel, okkernels);
   kernel = okkernels{knum==1};
end


U = (repmat(x,1,length(u)))'-repmat(u,1,length(x)); Ut = U/h;
Y = repmat(y,1,length(u))';
W = feval(kernel, Ut); 

Epsilon = 1E-8;
A = sum(W,2) + Epsilon;
B = sum(W.*U,2);
D = sum(W.*U.*U,2) + Epsilon;

% E = sum(W.*Y,2) + Epsilon;
% F = sum(W.*U.*Y,2);

deno =  A .* D - B.^2;
nume1 = D .* sum(W.*Y,2) - B .* sum(W.*U.*Y,2);
nume2 = -B .* sum(W.*Y,2) + A .* sum(W.*U.*Y,2);


ffit = (deno.^-1) .* nume1;
fdev = (deno.^-1) .* nume2;

% -----------------------------
% The following are functions that define smoothing kernels k(z).
% Each function takes a single input Z and returns the value of
% the smoothing kernel.  These sample kernels are designed to
% produce outputs that are somewhat comparable (differences due
% to shape rather than scale), so they are all probability
% density functions with unit variance.
%
% The density estimate has the form
%    f(x;k,h) = mean over i=1:n of k((x-y(i))/h) / h

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
f = (abs(z)<=a) ./ (2 * a);

function f = triangle(z)
%TRIANGLE Triangular kernel.
a = sqrt(6);
z = abs(z);
f = (z<=a) .* (1 - z/a) / a; 
