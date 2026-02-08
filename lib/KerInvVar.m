function InvVar = KerInvVar(X,Y,beta,m,n,varargin)
% InvVar = KerInvVar(X,Y,beta,m,n,varargin)
% InvVar: m-cell with n-by-1 (w(xi))
% beta: has to be the real beta.


% if iscell(X); modeX = 'cell'; elseif ismatrix(X); modeX = 'matrix'; end
options = varargin{1};
if ~isempty(varargin); mode = options{1};else; mode = 'Kernel';end


switch mode
    case 'Constant'
        w = options{2};
        InvVar = InvVarCons(w,m,n);
    case 'True'
        w = options{2};
        InvVar = InvVarTure(X,beta,m,w);
    case 'Kernel'
        InvVar = InvVarKernel(X,Y,beta,m,n);
end
end

function InvVar = InvVarKernel(X,Y,beta,m,n)
InvVar = cell(m,1);
w = KerWeight(X,Y,beta,m,'nope','bandwidth',[],'kernel','normal');
for i = 1:m
    IV = zeros(n,1);
    ww = w{i}; x=X{i};y=Y{i};b=beta{i};
    [M,~]=KerMeanFuncion(x,y,b,-1,'reduction');
    hy2  = (y - M).^(2);
    for j = 1:n
        IV(j,1) = sum(ww(:,j))/sum(ww(:,j).*hy2);
    end
    InvVar{i} = IV;
end

end

function InvVar = InvVarTure(X,beta,m,varargin)

if ~isempty(varargin); mode = varargin{1}; else; mode = 'setting1'; end
InvVar = cell(m,1);
switch mode
    case 'setting1'
        for i = 1:m
            X1 = X{i};
            xb = X1'* beta{i};
            sigma = (log(10+xb)).^(-1);
            InvVar{i} = sigma;
        end
    case 'setting2'
        for i = 1:m
            X1 = X{i};
            sigma = (exp(X1(1,:)')).^(-1);
            InvVar{i} = sigma;
        end
        
end


end

function InvVar = InvVarCons(w,m,n)

if max(size(w)) == m
    mode = 'matrix';
else
    mode = 'scalar';
end

InvVar = cell(m,1);
switch mode
    case 'matrix'
        for i = 1:m
            InvVar{i} = w(i)*ones(n(i),1);
        end
    case 'scalar'
        for i = 1:m
            InvVar{i} = w*ones(n(i),1);
        end
end

end







