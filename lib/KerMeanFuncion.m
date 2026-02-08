function [M,M1]=KerMeanFuncion(X,Y,beta,h,varargin)

if iscell(X)
    mode = 'cell';
elseif ismatrix(X)
    mode = 'matrix';
end

if ~isempty(varargin)
    com = varargin{1};
else
    com = 'nope';
end

switch mode
    case 'cell'
        if strcmp(com, 'reduction')
            for i =1:max(size(X))
                X1 = X{i}'* beta{i};
                X{i} = X1';
            end
        end
        [M,M1] = kercell(X,Y,h);
    case 'matrix'
        if strcmp(com, 'reduction')
            X = (X'*beta)';
        end
        [M,M1] = LocalLinear(X',Y,X',h);
end


end



function [M,M1] = kercell(X,Y,h)
m =  max(size(X));
M = cell(m,1);M1 = cell(m,1);
for i = 1:m
    
    X1 = X{i};
    Y1 = Y{i};
    [M{i},M1{i}] = LocalLinear(X1',Y1,X1',h);
    
end

end





