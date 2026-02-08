function w = KerWeight(X,Y,beta,m,varargin)

if ~isempty(varargin)
com = varargin{1};
else 
    com = 'nope';
end




if iscell(X)
    mode = 'cell';
elseif ismatrix(X)
    mode = 'matrix';
end

switch mode
    case 'cell'
        if strcmp(com, 'reduction')
            for i =1:m
                X1 = X{i}'* beta{i};
                X{i} = X1';
            end
        end
        w = kercell(X,Y,m,varargin{2:end});
    case 'matrix'
        if strcmp(com, 'reduction')
            X = (X'*beta)';
        end
        w = kermatrix(X,Y,varargin{2:end});
end



end


function w = kermatrix(X,Y,varargin)


    [~, ~, w, ~] = KerReg(X',Y,X',varargin{1:end});
    w = w';


end

function w = kercell(X,Y,m,varargin)

w = cell(m,1);



for i = 1:m
    Y1 = Y{i};
    X1 = X{i};
    [~, ~, w1, ~] = KerReg(X1',Y1,X1',varargin{1:end});
    w1 = w1';
    w{i} = w1;
end







end
