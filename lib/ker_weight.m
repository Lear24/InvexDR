function w = ker_weight(X,beta,h,m,n,p,d,varargin)

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
        w = kercell(X,h,m,n,p,d);
    case 'matrix'
        if strcmp(com, 'reduction')
            X = (X'*beta)';
        end
        w = kermatrix(X,h,m,n,p,d);
end



end


function w = kermatrix(X,h,m,n,p,d)


    w = zeros(n,n);
    for j = 1:size(X,2)
    x = X(:,j);
    diff = X - repmat(x,1,n);
    diagd = diag(diff'*diff);
    
    ww = ((1/sqrt(2*pi))^n)*(1/h^n)*exp((-1/2)*diagd/h^2);
    w(:,j) = ww/sum(ww);
    
    end
    


end

function w = kercell(X,h,m,n,p,d)

w = cell(m,1);



for i = 1:m
    
    X1 = X{i};
    w1 = zeros(n(i),n(i));
    for j = 1:size(X1,2)
    x = X1(:,j);
    diff = X1 - repmat(x,1,n(i));
    diagd = diag(diff'*diff);
    
    ww = ((1/sqrt(2*pi))^n(i))*(1/h^n(i))*exp((-1/2)*diagd/h^2);
    w1(:,j) = ww/sum(ww);
    
    end
    
    w{i} = w1;
end







end
