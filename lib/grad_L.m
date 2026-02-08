function dL = grad_L(alpha,n,hatY,hatX,varargin)
% mode: vector  matrix
% method: grad  diff

[p,d] = size(alpha);
mode = 'vector'; %默认为向量计算形式
method = 'grad';
if ~isempty(varargin)                                                      % Check if there are any optional parameters to input
    input_mode = varargin{1};
    input_method = varargin{2};                                              % Extract the first optional parameter
    if ischar(input_mode)                                                  % Verify whether the input is a string
        mode = lower(input_mode);                                          % Convert to lowercase to avoid case sensitivity issues
        method = lower(input_method);
    else
        error('Optional input must be a string');
    end
end

switch method
    case 'grad'
        dL = gradL(alpha,n,hatY,hatX,mode);
    case 'diff'
        dL = diffL(alpha,n,hatY,hatX,mode);
    otherwise
        error('Invalid mode. Use ''grad'' or ''diff''');
end

end


function dL = diffL(alpha,n,hatY,hatX,varargin)
mode = varargin{1};
[p,d] = size(alpha);
e = 0.0001;

dL=zeros(p,d);
for i=1:p
    for j = 1:d
        I=zeros(p,d);
        I(i,j)=e;
        f1=fun_L(alpha+I,n,hatY,hatX,mode);
        f2=fun_L(alpha-I,n,hatY,hatX,mode);
        dL(i,j)=(f1-f2)/(2*e); 
    end
end

end



function dL = gradL(alpha,n,hatY,hatX,varargin)
mode = varargin{1};
[p,d] = size(alpha);
switch mode
    case 'vector'
        dL = -2*(hatY - hatX'* alpha(:))'* hatX';
        dL = dL'/n ;
        dL = reshape(dL,p,d);
    case 'matrix'
        dL = 0;
        for i = 1:n
            dL = dL + (-2) * (hatY(i) - trace(hatX(:,(i-1)*d+1:i*d)'*alpha))...
                * hatX(:,(i-1)*d+1:i*d);
        end
        dL = dL/n;
    otherwise
        error('Invalid mode. Use ''vector'' or ''matrix''');
end

end