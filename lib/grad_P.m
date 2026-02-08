function dP = grad_P(alpha,S,m,varargin)

[p,d] = size(alpha);
mode = 'project'; %默认为向量计算形式
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
        dP = gradP(alpha,S,m,mode);
    case 'diff'
        dP = diffP(alpha,S,m,mode);
    otherwise
        error('Invalid mode. Use ''grad'' or ''diff''');
end




end



function dP = gradP(alpha,S,m,varargin)
mode = varargin{1};
[p,d] = size(alpha);


switch mode
    case 'project'
        dP = 0;
        A = util_projection(alpha);
        for i = 1:m
            B = util_projection(cell2mat(S(i,1)));
            dP = dP + (eye(p) - A) * B * alpha / (alpha'*alpha);
        end
        dP = -dP/m/d;
        
    case 'cosine'
        dP = 0;
        A = util_projection(alpha(:));
        for i = 1:m
            B = util_projection(S(:,i));
            dP = dP + ( eye(p*d) - A) * B* alpha(:)/(alpha(:)'*alpha(:));
        end
        dP = -dP/m/d;
        dP = reshape(dP,p,d);
    otherwise
        error('Invalid mode. Use ''vector'' or ''matrix''');
end

end









function  dP = diffP(alpha,S,m,varargin)

mode = varargin{1};
[p,d] = size(alpha);
e = 0.0001;

dP=zeros(p,d);
for i=1:p
    for j = 1:d
        I=zeros(p,d);
        I(i,j)=e;
        f1=fun_P(alpha+I,S,m,mode);
        f2=fun_P(alpha-I,S,m,mode);
        dP(i,j)=(f1-f2)/(2*e);
    end
end


end










