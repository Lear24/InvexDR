function P = fun_P(alpha,S,m,varargin)
[p,d] = size(alpha);
mode = 'project'; %默认为向量计算形式
if ~isempty(varargin)                                                      % Check if there are any optional parameters to input
    input_mode = varargin{1};                                              % Extract the first optional parameter
    if ischar(input_mode)                                                  % Verify whether the input is a string
        mode = lower(input_mode);                                          % Convert to lowercase to avoid case sensitivity issues
    else
        error('Optional input must be a string');
    end
end


switch mode                                                                % Select the calculation method based on the mode
    case 'project'
        
        P=0;
        A = util_projection(alpha);
        for i = 1:m
            B = util_projection(cell2mat(S(i,1)));
            P = P + trace(A*B);
        end

        P = (1- P)/2/m/d;                                % Computation in vector mode
    
    
    case 'cosine'
        % in this case S is an pd-by-m matrix including the alpha 
        P = 0;
        for i = 1:m
            
            P = P + util_squredcos(S(:,i),alpha(:));
            
        end
        P = (1-P)/2/m/d;
        
    otherwise
        error('Invalid mode. Use ''project''or ''cosine''');
end










end