function L = fun_L(alpha,n,hatY,hatX,varargin)
[p,d] = size(alpha);

mode = 'vector'; %默认为向量计算形式


if ~isempty(varargin)                                                      % Check if there are any optional parameters to input
    
    input_mode = varargin{1};                                              % Extract the first optional parameter
    
    
    if ischar(input_mode)                                                  % Verify whether the input is a string
        
        mode = lower(input_mode);                                          % Convert to lowercase to avoid case sensitivity issues
        
    else
        error('Optional input must be a string');
    end
end


switch mode                                                                % Select the calculation method based on the mode
    case 'vector'
        L = norm(hatY - hatX'* alpha(:)).^2/n;                                  % Computation in vector mode
    case 'matrix'
        L = 0;
        for i = 1:n
            L = L + (hatY(i) - trace(hatX(:,(i-1)*d+1:i*d)'*alpha)).^2/n;
        end
    otherwise
        error('Invalid mode. Use ''vector'' or ''matrix''');
end



end