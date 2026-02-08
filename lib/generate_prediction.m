function [PredError,sumErr,hatY,Y]=generate_prediction(beta,X,Y,m,n,varargin)

if max(size(n)) ==1; n = n*ones(m,1);end

if iscell(beta)
   mode = 'cell';
elseif ismatrix(beta)
    mode = 'matrix';
end

switch mode
    case 'cell'
        hatY=KerMeanFuncion(X,Y,beta,-1,'reduction');
        PredError = zeros(m,1);
        for i = 1:m
           PredError(i) = norm(hatY{i}-Y{i})/n(i);
        end
        sumErr = sum(PredError)/m;
    case 'matrix'
        hatY=KerMeanFuncion(X,Y,beta,-1,'reduction');
        PredError = norm(hatY-Y);
        sumErr = 0;
end

end