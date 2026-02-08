function [obj,grad] = obj_DHDR(alpha,S,m,n,hatY,hatX,lambda,varargin)

[p,d] = size(alpha);
mode = varargin{1};
modeL = mode{1};
modeP = mode{2};
modegrad = mode{3};




switch modeP 
    case  'cosine'
    S = util_CollVecBeta(S,m,p,d);
end

L = fun_L(alpha,n,hatY,hatX,modeL);
P = fun_P(alpha,S,m,modeP);

dL = grad_L(alpha,n,hatY,hatX,modeL,modegrad);
dP = grad_P(alpha,S,m,modeP,modegrad);

obj = L + lambda*P;
grad = dL + lambda*dP;


end