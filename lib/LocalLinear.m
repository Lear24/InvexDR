function [f,dev,SE,S,h]=LocalLinear(x,y,x0,h0)

% input
% x: the observed predictors [n,p]
% y: the observed response [n,1]
% x0: the grid points [nn,p]
% h0: the bandwidth

% output
% f: the fitted value [nn,1]
% S: the variance of the error term [nn,1]
% SE: the standard deviation of the estimated regression function [nn,1]
% dev: the derivative of the regression function [nn,p]
% h: the bandwidth selected by GCV

    [ss,ncov]=size(x);                                                  %obtain the sample size and data dimension
    D=diag(ones(1+ncov,1))*1.0E-9;                                      %the ridge matrix
    stde=sqrt(mean((y-x*inv(x'*x)*(x'*y)).^2));                         %a rough estimated STDE
    
    if h0<0
        ngrids=11;                                                      %number of grid points
        cv=zeros(ngrids,2);                                             %record the pred. err.
        for k=1:ngrids                                                  %do for each grid point
            h=-5+(k-1)/(ngrids-1)*10;h=exp(h)*stde*ss^(-0.2);           %different bandwidth
            f=zeros(ss,1);                                              %record f estimates
            for i=1:ss                                                  %update f est. at each obs.
                xx=x-ones(ss,1)*x(i,:);                                 %take difference in X-variable
                ww=sqrt((1-sum(xx.^2,2)/(h*h)).*(1>sum(xx.^2,2)/(h*h)));%Epanechnikov kernel
                ww(i,:)=0;                                              %leaving the i'th obs out
                X=[ones(ss,1),xx];                                      %the expanded X matrix
                X=X.*repmat(ww,1,1+ncov);                               %create weighted predictors
                Y=y.*ww;                                                %create weighted response
                b=(Y'*X/ss)*inv(X'*X/ss+D);                             %update the local estimate
                f(i,:)=b(1);                                            %get function estimate
            end
            cv(k,1)=h;                                                  %record the bandwidth
            cv(k,2)=mean((y-f).^2);                                     %record pred. error.
        end
        pos=find(cv(:,2)==min(cv(:,2)));                                %find the minimal cvndwidth
        pos=pos(1);                                                     %just in case there is a tie
        h=cv(pos,1);                                                    %the optimal bandwidth
        flag=1*(pos>1)*(pos<ngrids);                                    %whether on the boundary
    else
        h=h0;flag=1;cv=0;                                               %if user specified the bandwidth
    end
     
    ss0=length(x0);                                                     %the number of gridpoints
    f=zeros(ss0,1);dev=zeros(ss0,ncov);                                 %the output B  matrix
    SE=zeros(ss0,1);                                                    %SE of the est. of reg. func.
    SD2=zeros(ss0,1); S = SD2;                                          %the estimated resid. var.
    for i=1:ss0                                                         %update reg est. at each obs.
        xx=x-ones(ss,1)*x0(i,:);                                        %take difference in X-variable
        ww=sqrt((1-sum(xx.^2,2)/(h*h)).*(1>sum(xx.^2,2)/(h*h)));        %Epanechnikov kernel
        xx=[ones(ss,1),xx];                                             %the expanded X matrix
        X=xx.*repmat(ww,1,1+ncov);                                      %create weighted predictors
        X2=xx.*repmat(ww.^2,1,1+ncov);                                  %create weighted predictors
        Y=y.*ww;                                                        %create weighted response
        b=(Y'*X/ss)*inv(X'*X/ss+D);                                     %update the local estimate
        f(i,:)=b(1);                                                    %get reg. func. estimate
        dev(i,:)=b([2:ncov+1]);                                         %get derivative
        wresid=Y-X*b';                                                  %weighted resid: X&W are weighted pred.
        SD2(i)=mean(wresid.^2);                                         %the density weighted cond. var.
        S(i)=SD2(i)/mean(ww.^2);                                        %the local var.
        s=sqrt(S(i));                                                   %the local standard error
        COV=inv(X'*X/ss+D)*(X2'*X2/ss)*inv(X'*X/ss+D)/ss;               %the sandwitch COV matrix
        se=sqrt(diag(COV));SE(i,:)=se(1)*s;                             %get the reg. func. SE. est.
    end
    SD2=mean(SD2);                                                      %the final SD2: key input for BIC 
    