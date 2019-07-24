function [B,Wstar,T,P,Q,W,R2X,R2Y, RMSE]=plsnipals(X,Y, vl) %
% The NIPALS (Non-linear iterative partial least squares) algorithm 
% for PLS with one or more Y
% Inputs Arguments:
%       X: n x p matrix
%       Y: n x m matrix
%       vl: number of latent variables

% If the user do not choose the number of latent variables
if nargin<3
    vl=0;
end

% Total variability
varX=sum(sum(X.^2));
varY=sum(sum(Y.^2));

% PLS with 15 latents variables
for i=1:15
    error=1;
    u=Y(:,1);
    niter=0;
    while (error>1e-8 && niter<1000)  % Convergence test
        w=X'*u/(u'*u);
        w=w/norm(w);
        t=X*w;
        q=Y'*t/(t'*t);  % regress Y against t;
        u1=Y*q/(q'*q); % regress q against Y
        error=norm(u1-u)/norm(u);
        u=u1;
        niter=niter+1; % number of iteration
    end
    p=X'*t/(t'*t);
    X=X-t*p';
    Y=Y-t*q';
    
    % Store the model parameters
    W(:,i)=w;
    T(:,i)=t;
    P(:,i)=p;
    Q(:,i)=q;
    
    % Calculate explained variance
    R2X=diag(T'*T*P'*P)/varX;
    R2Y=diag(T'*T*Q'*Q)/varY;    
    
    Wstar1=W*(P'*W)^(-1); 
    % RMSEC
    Y0 = X*Wstar1*Q';
    [ac,~]=size(Y0);
    RMSE(i,:)=sqrt(sum((Y0-Y).^2)/(ac-i-1));

end

if vl<2 % If the user do not choose the number of latents variables
    R2oldX=Inf;
    R2newX=0;
    R2oldY=Inf;
    R2newY=0;
    i=1;
    % Optimal number of latents variables
    while abs(R2newX-R2oldX)>0.02||abs(R2newY-R2oldY)>0.02&&i>2&&i<16
        R2oldX=R2newX;
        R2oldY=R2newY;
        
        R2newX=R2X(i,:);
        R2newY=R2Y(i,:);

        i=i+1;
    end
    vl=i-1;
    % Store the optimal number of latents variables
    W=W(:,1:vl);
    T=T(:,1:vl);
    P=P(:,1:vl);
    Q=Q(:,1:vl);
else % If the user choose the number of latents variables
    W=W(:,1:vl);
    T=T(:,1:vl);
    P=P(:,1:vl);
    Q=Q(:,1:vl);
end

Wstar=W*(P'*W)^(-1); 
B=Wstar*Q'; % Calculate the regressions coefficients
Q=Q';