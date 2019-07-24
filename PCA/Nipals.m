function [t, p, i] = Nipals(X)
% Function NIPALS (Non-linear iterative partial least squares) 
% to calculate components for PCA.

    % Variable of high variability
    [~, j]=find(var(X)==max(var(X)));
    
    % Initialisation of t:
    t=X(:,j(1));
    oldlambda2=inf; % define the value n-1 at infinity
    lambda2=t'*t; % define the value n
    i=0;
    
    % NIPALS:
    while abs(lambda2-oldlambda2) > eps&&i<5000 % On cherche la convergence de t, que l'écart soit proche de zéro (epsilon)
        oldlambda2=lambda2; % on défini la valeur n-1
        p=X'*t/norm(X'*t);
        t=X*p;
        lambda2=t'*t; % on défini la valeur n
        i=i+1;
    end
    
end