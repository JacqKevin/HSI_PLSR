function [ACP,outliers] = ACPauto(X, fig, l, nbpc)
% PCA that work with corected spectra (centered, normed)
% Find the optimal number of components or the number of component 
% choosen by the user by the algorithm NIPALS
% (Non-linear iterative partial least squares)
%  
%       X: Centered Matrix
%       nbpc: number of components choosen by the user
%       fig: number of the figure if you want to see the PCA
%       l: wavelenght or wavenumber

if nargin<4
    nbpc=0;
end

% X Dimension
[n, m]=size(X);

Xd=X; % Save X

% If there is not wavelenght or wavenumber, creation of a number vector
if nargin<3
    l=1:m; l=l';
end

nspec=1:n; nspec=nspec';

% ============================
%            PCA
% ============================

% Scores and loadings calculate with NIPALS algorithm
R2old=inf;
R2new=0;
i=1;

if nbpc==0 % Optimal number of components
    while abs(R2new-R2old)>0.015
        R2old=R2new;

        [T(:,i), P(:,i), A] = Nipals(X);
        X=X-T(:,i)*P(:,i)'; % X that is not explained by the component i = résidus

    % Calculation of the DModX, residu
            Qres(:,i)=diag(X*X'); % variability on the diagonale of the squared, X'*X: variable, X*X':individu 
            sres(:,i)=sqrt(Qres(:,i)/(m-i)); % Residual standard deviation
            s0=sqrt(sum(Qres(:,i))/(m-i)*(n-i-1));
            DModX(:,i)=sres(:,i)/s0;

    % Calculation of the Hotelling T², distance for the middle of the model
        S2=(n-1)*var(T)/n;
        T2=T.*T./(S2(ones(n,1),:));
        T2Hot=cumsum(T2,2);

    % Calculation of the leverage H
        H=(T2Hot+1)/n;

    % Calculation of the determination coefficient R²
        Qtot=sum(diag(Xd'*Xd)); 
        Qreg=cumsum(diag(T'*T));
        R2=Qreg/Qtot;

    R2new=R2(i,:);

    i=i+1;

    end
    
ncomp=i-1;

else % Number of components choosen by the user
    T=zeros(n,nbpc); P=zeros(m,nbpc);
    Qres=zeros(n,nbpc); sres=zeros(n,nbpc); DModX=zeros(n,nbpc);
    
    for j=1:nbpc

        [Ti, Pi, A] = Nipals(X);
        T(:,j)=Ti; 
        P(:,j)=Pi;
        X=X-T(:,j)*P(:,j)'; % X not explained by the component i = résidus

    % Calculation of the DModX, residu
            Qres(:,j)=diag(X*X'); 
            sres(:,j)=sqrt(Qres(:,j)/(m-j));  
            s0=sqrt(sum(Qres(:,j))/(m-i)*(n-j-1));
            DModX(:,j)=sres(:,j)/s0;

    % Calculation of the Hotelling T², distance from the middle of the model
        S2=(n-1)*var(T)/n;
        T2=T.*T./(S2(ones(n,1),:));
        T2Hot=cumsum(T2,2);

    % Calculation of the leverage H
        H=(T2Hot+1)/n;

    % Calculation of the determination coefficient R²
        Qtot=sum(diag(Xd'*Xd));
        Qreg=cumsum(diag(T'*T));
        R2=Qreg/Qtot;

    end
    
    ncomp=nbpc;
end

% Calculation of the critics values for the DModX and the Hotelling T²
DModX1_crit(:,1:ncomp)=2*std(DModX(:,1:ncomp))+mean(DModX(:,1:ncomp));
T2Hot1_crit(:,1:ncomp)=2*std(T2Hot(:,1:ncomp))+mean(T2Hot(:,1:ncomp));

if ncomp==2||nbpc>0
    a=ncomp;
else
    a=ncomp-1;
end

% Save all the variables in a structure
ACP.T=T(:,1:a);
ACP.P=P(:,1:a);
ACP.T2Hot=T2Hot(:,1:a);
ACP.T2Hotcrit=T2Hot1_crit(:,1:a);
ACP.H=H(:,1:a);
ACP.R2=R2(1:a,:);
ACP.DModX=DModX(:,1:a);
ACP.DModXcrit=DModX1_crit(:,1:a);
ACP.A=A;
ACP.PC=a;

% ============================
% 	Outliers tests
% ============================

[outliers] = Outliers(ncomp,DModX,T2Hot,DModX1_crit,T2Hot1_crit);

% ============================
% 	Graphical figure
% ============================
if nargin>1&&fig>0
% 
% figure(fig+3) % Affichage des spectres
% plot(l, Xd(:,:))
% grid on
% xlabel('Longueur d''onde')
% title('Spectres')
    
% DModX
figure(fig+2)
subplot(2,1,1)
plot(nspec, DModX(:,ncomp),'b')
hold on
line([1 n],[DModX1_crit(:,ncomp) DModX1_crit(:,ncomp)],'Color','r')
grid on
xlabel ('Numéro spectre')
ylabel ('DModX')
title('Graphique des DModX')

% Hotelling T²
subplot(2,1,2)
plot(nspec, T2Hot(:,ncomp),'b')
hold on
line([1 n],[T2Hot1_crit(:,ncomp) T2Hot1_crit(:,ncomp)],'Color','r')
grid on
xlabel ('Numéro spectre')
ylabel ('T² d''Hotelling')
title('T² d''Hotelling')
   
% Scores
figure(fig+1)
subplot(2,1,1) 
plot(T(:,1),T(:,2),'b.','Markersize',15)
%hold on
%ellipse(0,0,2.8*std(sqrt(T2Hot(:,1))),2.8*std(sqrt(T2Hot(:,2)))) % A VERIFIER
grid on
xlabel (strcat('PC1 : ',num2str(R2(1,:)*100,3),' %'))
ylabel (strcat('PC2 : ',num2str((R2(2,:)-R2(1,:))*100,3),' %'))
title('Graphique des scores')

% Loadings
subplot(2,1,2)
plot(l, P(:,1)/std(P(:,1)),'b-')
hold on
plot(l, P(:,2)/std(P(:,2)), 'r-')
grid on
xlabel('Longueur d''onde')
ylabel('P[1],P[2]')
legend('P[1]','P[2]')
title('Graphique des loadings')

% R²
figure(fig)
bar(R2)
grid on
xlabel('PC')
ylabel('R²')
title('Graphique des R²')

end

end