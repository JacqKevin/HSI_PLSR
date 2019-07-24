function Lvip=VIP(X,Y)
%+++ Calculate the vip for each variable to the response;
%+++ T,W,Q are output from plsnipals.m
%+++ T: socres, which can be obtained by pls_nipals.m 
%+++ W: weight, which can be obtained by pls_nipals.m 
%+++ Q: Y-loadings
%+++ VIP=sqrt(p*q/s);
%+++      where, p is a constant, denoting the number of variables in x
%                q stands for the explained variance of Y by each variable
%                s represents the total variance explained by A components
%+++ Reference: Tahir Mehmood et al, Chemometrics and Intelligent Laboratory Systems 118 (2012)62?69
%+++ HDLi

[B,Wstar,T,P,Q,W,R2X,R2Y, RMSE]=plsnipals(Centrer(X),Y);

s=diag(T'*T*Q*Q');
%initializing
[~,p]=size(X);
[~,h]=size(T);
%+++ calculate VIP;
Lvip=[];
for i=1:p
    weight=[];
    for j=1:h
        weight(j,1)= (W(i,j)/norm(W(:,j)))^2;
    end
    q=s'*weight;  % explained variance by variable i
    Lvip(i)=sqrt(p*q/sum(s));   
end