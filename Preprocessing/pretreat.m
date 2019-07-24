function [X,para1,para2]=pretreat(X,method,para1,para2)
%+++   data pretreatment
%+++ HD Li, Central South University


if nargin==2
    if size(X,3)>1
        d=1;
        d1=size(X,1);
        d2=size(X,2);
        X=reshape(X,[],size(X,3));
    else
        d=0;
    end
  [Mx,Nx]=size(X);
   if strcmp(method,'autoscaling')
    para1=mean(X,1);para2=std(X);
   elseif strcmp(method,'center')
    para1=mean(X,1);para2=ones(1,Nx);
   elseif strcmp(method,'unilength')
    para1=mean(X,1);
    for j=1:size(X,2);
    para2(1,j)=norm(X(:,j)-para1(j));
    end
   elseif strcmp(method,'minmax')
    para1=min(X);maxv=max(X);
    para2=maxv-para1;  
   elseif strcmp(method,'pareto');
    para1=mean(X,1);para2=sqrt(std(X));
   else
    display('Wrong data pretreat method!');
   end
   
   for i=1:Nx
     X(:,i)=(X(:,i)-para1(i))/para2(i);
   end
   
   if d==1
      X=reshape(X,d1,d2,size(X,2)); 
   end
   
elseif nargin==4
   [Mx,Nx]=size(X);
   for i=1:Nx     
     X(:,i)=(X(:,i)-para1(i))/para2(i);
   end
end