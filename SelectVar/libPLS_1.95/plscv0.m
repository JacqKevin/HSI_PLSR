function CV=plscv(X,y,A,K,method,PROCESS,order)
%+++ K-fold Cross-validation for PLS
%+++ Input:  X: m x n  (Sample matrix)
%            y: m x 1  (measured property)
%            A: The maximal number of latent variables for cross-validation
%            K: fold. when K=m, it is leave-one-out CV
%       method: pretreatment method. Contains: autoscaling,
%               pareto,minmax,center or none.
%      PROCESS: =1 : print process.
%               =0 : don't print process.
%+++ Order: =0  sorted, default. For CV partition.
%           =1  random.
%           =2  original.
%+++ Output: Structural data: CV
%+++ Hongdong Li, Oct. 16, 2008.
%+++ Tutor: Prof. Yizeng Liang, yizeng_liang@263.net.
%+++ Contact: lhdcsu@gmail.com.

if nargin<7;order=0;end
if nargin<6;PROCESS=1;end
if nargin<5;method='center';end
if nargin<4;K=10;end
if nargin<3;A=3;end



if order==0
    [ydx,indexyy]=sort(sum(y,2));
    X=X(indexyy,:);
    y=y(indexyy,:);
elseif order==1
    indexyy=randperm(length(y));
    X=X(indexyy,:);
    y=y(indexyy,:);
elseif order==2
    indexyy=1:length(y);
    X=X(indexyy,:);
    y=y(indexyy,:);
end

[Mx,Nx]=size(X);
A=min([size(X) A]);
yytest=nan(Mx,size(y,2));
YR=nan(Mx,size(y,2));

groups = 1+rem(0:Mx-1,K);
for group=1:K
    
    calk = find(groups~=group);
    testk = find(groups==group);
    
    Xcal=X(calk,:);
    ycal=y(calk,:);
    
    Xtest=X(testk,:);
    ytest=y(testk,:);
    
    %   data pretreatment
    [Xs,xpara1,xpara2]=pretreat(Xcal,method);
    [ys,ypara1,ypara2]=pretreat(ycal,'center');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [B,Wstar,T,P,Q]=plsnipals(Xs,ys,A);   % no pretreatment.
    
    yp=[];
    
    %+++ calculate the coefficient linking Xcal and ycal.
    C=repmat(ypara2,size(B,1),1).*B./repmat(xpara2,size(B,2),1)';
    coef=[C;ypara1-xpara1*C;];
    
    %+++ predict
    Xteste=[Xtest ones(size(Xtest,1),1)];
    ypred=Xteste*coef;
    yp=[yp ypred];
    
    YR(testk,:)=yp;
    yytest(testk,:)=ytest;
    if PROCESS==1; 
        %fprintf('The %dth group finished.\n',group); 
    end
end

%+++ return the original order
YR(indexyy,:)=YR;
y(indexyy,:)=y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% mean and sd of squared error %%%%%%%%%%%%
error=YR-y;
error2=error.^2;
error2_MEAN=nansum(error2)/Mx;
error2_SD= sqrt(nansum((error2-repmat(nanmean(error2),Mx,1)).^2)/(Mx-1)); % unbiased estimator

%+++ calculate Q2
cv=sqrt(error2_MEAN);
RMSEP=min(cv);
SST=sumsqr(yytest-mean(y));
SSE=sumsqr(YR-y);
Q2=1-SSE/SST;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%+++ output  %%%%%%%%%%%%%%%%
CV.method=method;
CV.Ypred=YR;
CV.predError=error;
CV.RMSECV=cv;
CV.Q2=Q2;
CV.RMSECV_min=RMSEP;
CV.Q2_max=Q2;
CV.optLV=size(T,2);
CV.note='*** The following is based on global min MSE + 1SD';
CV.RMSECV_min_1SD=cv;
CV.Q2_max_1SD=Q2;
CV.optLV_1SD=size(T,2);
%%%%%%%%%%%%%%%%%%%%%%%%%s%%%%