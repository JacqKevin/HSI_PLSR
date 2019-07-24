function pred=PLScv(S,Y,Num,wl,cv,pretc,selectvar)
% Function that allows to create PLSR with a unique dataset and several
% preprocessing techniques.

% INPUT:
%       S: Spectra dataset
%       Y: Reference dataset
%       Num: If there is several reference from the same sample
%       wm: Wavelengths
%       cv: Cross-validation vector [mod nb nbt]
%               mod:    1:'KFold' (default)
%                       2:'HoldOut'
%                       validation
%                       3:'LeaveOut'
%               nb: Percentage or ratio for the validation set
%               nbt: Number of iteration
%       pretc: Preprocessing choice (0: All (default), 1:detrend,
%           2:SNV, 3:SNVD, 4:MSC, 5:D1, 6:D2, 7:SNV+D1, 8:SNVD+D1,
%           9:MSC+D1, 10:SNV+D2, 11:SNVD+D2, 12:MSC+D2, 13:CR)
%       selectvar: Wavelengths selected (by default all are used)

% Initialisation:
ptd1=7;
ptd2=9;
polyd1=2;
polyd2=2;
if nargin<5
    cv=[1 0.50 100];
end
contrainte=0;
fig=0;
if nargin<6
    pretc=0;
end
if nargin<7
    selectvar=1:length(wl);
end

% Set definition
[ical,ival]=SetDefinition(size(S,1),round(1/4*size(S,1)));
Scal=S(ical,:);
Sval=S(ival,:);
Ycal=Y(ical,:);
Yval=Y(ival,:);
Numcal=Num(ical,:);

%============
% Preprocessing
%============
if pretc==1
    Xpretcal = Scal;
    Xpretval = Sval;
else
    if pretc>1
        Xpretcal = AllPret(Scal,wl,pretc-1);
        Xpretval = AllPret(Sval,wl,pretc-1);
    else
        Xpretcal = AllPret(Scal,wl);
        Xpretval = AllPret(Sval,wl);
    end
    
end

pred.Scal=Xpretcal;
pred.Sval=Xpretval;
pred.Ycal=Ycal;
pred.Yval=Yval;

%============
% RAW
%============
if pretc==0||pretc==1
    predbrut=predPLS(Scal(:,selectvar), Sval(:,selectvar), Ycal, Yval, 0, fig,contrainte);
end
%============
% DETREND
%============
if pretc==0
    preddetrend=predPLS(Xpretcal.Detrend(:,selectvar),Xpretval.Detrend(:,selectvar),Ycal,Yval,0,fig,contrainte);
end
if pretc==2
    preddetrend=predPLS(Xpretcal(:,selectvar),Xpretval(:,selectvar),Ycal,Yval,0,fig,contrainte);
end
%============
% SNV
%============
if pretc==0
    predsnv=predPLS(Xpretcal.SNV(:,selectvar), Xpretval.SNV(:,selectvar), Ycal, Yval, 0, fig,contrainte);
end
if pretc==3
    predsnv=predPLS(Xpretcal(:,selectvar), Xpretval(:,selectvar), Ycal, Yval, 0, fig,contrainte);
end
%============
% MSC
%============
if pretc==0
    predmsc=predPLS(Xpretcal.MSC(:,selectvar), Xpretval.MSC(:,selectvar), Ycal, Yval, 0, fig,contrainte);
end
if pretc==4
    predmsc=predPLS(Xpretcal(:,selectvar), Xpretval(:,selectvar), Ycal, Yval, 0, fig,contrainte);
end
%============
% SNV + detrend
%============
if pretc==0
    predsnvd=predPLS(Xpretcal.SNVD(:,selectvar), Xpretval.SNVD(:,selectvar), Ycal, Yval, 0, fig,contrainte);
end
if pretc==5
    predsnvd=predPLS(Xpretcal(:,selectvar), Xpretval(:,selectvar), Ycal, Yval, 0, fig,contrainte);
end
%============
% D1
%============
if pretc==0
    predd1=predPLS(Xpretcal.D1(:,selectvar), Xpretval.D1(:,selectvar), Ycal, Yval, 0, fig,contrainte);
end
if pretc==6
    predd1=predPLS(Xpretcal(:,selectvar), Xpretval(:,selectvar), Ycal, Yval, 0, fig,contrainte);
end
%============
% D2
%============
if pretc==0
    predd2=predPLS(Xpretcal.D2(:,selectvar), Xpretval.D2(:,selectvar), Ycal, Yval, 0, fig,contrainte);
end
if pretc==7
    predd2=predPLS(Xpretcal(:,selectvar), Xpretval(:,selectvar), Ycal, Yval, 0, fig,contrainte);
end
%============
% MSC + D1
%============
if pretc==0
    predmscd1=predPLS(Xpretcal.MSCD1(:,selectvar), Xpretval.MSCD1(:,selectvar), Ycal, Yval, 0, fig,contrainte);
end
if pretc==8
    predmscd1=predPLS(Xpretcal(:,selectvar), Xpretval(:,selectvar), Ycal, Yval, 0, fig,contrainte);
end
%============
% MSC + D2
%============
if pretc==0
    predmscd2=predPLS(Xpretcal.MSCD2(:,selectvar), Xpretval.MSCD2(:,selectvar), Ycal, Yval, 0, fig,contrainte);
end
if pretc==9
    predmscd2=predPLS(Xpretcal(:,selectvar), Xpretval(:,selectvar), Ycal, Yval, 0, fig,contrainte);
end
%============
% SNV + D1
%============
if pretc==0
    predsnvd1=predPLS(Xpretcal.SNVD1(:,selectvar), Xpretval.SNVD1(:,selectvar), Ycal, Yval, 0, fig,contrainte);
end
if pretc==10
    predsnvd1=predPLS(Xpretcal(:,selectvar), Xpretval(:,selectvar), Ycal, Yval, 0, fig,contrainte);
end
%============
% SNV + D2
%============
if pretc==0
    predsnvd2=predPLS(Xpretcal.SNVD2(:,selectvar), Xpretval.SNVD2(:,selectvar), Ycal, Yval, 0, fig,contrainte);
end
if pretc==11
    predsnvd2=predPLS(Xpretcal(:,selectvar), Xpretval(:,selectvar), Ycal, Yval, 0, fig,contrainte);
end
%============
% SNV + detrend + D1
%============
if pretc==0
    predsnvdd1=predPLS(Xpretcal.SNVDD1(:,selectvar), Xpretval.SNVDD1(:,selectvar), Ycal, Yval, 0, fig,contrainte);
end
if pretc==12
    predsnvdd1=predPLS(Xpretcal(:,selectvar), Xpretval(:,selectvar), Ycal, Yval, 0, fig,contrainte);
end
%============
% SNV + detrend + D2
%============
if pretc==0
    predsnvdd2=predPLS(Xpretcal.SNVDD2(:,selectvar), Xpretval.SNVDD2(:,selectvar), Ycal, Yval, 0, fig,contrainte);
end
if pretc==13
    predsnvdd2=predPLS(Xpretcal(:,selectvar), Xpretval(:,selectvar), Ycal, Yval, 0, fig,contrainte);
end
%============
% CR
%============
if pretc==0
    predcr=predPLS(Xpretcal.CR(:,selectvar), Xpretval.CR(:,selectvar), Ycal, Yval, 0, fig,contrainte);
end
if pretc==14
    predcr=predPLS(Xpretcal(:,selectvar), Xpretval(:,selectvar), Ycal, Yval, 0, fig,contrainte);
end
%===================
%============
% CROSS VALIDATION
%============
%===================

[C,ia,ic] = unique(Numcal);
[nu, ~]=size(C);
[ny, ~]=size(Ycal);
if ny~=nu
    n=nu;
    
    [Mcv1]= CrossValDef(Scal(ia,:),cv(1,1),cv(1,2),cv(1,3));
    Mcv=zeros(ny,cv(1,3));
    
    for j=1:nu
        for k=1:cv(1,3)
            if Mcv1(j,k)==1
                [ni,~]=find(ic==j);
                Mcv(ni,k)=1;
            end
        end
    end
    
else
    [Mcv]= CrossValDef(Scal,cv(1,1),cv(1,2),cv(1,3));
end

[~, mcv]=size(Mcv);
if pretc==0||pretc==1
    brutSEP=zeros(mcv,size(Ycal,2));brutr2p=zeros(mcv,size(Ycal,2));brutrpdcv=zeros(mcv,size(Ycal,2));
end
if pretc==0||pretc==2
    detrendSEP=zeros(mcv,size(Ycal,2));detrendr2p=zeros(mcv,size(Ycal,2));detrendrpdcv=zeros(mcv,size(Ycal,2));
end
if pretc==0||pretc==3
    snvcvSEP=zeros(mcv,size(Ycal,2));snvcvr2p=zeros(mcv,size(Ycal,2));snvcvrpd=zeros(mcv,size(Ycal,2));
end
if pretc==0||pretc==4
    msccvSEP=zeros(mcv,size(Ycal,2));msccvr2p=zeros(mcv,size(Ycal,2));msccvrpd=zeros(mcv,size(Ycal,2));
end
if pretc==0||pretc==5
    snvdSEP=zeros(mcv,size(Ycal,2));snvdr2p=zeros(mcv,size(Ycal,2));snvdrpdcv=zeros(mcv,size(Ycal,2));
end
if pretc==0||pretc==6
    d1SEP=zeros(mcv,size(Ycal,2));d1r2p=zeros(mcv,size(Ycal,2));d1rpdcv=zeros(mcv,size(Ycal,2));
end
if pretc==0||pretc==7
    d2SEP=zeros(mcv,size(Ycal,2));d2r2p=zeros(mcv,size(Ycal,2));d2rpdcv=zeros(mcv,size(Ycal,2));
end
if pretc==0||pretc==8
    mscd1SEP=zeros(mcv,size(Ycal,2));mscd1r2p=zeros(mcv,size(Ycal,2));mscd1rpdcv=zeros(mcv,size(Ycal,2));
end
if pretc==0||pretc==9
    mscd2SEP=zeros(mcv,size(Ycal,2));mscd2r2p=zeros(mcv,size(Ycal,2));mscd2rpdcv=zeros(mcv,size(Ycal,2));
end
if pretc==0||pretc==10
    snvd1SEP=zeros(mcv,size(Ycal,2));snvd1r2p=zeros(mcv,size(Ycal,2));snvd1rpdcv=zeros(mcv,size(Ycal,2));
end
if pretc==0||pretc==11
    snvd2SEP=zeros(mcv,size(Ycal,2));snvd2r2p=zeros(mcv,size(Ycal,2));snvd2rpdcv=zeros(mcv,size(Ycal,2));
end
if pretc==0||pretc==12
    snvdd1SEP=zeros(mcv,size(Ycal,2));snvdd1r2p=zeros(mcv,size(Ycal,2));snvdd1rpdcv=zeros(mcv,size(Ycal,2));
end
if pretc==0||pretc==13
    snvdd2SEP=zeros(mcv,size(Ycal,2));snvdd2r2p=zeros(mcv,size(Ycal,2));snvdd2rpdcv=zeros(mcv,size(Ycal,2));
end
if pretc==0||pretc==14
    crSEP=zeros(mcv,size(Ycal,2));crr2p=zeros(mcv,size(Ycal,2));crrpdcv=zeros(mcv,size(Ycal,2));
end

if pretc==0
    parfor i=1:mcv
        setcv=Mcv(:,i);
        Scvcal=Scal(setcv==1,:);
        Scvval=Scal(setcv==0,:);
        Ycvcal=Ycal(setcv==1,:);
        Ycvval=Ycal(setcv==0,:);
        
        %============
        % PREPROCESSING
        %============
        Xpretcvcal = AllPret(Scvcal,wl);
        Xpretcvval = AllPret(Scvval,wl);
        
        %============
        % RAW
        %============
        predcv=predPLS(Scvcal(:,selectvar), Scvval(:,selectvar), Ycvcal, Ycvval, 0, 0, contrainte,predbrut.para.nbcomp);
        brutSEP(i,:)=predcv.coeff.SEC;
        brutr2p(i,:)=predcv.coeff.R2c;
        brutrpdcv(i,:)=predcv.coeff.RPDc;
        %============
        % DETREND
        %============
        predcv=predPLS(Xpretcvcal.Detrend(:,selectvar),Xpretcvval.Detrend(:,selectvar),Ycvcal,Ycvval,0,0,contrainte,preddetrend.para.nbcomp);
        detrendSEP(i,:)=predcv.coeff.SEC;
        detrendr2p(i,:)=predcv.coeff.R2c;
        detrendrpdcv(i,:)=predcv.coeff.RPDc;
        %============
        % SNV
        %============
        predcv=predPLS(Xpretcvcal.SNV(:,selectvar), Xpretcvval.SNV(:,selectvar), Ycvcal, Ycvval, 0, 0, contrainte,predsnv.para.nbcomp);
        snvcvSEP(i,:)=predcv.coeff.SEC;
        snvcvr2p(i,:)=predcv.coeff.R2c;
        snvcvrpd(i,:)=predcv.coeff.RPDc;
        %============
        % MSC
        %============
        predcv=predPLS(Xpretcvcal.MSC(:,selectvar), Xpretcvval.MSC(:,selectvar), Ycvcal, Ycvval, 0, 0, contrainte,predmsc.para.nbcomp);
        msccvSEP(i,:)=predcv.coeff.SEC;
        msccvr2p(i,:)=predcv.coeff.R2c;
        msccvrpd(i,:)=predcv.coeff.RPDc;
        %============
        % SNV + detrend
        %============
        predcv=predPLS(Xpretcvcal.SNVD(:,selectvar), Xpretcvval.SNVD(:,selectvar), Ycvcal, Ycvval, 0, 0, contrainte,predsnvd.para.nbcomp);
        snvdSEP(i,:)=predcv.coeff.SEC;
        snvdr2p(i,:)=predcv.coeff.R2c;
        snvdrpdcv(i,:)=predcv.coeff.RPDc;
        %============
        % D1
        %============
        predcv=predPLS(Xpretcvcal.D1(:,selectvar), Xpretcvval.D1(:,selectvar), Ycvcal, Ycvval, 0, 0, contrainte, predd1.para.nbcomp);
        d1SEP(i,:)=predcv.coeff.SEC;
        d1r2p(i,:)=predcv.coeff.R2c;
        d1rpdcv(i,:)=predcv.coeff.RPDc;
        %============
        % D2
        %============
        predcv=predPLS(Xpretcvcal.D2(:,selectvar), Xpretcvval.D2(:,selectvar), Ycvcal, Ycvval, 0, 0, contrainte, predd2.para.nbcomp);
        d2SEP(i,:)=predcv.coeff.SEC;
        d2r2p(i,:)=predcv.coeff.R2c;
        d2rpdcv(i,:)=predcv.coeff.RPDc;
        %============
        % MSC + D1
        %============
        predcv=predPLS(Xpretcvcal.MSCD1(:,selectvar), Xpretcvval.MSCD1(:,selectvar), Ycvcal, Ycvval, 0, 0, contrainte, predmscd1.para.nbcomp);
        mscd1SEP(i,:)=predcv.coeff.SEC;
        mscd1r2p(i,:)=predcv.coeff.R2c;
        mscd1rpdcv(i,:)=predcv.coeff.RPDc;
        %============
        % MSC + D2
        %============
        predcv=predPLS(Xpretcvcal.MSCD2(:,selectvar), Xpretcvval.MSCD2(:,selectvar), Ycvcal, Ycvval, 0, 0, contrainte, predmscd2.para.nbcomp);
        mscd2SEP(i,:)=predcv.coeff.SEC;
        mscd2r2p(i,:)=predcv.coeff.R2c;
        mscd2rpdcv(i,:)=predcv.coeff.RPDc;
        %============
        % SNV + D1
        %============
        predcv=predPLS(Xpretcvcal.SNVD1(:,selectvar), Xpretcvval.SNVD1(:,selectvar), Ycvcal, Ycvval, 0, 0, contrainte, predsnvd1.para.nbcomp);
        snvd1SEP(i,:)=predcv.coeff.SEC;
        snvd1r2p(i,:)=predcv.coeff.R2c;
        snvd1rpdcv(i,:)=predcv.coeff.RPDc;
        %============
        % SNV + D2
        %============
        predcv=predPLS(Xpretcvcal.SNVD2(:,selectvar), Xpretcvval.SNVD2(:,selectvar), Ycvcal, Ycvval, 0, 0, contrainte, predsnvd2.para.nbcomp);
        snvd2SEP(i,:)=predcv.coeff.SEC;
        snvd2r2p(i,:)=predcv.coeff.R2c;
        snvd2rpdcv(i,:)=predcv.coeff.RPDc;
        %============
        % SNV + detrend + D1
        %============
        predcv=predPLS(Xpretcvcal.SNVDD1(:,selectvar), Xpretcvval.SNVDD1(:,selectvar), Ycvcal, Ycvval, 0, 0, contrainte, predsnvdd1.para.nbcomp);
        snvdd1SEP(i,:)=predcv.coeff.SEC;
        snvdd1r2p(i,:)=predcv.coeff.R2c;
        snvdd1rpdcv(i,:)=predcv.coeff.RPDc;
        %============
        % SNV + detrend + D2
        %============
        predcv=predPLS(Xpretcvcal.SNVDD2(:,selectvar), Xpretcvval.SNVDD2(:,selectvar), Ycvcal, Ycvval, 0, 0, contrainte, predsnvdd2.para.nbcomp);
        snvdd2SEP(i,:)=predcv.coeff.SEC;
        snvdd2r2p(i,:)=predcv.coeff.R2c;
        snvdd2rpdcv(i,:)=predcv.coeff.RPDc;
        %============
        % CR
        %============
        predcv=predPLS(Xpretcvcal.CR(:,selectvar), Xpretcvval.CR(:,selectvar), Ycvcal, Ycvval, 0, 0, contrainte, predcr.para.nbcomp);
        crSEP(i,:)=predcv.coeff.SEC;
        crr2p(i,:)=predcv.coeff.R2c;
        crrpdcv(i,:)=predcv.coeff.RPDc;
    end
else if pretc>0
        for i=1:mcv
            setcv=Mcv(:,i);
            Scvcal=Scal(setcv==1,:);
            Scvval=Scal(setcv==0,:);
            Ycvcal=Ycal(setcv==1,:);
            Ycvval=Ycal(setcv==0,:);
            
            %============
            % PREPROCESSING
            %============
            if pretc>1
            Xpretcvcal = AllPret(Scvcal,wl,pretc-1);
            Xpretcvval = AllPret(Scvval,wl,pretc-1);            
            end
            %============
            % RAW
            %============
            if pretc==0||pretc==1
                predcv=predPLS(Scvcal(:,selectvar), Scvval(:,selectvar), Ycvcal, Ycvval, 0, 0, contrainte,predbrut.para.nbcomp);
                brutSEP(i,:)=predcv.coeff.SEC;
                brutr2p(i,:)=predcv.coeff.R2c;
                brutrpdcv(i,:)=predcv.coeff.RPDc;
            end
            %============
            % DETREND
            %============
            if pretc==2
                predcv=predPLS(Xpretcvcal(:,selectvar),Xpretcvval(:,selectvar),Ycvcal,Ycvval,0,0,contrainte,preddetrend.para.nbcomp);
                detrendSEP(i,:)=predcv.coeff.SEC;
                detrendr2p(i,:)=predcv.coeff.R2c;
                detrendrpdcv(i,:)=predcv.coeff.RPDc;
            end
            %============
            % SNV
            %============
            if pretc==3
                predcv=predPLS(Xpretcvcal(:,selectvar), Xpretcvval(:,selectvar), Ycvcal, Ycvval, 0, 0, contrainte,predsnv.para.nbcomp);
                snvcvSEP(i,:)=predcv.coeff.SEC;
                snvcvr2p(i,:)=predcv.coeff.R2c;
                snvcvrpd(i,:)=predcv.coeff.RPDc;
            end
            %============
            % MSC
            %============
            if pretc==4
                predcv=predPLS(Xpretcvcal(:,selectvar), Xpretcvval(:,selectvar), Ycvcal, Ycvval, 0, 0, contrainte,predmsc.para.nbcomp);
                msccvSEP(i,:)=predcv.coeff.SEC;
                msccvr2p(i,:)=predcv.coeff.R2c;
                msccvrpd(i,:)=predcv.coeff.RPDc;
            end
            %============
            % SNV + detrend
            %============
            if pretc==5
                predcv=predPLS(Xpretcvcal(:,selectvar), Xpretcvval(:,selectvar), Ycvcal, Ycvval, 0, 0, contrainte,predsnvd.para.nbcomp);
                snvdSEP(i,:)=predcv.coeff.SEC;
                snvdr2p(i,:)=predcv.coeff.R2c;
                snvdrpdcv(i,:)=predcv.coeff.RPDc;
            end
            %============
            % D1
            %============
            if pretc==6
                predcv=predPLS(Xpretcvcal(:,selectvar), Xpretcvval(:,selectvar), Ycvcal, Ycvval, 0, 0, contrainte, predd1.para.nbcomp);
                d1SEP(i,:)=predcv.coeff.SEC;
                d1r2p(i,:)=predcv.coeff.R2c;
                d1rpdcv(i,:)=predcv.coeff.RPDc;
            end
            %============
            % D2
            %============
            if pretc==7
                predcv=predPLS(Xpretcvcal(:,selectvar), Xpretcvval(:,selectvar), Ycvcal, Ycvval, 0, 0, contrainte, predd2.para.nbcomp);
                d2SEP(i,:)=predcv.coeff.SEC;
                d2r2p(i,:)=predcv.coeff.R2c;
                d2rpdcv(i,:)=predcv.coeff.RPDc;
            end
            %============
            % MSC + D1
            %============
            if pretc==8
                predcv=predPLS(Xpretcvcal(:,selectvar), Xpretcvval(:,selectvar), Ycvcal, Ycvval, 0, 0, contrainte, predmscd1.para.nbcomp);
                mscd1SEP(i,:)=predcv.coeff.SEC;
                mscd1r2p(i,:)=predcv.coeff.R2c;
                mscd1rpdcv(i,:)=predcv.coeff.RPDc;
            end
            %============
            % MSC + D2
            %============
            if pretc==9
                predcv=predPLS(Xpretcvcal(:,selectvar), Xpretcvval(:,selectvar), Ycvcal, Ycvval, 0, 0, contrainte, predmscd2.para.nbcomp);
                mscd2SEP(i,:)=predcv.coeff.SEC;
                mscd2r2p(i,:)=predcv.coeff.R2c;
                mscd2rpdcv(i,:)=predcv.coeff.RPDc;
            end
            %============
            % SNV + D1
            %============
            if pretc==10
                predcv=predPLS(Xpretcvcal(:,selectvar), Xpretcvval(:,selectvar), Ycvcal, Ycvval, 0, 0, contrainte, predsnvd1.para.nbcomp);
                snvd1SEP(i,:)=predcv.coeff.SEC;
                snvd1r2p(i,:)=predcv.coeff.R2c;
                snvd1rpdcv(i,:)=predcv.coeff.RPDc;
            end
            %============
            % SNV + D2
            %============
            if pretc==11
                predcv=predPLS(Xpretcvcal(:,selectvar), Xpretcvval(:,selectvar), Ycvcal, Ycvval, 0, 0, contrainte, predsnvd2.para.nbcomp);
                snvd2SEP(i,:)=predcv.coeff.SEC;
                snvd2r2p(i,:)=predcv.coeff.R2c;
                snvd2rpdcv(i,:)=predcv.coeff.RPDc;
            end
            %============
            % SNV + detrend + D1
            %============
            if pretc==12
                predcv=predPLS(Xpretcvcal(:,selectvar), Xpretcvval(:,selectvar), Ycvcal, Ycvval, 0, 0, contrainte, predsnvdd1.para.nbcomp);
                snvdd1SEP(i,:)=predcv.coeff.SEC;
                snvdd1r2p(i,:)=predcv.coeff.R2c;
                snvdd1rpdcv(i,:)=predcv.coeff.RPDc;
            end
            %============
            % SNV + detrend + D2
            %============
            if pretc==13
                predcv=predPLS(Xpretcvcal(:,selectvar), Xpretcvval(:,selectvar), Ycvcal, Ycvval, 0, 0, contrainte, predsnvdd2.para.nbcomp);
                snvdd2SEP(i,:)=predcv.coeff.SEC;
                snvdd2r2p(i,:)=predcv.coeff.R2c;
                snvdd2rpdcv(i,:)=predcv.coeff.RPDc;
            end
            %============
            % CR
            %============
            if pretc==14
                predcv=predPLS(Xpretcvcal(:,selectvar), Xpretcvval(:,selectvar), Ycvcal, Ycvval, 0, 0, contrainte, predcr.para.nbcomp);
                crSEP(i,:)=predcv.coeff.SEC;
                crr2p(i,:)=predcv.coeff.R2c;
                crrpdcv(i,:)=predcv.coeff.RPDc;
            end
        end
    end
end

% Mean cross-validation performance
if size(Y,2)==1
    if pretc==0||pretc==1
        mcvi=length(brutSEP(~isnan(brutr2p)));
        bruts=sqrt(sum(brutSEP(~isnan(brutr2p)).^2)/mcvi);
        brutr=mean(brutr2p(~isnan(brutr2p)));
        brutrpdcvc=mean(brutrpdcv(~isnan(brutr2p)));
        
    end
    if pretc==0||pretc==2
        mcvi=length(detrendSEP(~isnan(detrendr2p)));
        detrends=sqrt(sum(detrendSEP(~isnan(detrendr2p)).^2)/mcvi);
        detrendr=mean(detrendr2p(~isnan(detrendr2p)));
        detrendrpdcvc=mean(detrendrpdcv(~isnan(detrendr2p)));
        
    end
    if pretc==0||pretc==3
        mcvi=length(snvcvSEP(~isnan(snvcvr2p)));
        snvs=sqrt(sum(snvcvSEP(~isnan(snvcvr2p)).^2)/mcvi);
        snvr=mean(snvcvr2p(~isnan(snvcvr2p)));
        snvrpdcvc=mean(snvcvrpd(~isnan(snvcvr2p)));
        
    end
    if pretc==0||pretc==4
        mcvi=length(msccvSEP(~isnan(msccvr2p)));
        mscs=sqrt(sum(msccvSEP(~isnan(msccvr2p)).^2)/mcvi);
        mscr=mean(msccvr2p(~isnan(msccvr2p)));
        mscrpdcvc=mean(msccvrpd(~isnan(msccvr2p)));
        
    end
    if pretc==0||pretc==5
        mcvi=length(snvdSEP(~isnan(snvdr2p)));
        snvds=sqrt(sum(snvdSEP(~isnan(snvdr2p)).^2)/mcvi);
        snvdr=mean(snvdr2p(~isnan(snvdr2p)));
        snvdrpdcvc=mean(snvdrpdcv(~isnan(snvdr2p)));
        
    end
    if pretc==0||pretc==6
        mcvi=length(d1SEP(~isnan(d1r2p)));
        d1s=sqrt(sum(d1SEP(~isnan(d1r2p)).^2)/mcvi);
        d1r=mean(d1r2p(~isnan(d1r2p)));
        d1rpdcvc=mean(d1rpdcv(~isnan(d1r2p)));
        
    end
    if pretc==0||pretc==7
        mcvi=length(d2SEP(~isnan(d2r2p)));
        d2s=sqrt(sum(d2SEP(~isnan(d2r2p)).^2)/mcvi);
        d2r=mean(d2r2p(~isnan(d2r2p)));
        d2rpdcvc=mean(d2rpdcv(~isnan(d2r2p)));
        
    end
    if pretc==0||pretc==8
        mcvi=length(mscd1SEP(~isnan(mscd1r2p)));
        mscd1s=sqrt(sum(mscd1SEP(~isnan(mscd1r2p)).^2)/mcvi);
        mscd1r=mean(mscd1r2p(~isnan(mscd1r2p)));
        mscd1rpdcvc=mean(mscd1rpdcv(~isnan(mscd1r2p)));
        
    end
    if pretc==0||pretc==9
        mcvi=length(mscd2SEP(~isnan(mscd2r2p)));
        mscd2s=sqrt(sum(mscd2SEP(~isnan(mscd2r2p)).^2)/mcvi);
        mscd2r=mean(mscd2r2p(~isnan(mscd2r2p)));
        mscd2rpdcvc=mean(mscd2rpdcv(~isnan(mscd2r2p)));
        
    end
    if pretc==0||pretc==10
        mcvi=length(snvd1SEP(~isnan(snvd1r2p)));
        snvd1s=sqrt(sum(snvd1SEP(~isnan(snvd1r2p)).^2)/mcvi);
        snvd1r=mean(snvd1r2p(~isnan(snvd1r2p)));
        snvd1rpdcvc=mean(snvd1rpdcv(~isnan(snvd1r2p)));
        
    end
    if pretc==0||pretc==11
        mcvi=length(snvd2SEP(~isnan(snvd2r2p)));
        snvd2s=sqrt(sum(snvd2SEP(~isnan(snvd2r2p)).^2)/mcvi);
        snvd2r=mean(snvd2r2p(~isnan(snvd2r2p)));
        snvd2rpdcvc=mean(snvd2rpdcv(~isnan(snvd2r2p)));
        
    end
    if pretc==0||pretc==12
        mcvi=length(snvdd1SEP(~isnan(snvdd1r2p)));
        snvdd1s=sqrt(sum(snvdd1SEP(~isnan(snvdd1r2p)).^2)/mcvi);
        snvdd1r=mean(snvdd1r2p(~isnan(snvdd1r2p)));
        snvdd1rpdcvc=mean(snvdd1rpdcv(~isnan(snvdd1r2p)));
        
    end
    if pretc==0||pretc==13
        mcvi=length(snvdd2SEP(~isnan(snvdd2r2p)));
        snvdd2s=sqrt(sum(snvdd2SEP(~isnan(snvdd2r2p)).^2)/mcvi);
        snvdd2r=mean(snvdd2r2p(~isnan(snvdd2r2p)));
        snvdd2rpdcvc=mean(snvdd2rpdcv(~isnan(snvdd2r2p)));
        
    end
    if pretc==0||pretc==14
        mcvi=length(crSEP(~isnan(crr2p)));
        crs=sqrt(sum(crSEP(~isnan(crr2p)).^2)/mcvi);
        crr=mean(crr2p(~isnan(crr2p)));
        crrpdcvc=mean(crrpdcv(~isnan(crr2p)));
        
    end
    
else % TO IMPROVE!!!!
    for i=1:size(Y,2)
        if pretc==0||pretc==1
            mcvi=length(brutSEP(brutr2p( :,i)>0,i));
            bruts(i)=sqrt(sum(brutSEP(brutr2p( :,i)>0,i).^2)/mcvi);
            brutr(i)=mean(brutr2p(brutr2p( :,i)>0,i));
            brutrpdcvc(i)=mean(brutrpdcv(brutr2p( :,i)>0,i));
        end
        if pretc==0||pretc==2
            mcvi=length(detrendSEP(detrendr2p( :,i)>0,i));
            detrends(i)=sqrt(sum(detrendSEP(detrendr2p( :,i)>0,i).^2)/mcvi);
            detrendr(i)=mean(detrendr2p(detrendr2p( :,i)>0,i));
            detrendrpdcvc(i)=mean(detrendrpdcv(detrendr2p( :,i)>0,i));
        end
        if pretc==0||pretc==3
            mcvi=length(snvcvSEP(snvcvr2p( :,i)>0,i));
            snvs(i)=sqrt(sum(snvcvSEP(snvcvr2p( :,i)>0,i).^2)/mcvi);
            snvr(i)=mean(snvcvr2p(snvcvr2p( :,i)>0,i));
            snvrpdcvc(i)=mean(snvcvrpd(snvcvr2p( :,i)>0,i));
        end
        if pretc==0||pretc==4
            mcvi=length(msccvSEP(msccvr2p( :,i)>0,i));
            mscs(i)=sqrt(sum(msccvSEP(msccvr2p( :,i)>0,i).^2)/mcvi);
            mscr(i)=mean(msccvr2p(msccvr2p( :,i)>0,i));
            mscrpdcvc(i)=mean(msccvrpd(msccvr2p( :,i)>0,i));
        end
        if pretc==0||pretc==5
            mcvi=length(snvdSEP(snvdr2p( :,i)>0,i));
            snvds(i)=sqrt(sum(snvdSEP(snvdr2p( :,i)>0,i).^2)/mcvi);
            snvdr(i)=mean(snvdr2p(snvdr2p( :,i)>0 ,i));
            snvdrpdcvc(i)=mean(snvdrpdcv(snvdr2p( :,i)>0,i));
        end
        if pretc==0||pretc==6
            mcvi=length(d1SEP(d1r2p( :,i)>0,i));
            d1s(i)=sqrt(sum(d1SEP(d1r2p( :,i)>0,i).^2)/mcvi);
            d1r(i)=mean(d1r2p(d1r2p( :,i)>0,i));
            d1rpdcvc(i)=mean(d1rpdcv(d1r2p( :,i)>0,i));
        end
        if pretc==0||pretc==7
            mcvi=length(d2SEP(d2r2p( :,i)>0,i));
            d2s(i)=sqrt(sum(d2SEP(d2r2p( :,i)>0,i).^2)/mcvi);
            d2r(i)=mean(d2r2p(d2r2p( :,i)>0,i));
            d2rpdcvc(i)=mean(d2rpdcv(d2r2p( :,i)>0,i));
        end
        if pretc==0||pretc==8
            mcvi=length(mscd1SEP(mscd1r2p( :,i)>0));
            mscd1s(i)=sqrt(sum(mscd1SEP(mscd1r2p( :,i)>0,i).^2)/mcvi);
            mscd1r(i)=mean(mscd1r2p(mscd1r2p( :,i)>0,i));
            mscd1rpdcvc(i)=mean(mscd1rpdcv(mscd1r2p( :,i)>0,i));
        end
        if pretc==0||pretc==9
            mcvi=length(mscd2SEP(mscd2r2p( :,i)>0,i));
            mscd2s(i)=sqrt(sum(mscd2SEP(mscd2r2p( :,i)>0,i).^2)/mcvi);
            mscd2r(i)=mean(mscd2r2p(mscd2r2p( :,i)>0,i));
            mscd2rpdcvc(i)=mean(mscd2rpdcv(mscd2r2p( :,i)>0,i));
        end
        if pretc==0||pretc==10
            mcvi=length(snvd1SEP(snvd1r2p( :,i)>0,i));
            snvd1s(i)=sqrt(sum(snvd1SEP(snvd1r2p( :,i)>0,i).^2)/mcvi);
            snvd1r(i)=mean(snvd1r2p(snvd1r2p( :,i)>0,i));
            snvd1rpdcvc(i)=mean(snvd1rpdcv(snvd1r2p( :,i)>0,i));
        end
        if pretc==0||pretc==11
            mcvi=length(snvd2SEP(snvd2r2p( :,i)>0,i));
            snvd2s(i)=sqrt(sum(snvd2SEP(snvd2r2p( :,i)>0,i).^2)/mcvi);
            snvd2r(i)=mean(snvd2r2p(snvd2r2p( :,i)>0,i));
            snvd2rpdcvc(i)=mean(snvd2rpdcv(snvd2r2p( :,i)>0,i));
        end
        if pretc==0||pretc==12
            mcvi=length(snvdd1SEP(snvdd1r2p( :,i)>0,i));
            snvdd1s(i)=sqrt(sum(snvdd1SEP(snvdd1r2p( :,i)>0,i).^2)/mcvi);
            snvdd1r(i)=mean(snvdd1r2p(snvdd1r2p( :,i)>0,i));
            snvdd1rpdcvc(i)=mean(snvdd1rpdcv(snvdd1r2p( :,i)>0,i));
        end
        if pretc==0||pretc==13
            mcvi=length(snvdd2SEP(snvdd2r2p( :,i)>0,i));
            snvdd2s(i)=sqrt(sum(snvdd2SEP(snvdd2r2p( :,i)>0,i).^2)/mcvi);
            snvdd2r(i)=mean(snvdd2r2p(snvdd2r2p( :,i)>0,i));
            snvdd2rpdcvc(i)=mean(snvdd2rpdcv(snvdd2r2p( :,i)>0,i));
        end
        if pretc==0||pretc==14
            mcvi=length(crSEP(crr2p( :,i)>0,i));
            crs(i)=sqrt(sum(crSEP(crr2p( :,i)>0,i).^2)/mcvi);
            crr(i)=mean(crr2p(crr2p( :,i)>0,i));
            crrpdcvc(i)=mean(crrpdcv(crr2p( :,i)>0,i));
        end
        
    end
end

%============
% SUMMARY
%============
if pretc==0
    [nc, ~]=size(Scal);
    [nv, ~]=size(Sval);
    
    [~, yk]=size(Y);
    
    C1=[ones(1,yk)*nc ;min(Ycal) ;mean(Ycal) ;max(Ycal) ;std(Ycal)];
    P1=[ones(1,yk)*nv; min(Yval); mean(Yval); max(Yval); std(Yval)];
    pret={'Raw' 'Detrend' 'MSC' 'SNV' 'SNVD' 'D1' 'D2' 'MSC+D1' 'SNV+D1' 'SNVD+D1' 'MSC+D2' 'SNV+D2' 'SNVD+D2' 'CR'}';
    column={'Pret' 'NbVL' 'SEC' 'R²c' 'RPDc' 'SECV' 'R²cv' 'RPDcv' 'SEP' 'RMSEP' 'R²p' 'RPDp'};
    Pretc=[predbrut.para.nbcomp mean(predbrut.coeff.SEC) mean(predbrut.coeff.R2c) mean(predbrut.coeff.RPDc);
        preddetrend.para.nbcomp mean(preddetrend.coeff.SEC) mean(preddetrend.coeff.R2c) mean(preddetrend.coeff.RPDc);
        predmsc.para.nbcomp mean(predmsc.coeff.SEC) mean(predmsc.coeff.R2c) mean(predmsc.coeff.RPDc);
        predsnv.para.nbcomp mean(predsnv.coeff.SEC) mean(predsnv.coeff.R2c) mean(predsnv.coeff.RPDc);
        predsnvd.para.nbcomp mean(predsnvd.coeff.SEC) mean(predsnvd.coeff.R2c) mean(predsnvd.coeff.RPDc);
        predd1.para.nbcomp mean(predd1.coeff.SEC) mean(predd1.coeff.R2c) mean(predd1.coeff.RPDc);
        predd2.para.nbcomp mean(predd2.coeff.SEC) mean(predd2.coeff.R2c) mean(predd2.coeff.RPDc);
        predmscd1.para.nbcomp mean(predmscd1.coeff.SEC) mean(predmscd1.coeff.R2c) mean(predmscd1.coeff.RPDc);
        predsnvd1.para.nbcomp mean(predsnvd1.coeff.SEC) mean(predsnvd1.coeff.R2c) mean(predsnvd1.coeff.RPDc);
        predsnvdd1.para.nbcomp mean(predsnvdd1.coeff.SEC) mean(predsnvdd1.coeff.R2c) mean(predsnvdd1.coeff.RPDc);
        predmscd2.para.nbcomp mean(predmscd2.coeff.SEC) mean(predmscd2.coeff.R2c) mean(predmscd2.coeff.RPDc);
        predsnvd2.para.nbcomp mean(predsnvd2.coeff.SEC) mean(predsnvd2.coeff.R2c) mean(predsnvd2.coeff.RPDc);
        predsnvdd2.para.nbcomp mean(predsnvdd2.coeff.SEC) mean(predsnvdd2.coeff.R2c) mean(predsnvdd2.coeff.RPDc);
        predcr.para.nbcomp mean(predcr.coeff.SEC) mean(predcr.coeff.R2c) mean(predcr.coeff.RPDc)];
    Pretcv=[mean(bruts) mean(brutr) mean(brutrpdcvc);
        mean(detrends) mean(detrendr) mean(detrendrpdcvc);
        mean(mscs) mean(mscr) mean(mscrpdcvc);
        mean(snvs) mean(snvr) mean(snvrpdcvc);
        mean(snvds) mean(snvdr) mean(snvdrpdcvc);
        mean(d1s) mean(d1r) mean(d1rpdcvc);
        mean(d2s) mean(d2r) mean(d2rpdcvc);
        mean(mscd1s) mean(mscd1r) mean(mscd1rpdcvc);
        mean(snvd1s) mean(snvd1r) mean(snvd1rpdcvc);
        mean(snvdd1s) mean(snvdd1r) mean(snvdd1rpdcvc);
        mean(mscd2s) mean(mscd2r) mean(mscd2rpdcvc);
        mean(snvd2s) mean(snvd2r) mean(snvd2rpdcvc);
        mean(snvdd2s) mean(snvdd2r) mean(snvdd2rpdcvc);
        mean(crs) mean(crr) mean(crrpdcvc)];
    Pretp=[mean(predbrut.coeff.SEP) mean(predbrut.coeff.RMSEP) mean(predbrut.coeff.R2p) mean(predbrut.coeff.RPDp);
        mean(preddetrend.coeff.SEP) mean(preddetrend.coeff.RMSEP) mean(preddetrend.coeff.R2p) mean(preddetrend.coeff.RPDp);
        mean(predmsc.coeff.SEP) mean(predmsc.coeff.RMSEP) mean(predmsc.coeff.R2p) mean(predmsc.coeff.RPDp);
        mean(predsnv.coeff.SEP) mean(predsnv.coeff.RMSEP) mean(predsnv.coeff.R2p) mean(predsnv.coeff.RPDp);
        mean(predsnvd.coeff.SEP) mean(predsnvd.coeff.RMSEP) mean(predsnvd.coeff.R2p) mean(predsnvd.coeff.RPDp);
        mean(predd1.coeff.SEP) mean(predd1.coeff.RMSEP) mean(predd1.coeff.R2p) mean(predd1.coeff.RPDp);
        mean(predd2.coeff.SEP) mean(predd2.coeff.RMSEP) mean(predd2.coeff.R2p) mean(predd2.coeff.RPDp);
        mean(predmscd1.coeff.SEP) mean(predmscd1.coeff.RMSEP) mean(predmscd1.coeff.R2p) mean(predmscd1.coeff.RPDp);
        mean(predsnvd1.coeff.SEP) mean(predsnvd1.coeff.RMSEP) mean(predsnvd1.coeff.R2p) mean(predsnvd1.coeff.RPDp);
        mean(predsnvdd1.coeff.SEP) mean(predsnvdd1.coeff.RMSEP) mean(predsnvdd1.coeff.R2p) mean(predsnvdd1.coeff.RPDp);
        mean(predmscd2.coeff.SEP) mean(predmscd2.coeff.RMSEP) mean(predmscd2.coeff.R2p) mean(predmscd2.coeff.RPDp);
        mean(predsnvd2.coeff.SEP)  mean(predsnvd2.coeff.RMSEP) mean(predsnvd2.coeff.R2p) mean(predsnvd2.coeff.RPDp);
        mean(predsnvdd2.coeff.SEP) mean(predsnvdd2.coeff.RMSEP) mean(predsnvdd2.coeff.R2p) mean(predsnvdd2.coeff.RPDp);
        mean(predcr.coeff.SEP) mean(predcr.coeff.RMSEP) mean(predcr.coeff.R2p) mean(predcr.coeff.RPDp)];
    
    recap(2:15,1)=pret;
    recap(1,1:12)=column;
    recap(2:15,2:5)=mat2cell(Pretc,ones(1,size(Pretc,1)),ones(1,size(Pretc,2)));
    recap(2:15,6:8)=mat2cell(Pretcv,ones(1,size(Pretcv,1)),ones(1,size(Pretcv,2)));
    recap(2:15,9:12)=mat2cell(Pretp,ones(1,size(Pretp,1)),ones(1,size(Pretp,2)));
    
    if size(Y,2)>1
        for i=1:size(Y,2)
            Pretc=[predbrut.para.nbcomp predbrut.coeff.SEC(1,i) predbrut.coeff.R2c(1,i) predbrut.coeff.RPDc(1,i);
                preddetrend.para.nbcomp preddetrend.coeff.SEC(1,i) preddetrend.coeff.R2c(1,i) preddetrend.coeff.RPDc(1,i);
                predmsc.para.nbcomp predmsc.coeff.SEC(1,i) predmsc.coeff.R2c(1,i) predmsc.coeff.RPDc(1,i);
                predsnv.para.nbcomp predsnv.coeff.SEC(1,i) predsnv.coeff.R2c(1,i) predsnv.coeff.RPDc(1,i);
                predsnvd.para.nbcomp predsnvd.coeff.SEC(1,i) predsnvd.coeff.R2c(1,i) predsnvd.coeff.RPDc(1,i);
                predd1.para.nbcomp predd1.coeff.SEC(1,i) predd1.coeff.R2c(1,i) predd1.coeff.RPDc(1,i);
                predd2.para.nbcomp predd2.coeff.SEC(1,i) predd2.coeff.R2c(1,i) predd2.coeff.RPDc(1,i);
                predmscd1.para.nbcomp predmscd1.coeff.SEC(1,i) predmscd1.coeff.R2c(1,i) predmscd1.coeff.RPDc(1,i);
                predsnvd1.para.nbcomp predsnvd1.coeff.SEC(1,i) predsnvd1.coeff.R2c(1,i) predsnvd1.coeff.RPDc(1,i);
                predsnvdd1.para.nbcomp predsnvdd1.coeff.SEC(1,i) predsnvdd1.coeff.R2c(1,i) predsnvdd1.coeff.RPDc(1,i);
                predmscd2.para.nbcomp predmscd2.coeff.SEC(1,i) predmscd2.coeff.R2c(1,i) predmscd2.coeff.RPDc(1,i);
                predsnvd2.para.nbcomp predsnvd2.coeff.SEC(1,i) predsnvd2.coeff.R2c(1,i) predsnvd2.coeff.RPDc(1,i);
                predsnvdd2.para.nbcomp predsnvdd2.coeff.SEC(1,i) predsnvdd2.coeff.R2c(1,i) predsnvdd2.coeff.RPDc(1,i);
                predcr.para.nbcomp predcr.coeff.SEC(1,i) predcr.coeff.R2c(1,i) predcr.coeff.RPDc(1,i)];
            Pretcv=[bruts(1,i) brutr(1,i) brutrpdcvc(1,i);
                detrends(1,i) detrendr(1,i) detrendrpdcvc(1,i);
                mscs(1,i) mscr(1,i) mscrpdcvc(1,i);
                snvs(1,i) snvr(1,i) snvrpdcvc(1,i);
                snvds(1,i) snvdr(1,i) snvdrpdcvc(1,i);
                d1s(1,i) d1r(1,i) d1rpdcvc(1,i);
                d2s(1,i) d2r(1,i) d2rpdcvc(1,i);
                mscd1s(1,i) mscd1r(1,i) mscd1rpdcvc(1,i);
                snvd1s(1,i) snvd1r(1,i) snvd1rpdcvc(1,i);
                snvdd1s(1,i) snvdd1r(1,i) snvdd1rpdcvc(1,i);
                mscd2s(1,i) mscd2r(1,i) mscd2rpdcvc(1,i);
                snvd2s(1,i) snvd2r(1,i) snvd2rpdcvc(1,i);
                snvdd2s(1,i) snvdd2r(1,i) snvdd2rpdcvc(1,i);
                crs(1,i) crr(1,i) crrpdcvc(1,i)];
            Pretp=[predbrut.coeff.SEP(1,i) predbrut.coeff.RMSEP(1,i) predbrut.coeff.R2p(1,i) predbrut.coeff.RPDp(1,i);
                preddetrend.coeff.SEP(1,i) preddetrend.coeff.RMSEP(1,i) preddetrend.coeff.R2p(1,i) preddetrend.coeff.RPDp(1,i);
                predmsc.coeff.SEP(1,i) predmsc.coeff.RMSEP(1,i) predmsc.coeff.R2p(1,i) predmsc.coeff.RPDp(1,i);
                predsnv.coeff.SEP(1,i) predsnv.coeff.RMSEP(1,i) predsnv.coeff.R2p(1,i) predsnv.coeff.RPDp(1,i);
                predsnvd.coeff.SEP(1,i) predsnvd.coeff.RMSEP(1,i) predsnvd.coeff.R2p(1,i) predsnvd.coeff.RPDp(1,i);
                predd1.coeff.SEP(1,i) predd1.coeff.RMSEP(1,i) predd1.coeff.R2p(1,i) predd1.coeff.RPDp(1,i);
                predd2.coeff.SEP(1,i) predd2.coeff.RMSEP(1,i) predd2.coeff.R2p(1,i) predd2.coeff.RPDp(1,i);
                predmscd1.coeff.SEP(1,i) predmscd1.coeff.RMSEP(1,i) predmscd1.coeff.R2p(1,i) predmscd1.coeff.RPDp(1,i);
                predsnvd1.coeff.SEP(1,i) predsnvd1.coeff.RMSEP(1,i) predsnvd1.coeff.R2p(1,i) predsnvd1.coeff.RPDp(1,i);
                predsnvdd1.coeff.SEP(1,i) predsnvdd1.coeff.RMSEP(1,i) predsnvdd1.coeff.R2p(1,i) predsnvdd1.coeff.RPDp(1,i);
                predmscd2.coeff.SEP(1,i) predmscd2.coeff.RMSEP(1,i) predmscd2.coeff.R2p(1,i) predmscd2.coeff.RPDp(1,i);
                predsnvd2.coeff.SEP(1,i) predsnvd2.coeff.RMSEP(1,i) predsnvd2.coeff.R2p(1,i) predsnvd2.coeff.RPDp(1,i);
                predsnvdd2.coeff.SEP(1,i) predsnvdd2.coeff.RMSEP(1,i) predsnvdd2.coeff.R2p(1,i) predsnvdd2.coeff.RPDp(1,i);
                predcr.coeff.SEP(1,i) predcr.coeff.RMSEP(1,i) predcr.coeff.R2p(1,i) predcr.coeff.RPDp(1,i)];
            
            recapi(2:15,1)=pret;
            recapi(1,1:12)=column;
            recapi(2:15,2:5)=mat2cell(Pretc,ones(1,size(Pretc,1)),ones(1,size(Pretc,2)));
            recapi(2:15,6:8)=mat2cell(Pretcv,ones(1,size(Pretcv,1)),ones(1,size(Pretcv,2)));
            recapi(2:15,9:12)=mat2cell(Pretp,ones(1,size(Pretp,1)),ones(1,size(Pretp,2)));
            
            RECAP{i,1}=recapi;
        end
    end
end

%============
% SAVE
%============
if pretc==0
    pred.recap=recap;
end
if size(Y,2)>1&&pretc==0
    pred.RECAP=RECAP;
end
if pretc==0||pretc==1
    pred.brut=predbrut;
    pred.brut.cv.SECV=bruts;
    pred.brut.cv.r2cv=brutr;
    pred.brut.cv.rpdcv=brutrpdcvc;
end

if pretc==0
    pred.detrend=preddetrend;
    pred.detrend.cv.SECV=detrends;
    pred.detrend.cv.r2cv=detrendr;
    pred.detrend.cv.rpdcv=detrendrpdcv;
end
if pretc==2
    pred.brut=preddetrend;
    pred.brut.cv.SECV=detrends;
    pred.brut.cv.r2cv=detrendr;
    pred.brut.cv.rpdcv=detrendrpdcv;
end

if pretc==0
    pred.msc=predmsc;
    pred.msc.cv.SECV=mscs;
    pred.msc.cv.r2cv=mscr;
    pred.msc.cv.rpdcv=mscrpdcvc;
end
if pretc==4
    pred.brut=predmsc;
    pred.brut.cv.SECV=mscs;
    pred.brut.cv.r2cv=mscr;
    pred.brut.cv.rpdcv=mscrpdcvc;
end

if pretc==0
    pred.snv=predsnv;
    pred.snv.cv.SECV=snvs;
    pred.snv.cv.r2cv=snvr;
    pred.snv.cv.rpdcv=snvrpdcvc;
end
if pretc==3
    pred.brut=predsnv;
    pred.brut.cv.SECV=snvs;
    pred.brut.cv.r2cv=snvr;
    pred.brut.cv.rpdcv=snvrpdcvc;
end

if pretc==0
    pred.snvd=predsnvd;
    pred.snvd.cv.SECV=snvds;
    pred.snvd.cv.r2cv=snvdr;
    pred.snvd.cv.rpdcv=snvdrpdcv;
end
if pretc==5
    pred.brut=predsnvd;
    pred.brut.cv.SECV=snvds;
    pred.brut.cv.r2cv=snvdr;
    pred.brut.cv.rpdcv=snvdrpdcv;
end

if pretc==0
    pred.d1=predd1;
    pred.d1.cv.SECV=d1s;
    pred.d1.cv.r2cv=d1r;
    pred.d1.cv.rpdcv=d1rpdcv;
end
if pretc==6
    pred.brut=predd1;
    pred.brut.cv.SECV=d1s;
    pred.brut.cv.r2cv=d1r;
    pred.brut.cv.rpdcv=d1rpdcv;
end

if pretc==0
    pred.d2=predd2;
    pred.d2.cv.SECV=d2s;
    pred.d2.cv.r2cv=d2r;
    pred.d2.cv.rpdcv=d2rpdcv;
end
if pretc==7
    pred.brut=predd2;
    pred.brut.cv.SECV=d2s;
    pred.brut.cv.r2cv=d2r;
    pred.brut.cv.rpdcv=d2rpdcv;
end

if pretc==0
    pred.mscd1=predmscd1;
    pred.mscd1.cv.SECV=mscd1s;
    pred.mscd1.cv.r2cv=mscd1r;
    pred.mscd1.cv.rpdcv=mscd1rpdcv;
end
if pretc==8
    pred.brut=predmscd1;
    pred.brut.cv.SECV=mscd1s;
    pred.brut.cv.r2cv=mscd1r;
    pred.brut.cv.rpdcv=mscd1rpdcv;
end

if pretc==0
    pred.mscd2=predmscd2;
    pred.mscd2.cv.SECV=mscd2s;
    pred.mscd2.cv.r2cv=mscd2r;
    pred.mscd2.cv.rpdcv=mscd2rpdcv;
end
if pretc==9
    pred.brut=predmscd2;
    pred.brut.cv.SECV=mscd2s;
    pred.brut.cv.r2cv=mscd2r;
    pred.brut.cv.rpdcv=mscd2rpdcv;
end

if pretc==0
    pred.snvd1=predsnvd1;
    pred.snvd1.cv.SECV=snvd1s;
    pred.snvd1.cv.r2cv=snvd1r;
    pred.snvd1.cv.rpdcv=snvd1rpdcv;
end
if pretc==10
    pred.brut=predsnvd1;
    pred.brut.cv.SECV=snvd1s;
    pred.brut.cv.r2cv=snvd1r;
    pred.brut.cv.rpdcv=snvd1rpdcv;
end

if pretc==0
    pred.snvd2=predsnvd2;
    pred.snvd2.cv.SECV=snvd2s;
    pred.snvd2.cv.r2cv=snvd2r;
    pred.snvd2.cv.rpdcv=snvd2rpdcv;
end
if pretc==11
    pred.brut=predsnvd2;
    pred.brut.cv.SECV=snvd2s;
    pred.brut.cv.r2cv=snvd2r;
    pred.brut.cv.rpdcv=snvd2rpdcv;
end

if pretc==0
    pred.snvdd1=predsnvdd1;
    pred.snvdd1.cv.SECV=snvdd1s;
    pred.snvdd1.cv.r2cv=snvdd1r;
    pred.snvdd1.cv.rpdcv=snvdd1rpdcv;
end
if pretc==12
    pred.brut=predsnvdd1;
    pred.brut.cv.SECV=snvdd1s;
    pred.brut.cv.r2cv=snvdd1r;
    pred.brut.cv.rpdcv=snvdd1rpdcv;
end

if pretc==0
    pred.snvdd2=predsnvdd2;
    pred.snvdd2.cv.SECV=snvdd2s;
    pred.snvdd2.cv.r2cv=snvdd2r;
    pred.snvdd2.cv.rpdcv=snvdd2rpdcv;
end
if pretc==13
    pred.brut=predsnvdd2;
    pred.brut.cv.SECV=snvdd2s;
    pred.brut.cv.r2cv=snvdd2r;
    pred.brut.cv.rpdcv=snvdd2rpdcv;
end

if pretc==0
    pred.cr=predcr;
    pred.cr.cv.SECV=crs;
    pred.cr.cv.r2cv=crr;
    pred.cr.cv.rpdcv=crrpdcvc;
end
if pretc==14
    pred.brut=predcr;
    pred.brut.cv.SECV=crs;
    pred.brut.cv.r2cv=crr;
    pred.brut.cv.rpdcv=crrpdcvc;
end

end