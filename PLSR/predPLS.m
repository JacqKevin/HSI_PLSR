function pred=predPLS(Scal, Sval, ref_cal, ref_val, model, fig, constraint, vl)
% PLSR estimation thanks to calibration and validation datasets.

% INPUT:    Scal,Sval: Calibration and validation spectra
%           ref_cal, ref_val: Calibration and validation reference
%           model: 0 if a model has to be created, 1 if a model already
%           exist
%           fig: figure output (0 none, 1 to display)
%           constraint: if there is prediction range to to constrain (0:
%           none, [min max])
%           vl: to constrain the number of latent variables

% OUTPUT:   pred with:
%               .para.T,P,U,Q,B,W: PLSR parameters
%               .coeff : PLSR performances
%               .Pred : Prediction and reference vectors

if nargin<8
    vl=0;
end

Sc_cal=Centrer(Scal);
Sc_val=Centrerval(Sval,Scal);

Y_cal = ref_cal;
Y_val = ref_val;

Y_calb=Y_cal;
Y_valb=Y_val;

Y_cal=Centrer(Y_cal);

if nargin<7
    constraint=0;
end

if model==1
    T=para.PLSDA.T;
    P=para.PLSDA.P;
    Q=para.PLSDA.Q;
    B=para.PLSDA.B;
    W=para.PLSDA.W;
    Wstar=para.PLSDA.Wstar;
    R2X=para.PLSDA.R2X;
    R2Y=para.PLSDA.R2Y;
else
    [B,Wstar,T,P,Q,W,R2X,R2Y,RMSE]=plsnipals(Sc_cal,Y_cal,vl);
    R2X=1-R2X;
    R2Y=1-R2Y;
end

% CALIBRATION
[ac, ~]=size(Y_calb);
Y0 = Sc_cal*Wstar*Q;
Y0 = Decentrerval( Y0, Y_calb );
if constraint(1,1)>0
    [n, ~]=size(Y0);
    for i=1:n
        if constraint(1,1)==1||constraint(1,1)==3
            if Y0(i)<constraint(1,2)
                Y0(i)=constraint(1,2);
            else
                Y0(i)=Y0(i);
            end
        end
        if constraint(1,1)==2||constraint(1,1)==3
            if Y0(i)>constraint(1,3)
                Y0(i)=constraint(1,3);
            else
                Y0(i)=Y0(i);
            end
        end
    end
end
SEC=sqrt(sum((Y0-Y_calb).^2)/ac);
[~, s2]=size(Y_calb);
z1=zeros(1,s2);
for i=1:s2
    z1(:,i)=sum((Y_calb(:,i)-mean(Y_calb(:,i))).^2);
end
z2=sum((Y0-Y_calb).^2);
SDc=std(Y_calb);
R2c=1-(z2./z1);
RPDc=1./sqrt(1-R2c);

% PREDICTION
[av, ~]=size(Y_valb);
Y1 = Sc_val*Wstar*Q;
Y1 = Decentrerval( Y1, Y_calb );
if constraint(1,1)>0
    [n, ~]=size(Y1);
    for i=1:n
        if constraint(1,1)==1||constraint(1,1)==3
            if Y1(i)<constraint(1,2)
                Y1(i)=constraint(1,2);
            else
                Y1(i)=Y1(i);
            end
        end
        if constraint(1,1)==2||constraint(1,1)==3
            if Y1(i)>constraint(1,3)
                Y1(i)=constraint(1,3);
            else
                Y1(i)=Y1(i);
            end
        end
    end
end
SEP=sqrt(sum((Y1-Y_valb).^2)/av);
[~, s2]=size(Y_valb);
z1=zeros(1,s2);
if size(Y_valb,1)>1
    for i=1:s2
        z1(:,i)=sum((Y_valb(:,i)-mean(Y_valb(:,i))).^2);
    end
else
    for i=1:s2
        z1(:,i)=sum((Y_valb(:,i)-mean(Y_calb(:,i))).^2);
    end
end
[n, ~]=size(Y1);
z2=sum((Y1-Y_valb).^2);
R2p=1-(z2./z1);
SDp=std(Y_valb);
Bias=sum(Y1-Y_valb)/n;
RMSEP=sqrt(SEP.^2+Bias.^2);
RSD=z2/(n-2);
RPDp=SDp./SEP;
[~, yi]=size(Y_valb);
Intercept=zeros(1,yi);
Slope=zeros(1,yi);
for i=1:yi
    mdl = fitlm(Y_valb(:,yi),Y1(:,yi));
    warning('off','all')
    Coeff=table2array(mdl.Coefficients);
    Intercept(1,i)=round(Coeff(1,1));
    Slope(1,i)=Coeff(2,1);
end

% GRAPHICAL OUTPUT
if fig>0
    figure(fig);
    subplot(2,1,1)
    plot(ref_cal,Y0,'.')
    title('Calibration')
    grid on
    xlabel('Observed')
    ylabel('Predicted')
    subplot(2,1,2)
    plot(ref_val,Y1,'.')
    title('Prediction')
    grid on
    xlabel('Observed')
    ylabel('Predicted')
end

% OUTPUT
pred.coeff.SEC=SEC;
pred.coeff.SEP=SEP;
pred.coeff.RMSEP=RMSEP;
pred.coeff.R2c=R2c;
pred.coeff.R2p=R2p;
pred.coeff.SDc=SDc;
pred.coeff.SDp=SDp;
pred.coeff.RSD=RSD;
pred.coeff.RPDc=RPDc;
pred.coeff.RPDp=RPDp;
pred.coeff.Bias=Bias;
pred.coeff.Intercept=Intercept;
pred.coeff.Slope=Slope;
[~, m]=size(T);
pred.coeff.Recap={'Ymin' min([ref_cal; ref_val]); 'Ymax' max([ref_cal; ref_val]);'VL' m; 'SEC' SEC; 'SEP' SEP; 'RMSEP' RMSEP; 'R²c' R2c; 'R²p' R2p; 'SDc' SDc; 'SDp' SDp; 'RSD' RSD; 'RPDc' RPDc; 'RPDp' RPDp; 'Bias' Bias; 'Intercept' Intercept; 'Slope' Slope};
pred.para.T=T;
pred.para.P=P;
pred.para.W=W;
pred.para.Wstar=Wstar;
pred.para.Q=Q;
pred.para.B=B;
pred.para.R2X=R2X;
pred.para.R2Y=R2Y;
pred.para.RMSE=RMSE;
pred.Pred.Ycal=[Y_calb Y0];
pred.Pred.Yval=[Y_valb Y1];

pred.para.nbcomp=m;
end