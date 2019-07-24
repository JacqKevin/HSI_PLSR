function CreateModel(M,dm,wl,Y,dy,Yn)
% Function to create a PLSR model to predict a destructive variable with a
% hyperpectral image.
% INPUT :   
%           M: Hyperspectral datacube (n*m*p)
%           dm: Associated depth (1*m)
%           wl: Associated wavelengths (1*p)
%           Y: Reference values
%           dy: Reference depth
%           Yn: Name of the variable
%
% This is the software from the papers:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Jacq, K., Perrette, Y., Fanget, B., Sabatier, P., Coquin, D., 
% Martinez-Lamas, R., Debret, M., Arnaud, F., 2019. High-resolution 
% prediction of organic matter concentration with hyperspectral imaging 
% on a sediment core. Sci. Total Environ. 663, 236–244. 
% https://doi.org/10.1016/j.scitotenv.2019.01.320

% Please cite our papers if you use our code for your research.

%% Input verification
if size(M,2)~=length(dm)
    error(['M and dm do not have the same size'])
end
if size(M,3)~=length(wl)
    error(['M and wl do not have the same size'])
end

if size(Y,1)~=length(dy)
    error('Y and dy do not have the same dimension')
end
if iscell(Yn)||isnumeric(Yn)
    if size(Y,2)~=length(Yn)
        error('Y and Yn do not have the same dimension')
    end
else
    if size(Y,2)>1
        error('Y and Yn do not have the same dimension')
    end
end

if size(M,3)==98||size(M,3)==144
    M=M(:,:,15:end-10);
    wl=wl(15:end-10);
end

%% GUI at the screen size
set(gcf,'Units','pixels')
size_pixels=get(gcf,'Position');
set(gcf,'Units','characters')
size_characters=get(gcf,'Position');
f=size_pixels(3:4)./size_characters(3:4);
Pix_SS = get(0,'screensize');
a=Pix_SS(:,3:4);
a_characters=a./f;
close all

%% Conversion in pseudo-absorbance
pseudoabs = questdlg('Do you want to convert you data in pseudo-absorbance ?','Pseudo-Absorbance','Yes','No','No');
if strcmp(pseudoabs,'Yes')
    S=reshape(M,size(M,1)*size(M,2),size(M,3));
    [a,b]=find(S==0);
    S(a,b)=eps;
    S=log(1./S);
    M=reshape(S,size(M,1),size(M,2),size(M,3));
end

%% Subsampling of the hyperspectral data at the reference resolution 
if size(dy,2)==1
    % Zone ref
    dataref = questdlg('How your thickness are calculated?','Reference thickness','Min','Center','Max','Max');
    datadref = inputdlg('What is the thickness:');
    datadref=str2double(datadref);
    
    dyb=zeros(2,length(dy));
    if strcmp(dataref,'Max')
        for i=1:length(dy)
            dyb(1,i)=dy(i)-2/3*datadref;
            dyb(2,i)=dy(i)-1/3*datadref;
        end
    else if strcmp(dataref,'Min')
            for i=1:length(dy)
                dyb(1,i)=dy(i)+1/3*datadref;
                dyb(2,i)=dy(i)+2/3*datadref;
            end
        else if strcmp(dataref,'Center')
                for i=1:length(dy)
                    dyb(1,i)=dy(i)-1/6*datadref;
                    dyb(2,i)=dy(i)+1/6*datadref;
                end
            end
        end
    end
    
else
    for i=1:size(dy,1)
        dd=dy(i,2)-dy(i,1);
        dyb(1,i)=dy(i,1)+1/3*dd;
        dyb(2,i)=dy(i,1)+2/3*dd;
    end
end

datapretreat = questdlg('Do you want to preprocess your data?','Preprocessing','No','center','autoscaling','No');
gaplsdo='No';
% Depth registration
iter=1;
diy=zeros(2,length(dy));
for i=1:length(dy)
    for j=1:2
        if i==54
            a=1;
        end
        [~, a]=find(abs(dm-dyb(j,i))==min(abs(dm-dyb(j,i))));
        if abs(dm(a(1))-dyb(j,i))<0.5
            dis(j,iter)=a(1);
            Ye(iter,:)=Y(i,:);
        else
            diy(j)=1;
        end
    end
    if abs(dm(a(1))-dyb(j,i))<0.5
        iter=iter+1;
    end
end
if dis(2,end)==0
    dis=dis(:,1:end-1);
    Ye=Ye(1:end-1,:);
end

% Keep the central third 
Mi=M(round(1/3*size(M,1):2/3*size(M,1)),:,:);
Mia=M(round(1/3*size(M,1):2/3*size(M,1)),:,:);

% Preprocessing
Si=reshape(Mi,size(Mi,1)*size(Mi,2),size(Mi,3));
if strcmp(datapretreat,'center')
    [Si,para1,para2]=pretreat(Si,'center');
else if strcmp(datapretreat,'autoscaling')
        [Si,para1,para2]=pretreat(Si,'autoscaling');
    end
end
Mi=reshape(Si,size(Mi,1),size(Mi,2),size(Mi,3));

% Boostrap:
Sw=[];

disx=dis(2,:)-dis(1,:);
if sum(disx==0)>0
    a=disx==0;
    disx=disx(1,a==0);
    Ye=Ye(a==0,:);
end
if sum(disx<0)>0
    [a,b]=find(disx<0);
    for i=1:length(b)
        disxi= dis(1,b(i));
        dis(1,b(i))=dis(2,b(i));
        dis(2,b(i))=disxi;
    end
end
% coordonnée y représente le nombre de colonne du plan :
disy=ones(1,size(Ye,1))*size(Mi,1);

%% PLS creation
if size(Y,2)>1
    plst='PLS2';
else
    plst='PLS1';
end
% cvt = questdlg('Avec quelle algorithme de cross-validation voulez vous travailler ?','Cross-Validation','K-Fold','HoldOut','LeaveOut','K-Fold');
cvt='LeaveOut';
% Sub-sampling
ech = questdlg('How would you like to subsample?','Subsampling','Bootstrapping','Mean','Median','Bootstrapping');

% Wavelengths selection
chxwl = questdlg('How do you choose the number of wavelengths?','Wavelength selection','Manual','Automatic','Automatic');

% Subsampling methods
if strcmp(ech,'Bootstrapping')
    nbiterboostrap=100;
    for i = 1 : nbiterboostrap
        tiragex=ceil(rand(1,size(Ye,1)).*disx);
        tiragey=ceil(rand(1,size(Ye,1)).*disy);
        
        for j = 1:size(Ye,1)
            Sw(i,j,:)=Mi(tiragey(j),dis(j)-1+tiragex(j),:);
            Sboot(i,j,:)=Mia(tiragey(j),dis(j)-1+tiragex(j),:);
        end
    end
else if strcmp(ech,'Mean')
        for j = 1:size(Ye,1)
            Mii=[];
            Mii=Mi(:,dis(1,j):dis(2,j),:);
            Si=reshape(Mii,[size(Mii,1)*size(Mii,2) size(Mii,3)]);
            % Remove outliers
            [ACP,outliers] = ACPauto(Si, 0, wl);
            Si=Si(outliers(:,7)==0,:);
            % Remove outliers
            [ACP,outliers] = ACPauto(Si, 0, wl);
            Si=Si(outliers(:,7)==0,:);
            Sw(1,j,:)=mean(Si);
            Yee(j,:)=Ye(j,:);
        end
        if size(Ye,1)<100
            nbiterboostrap=size(Ye,1);
        else
            nbiterboostrap=100;
        end
        for k=2:nbiterboostrap
            Sw(k,:,:)=Sw(1,:,:);
        end
    else if strcmp(ech,'Median')
            for j = 1:size(Ye,1)
                Mii=[];
                Mii=Mi(:,dis(1,j):dis(2,j),:);
                Si=reshape(Mii,[size(Mii,1)*size(Mii,2) size(Mii,3)]);
                % Remove outliers
                [ACP,outliers] = ACPauto(Si, 0, wl);
                Si=Si(outliers(:,7)==0,:);
                % Remove outliers
                [ACP,outliers] = ACPauto(Si, 0, wl);
                Si=Si(outliers(:,7)==0,:);
                Sw(1,j,:)=median(Si);
                Yee(j,:)=Ye(j,:);
            end
            if size(Ye,1)<100
                nbiterboostrap=size(Ye,1);
            else
                nbiterboostrap=100;
            end
            for k=2:nbiterboostrap
                Sw(k,:,:)=Sw(1,:,:);
            end
        end
    end
end
Sboot=Sw;
Num=(1:size(Ye,1))';

% Saving subsampling data
DataPLS.S=Sboot;
DataPLS.Y=Ye;
DataPLS.Num=Num;
DataPLS.Pret=datapretreat;
if strcmp(datapretreat,'center')||strcmp(datapretreat,'autoscaling')
    DataPLS.Pret1=para1;
    DataPLS.Pret2=para2;
else
    DataPLS.Pret1=[1 1];
    DataPLS.Pret2=[1 1];
end
assignin('base', 'DataPLS', DataPLS);

%% Histogram
if size(Ye,2)<12
    fig = figure;
    fig.PaperPositionMode = 'auto';
    fig.InvertHardcopy = 'off';
    set(fig,'units','characters', 'position',[0 0 a_characters(1) a_characters(2)], 'name', 'Histogram');
    na=floor(sqrt(size(Ye,2)));
    nb=ceil(sqrt(size(Ye,2)));
    if na*nb<size(Ye,2)
        na=nb;
    end
    for i=1:size(Ye,2)
        subplot(na,nb,i);hist(Ye(:,i),100);grid on
    end
else
    nby=1;
    nbiter=floor(size(Ye,2)/12);
    for i=1:nbiter
        fig = figure;
        fig.PaperPositionMode = 'auto';
        fig.InvertHardcopy = 'off';
        set(fig,'units','characters', 'position',[0 0 a_characters(1) a_characters(2)], 'name', 'Histogram');
        for j=1:12
            subplot(3,4,j);hist(Ye(:,nby),100);grid on;title(Yn(nby))
            nby=nby+1;
        end
    end
    na=floor(sqrt(size(Ye,2)-nby+1));
    nb=ceil(sqrt(size(Ye,2)-nby+1));
    if na*nb<size(Ye,2)-nby+1
        na=nb;
    end
    fig = figure;
    fig.PaperPositionMode = 'auto';
    fig.InvertHardcopy = 'off';
    set(fig,'units','characters', 'position',[0 0 a_characters(1) a_characters(2)], 'name', 'Histogram');
    for i=1:size(Ye,2)-nby+1
        subplot(na,nb,i);hist(Ye(:,nby),100);grid on;title(Yn(nby))
        nby=nby+1;
    end
end

%% Cros-Validation
if strcmp(cvt,'K-Fold')
    cv=zeros(1,3);
    cv(1,1)=2;
    prompt = {'Number of group :','Number of iteration :'};
    dlg_title = 'Cross-validation';
    num_lines = 1;
    defaultans = {'10','100'};
    cvc = inputdlg(prompt,dlg_title,num_lines,defaultans);
    cv(1,2)=str2num(cvc{1,1});
    cv(1,3)=str2num(cvc{2,1});
else if strcmp(cvt,'HoldOut')
        cv=zeros(1,3);
        cv(1,1)=1;
        prompt = {'Ratio in validation (between 0 et 1) :','Number of iteration :'};
        dlg_title = 'Cross-validation';
        num_lines = 1;
        defaultans = {'0.5','100'};
        cvc = inputdlg(prompt,dlg_title,num_lines,defaultans);
        cv(1,2)=str2num(cvc{1,1});
        cv(1,3)=str2num(cvc{2,1});
    else if strcmp(cvt,'LeaveOut')
            cv=zeros(1,3);
            cv(1,1)=3;
        end
    end
end

if strcmp(plst,'PLS1')||strcmp(plst,'PLS1 et PLS2')
    % Bootstrap
    R2cal=zeros(size(Ye,2),nbiterboostrap,14);
    R2cv=zeros(size(Ye,2),nbiterboostrap,14);
    R2val=zeros(size(Ye,2),nbiterboostrap,14);
    Nbvl=zeros(size(Ye,2),nbiterboostrap,14);
    
    % Verify that there is several different values
    if size(Ye,1)/6>50
        lim=20;
    else
        lim=0.5*size(Ye,1)/6;
    end
    
    chxpret = questdlg('Which preprocessing would you like to do?','Preprocessing','Raw','All','Raw');
    
    if strcmp(chxpret,'Raw')
        h = waitbar(0,'PLS-CV-boostrap');
        iter=1;
        for k=1:nbiterboostrap
            % PLS-CV
            for i=1:size(Ye,2)
                if sum(abs(Ye(:,i)))>0&&length(unique(Ye(:,i)))>lim
                    pred=PLScv(squeeze(Sw(nbiterboostrap,:,:)),Ye(:,i),Num,wl,cv);
                    PRED{k,i}=pred;
                    
                    R2cal(i,k,:)=cell2mat(PRED{k,i}.recap  (2:end,4));
                    R2cv(i,k,:)=cell2mat(PRED{k,i}.recap  (2:end,7));
                    R2val(i,k,:)=cell2mat(PRED{k,i}.recap  (2:end,11));
                    Nbvl(i,k,:)=cell2mat(PRED{k,i}.recap  (2:end,2));
                end
                h = waitbar(iter/(nbiterboostrap*size(Ye,2)));
                iter=iter+1;
            end
        end
        close(h)
        
        % Find the optimal preprocessing
        R2g=R2cal+R2cv+R2val;
        
        R2calopt=zeros(size(Ye,2),nbiterboostrap);
        R2cvopt=zeros(size(Ye,2),nbiterboostrap);
        R2valopt=zeros(size(Ye,2),nbiterboostrap);
        Nbvlopt=zeros(size(Ye,2),nbiterboostrap);
        Pret={'Raw' 'Detrend' 'MSC' 'SNV' 'SNVD' 'D1' 'D2' 'MSCD1' 'SNVD1' 'SNVDD1' 'MSCD2' 'SNVD2' 'SNVDD2' 'CR'};
        
        % Choose the optimal preprocessing
        for k=1:nbiterboostrap
            for i=1:size(Ye,2)
                [opt(k),~]=find(squeeze(R2g(i,k,:))==max(squeeze(R2g(i,k,:))));
                Pretoptc(i,k)=opt(k);
            end
        end
        
        % Optimal preprocessing
        Pretoptg=mode(Pretoptc);
    else
        Pretoptg=1;
        h = waitbar(0,'PLS-CV-boostrap');
        iter=1;
        for k=1:nbiterboostrap
            % PLS-CV
            for i=1:size(Ye,2)
                if sum(abs(Ye(:,i)))>0&&length(unique(Ye(:,i)))>lim
                    pred=PLScv(squeeze(Sw(nbiterboostrap,:,:)),Ye(:,i),Num,wl,cv,Pretoptg);
                    PRED{k,i}=pred;
                    
                    R2cal(i,k,:)=pred.brut.coeff.R2c;
                    R2cv(i,k,:)=pred.brut.cv.r2cv;
                    R2val(i,k,:)=pred.brut.coeff.R2p;
                    Nbvl(i,k,:)=pred.brut.para.nbcomp;
                end
                h = waitbar(iter/(nbiterboostrap*size(Ye,2)));
                iter=iter+1;
            end
        end
        close(h)
        
        % Find the optimal preprocessing
        R2g=R2cal+R2cv+R2val;
        
        R2calopt=zeros(size(Ye,2),nbiterboostrap);
        R2cvopt=zeros(size(Ye,2),nbiterboostrap);
        R2valopt=zeros(size(Ye,2),nbiterboostrap);
        Nbvlopt=zeros(size(Ye,2),nbiterboostrap);
        Pret={'Raw' 'Detrend' 'MSC' 'SNV' 'SNVD' 'D1' 'D2' 'MSCD1' 'SNVD1' 'SNVDD1' 'MSCD2' 'SNVD2' 'SNVDD2' 'CR'};
    end
    
    for i=1:size(Ye,2)
        R2calopt(i,:)=squeeze(R2cal(i,:,Pretoptg));
        R2cvopt(i,:)=squeeze(R2cv(i,:,Pretoptg));
        R2valopt(i,:)=squeeze(R2val(i,:,Pretoptg));
        Nbvlopt(i,:)=squeeze(Nbvl(i,:,Pretoptg));
        Pretopt{i}=Pret(Pretoptg);
    end
    % Saving
    plsPred.PRED=PRED;
    plsPred.R2cal=R2cal;
    plsPred.R2cv=R2cv;
    plsPred.R2val=R2val;
    plsPred.Nbvl=Nbvl;
    assignin('base', 'plsPredPLS1boostrap', plsPred);
    
    plsPredOpt.predopt=PRED;
    plsPredOpt.R2calopt=R2calopt;
    plsPredOpt.R2cvopt=R2cvopt;
    plsPredOpt.R2valopt=R2valopt;
    plsPredOpt.Nbvlopt=Nbvlopt;
    plsPredOpt.Pretopt=Pretopt;
    assignin('base', 'plsPredOptPLS1', plsPredOpt);
    
    % Display
    for i=1:size(Ye,2)
        fig = figure;
        fig.PaperPositionMode = 'auto';
        fig.InvertHardcopy = 'off';
        set(fig,'units','characters', 'position',[0 0 a_characters(1) a_characters(2)], 'name', 'Results');
        subplot(221),hist(R2calopt(i,:),100),grid on,title('R² cal opt')
        subplot(222),hist(R2cvopt(i,:),100),grid on,title('R² cv opt')
        subplot(223),hist(R2valopt(i,:),100),grid on,title('R² val opt')
        subplot(224),hist(Nbvlopt(i,:),100),grid on,title('Nv VL opt')
    end
    
    %% Optimal bootstrapp set
    R2optg=R2calopt+2*R2cvopt+2*R2valopt;
    [optsortcal,optcal] = sort(R2calopt,'descend');
    [optsortcv,optcv] = sort(R2cvopt,'descend');
    [optsortval,optval] = sort(R2valopt,'descend');
    [optsort,opt] = sort(R2optg,'descend');
    %     [~, opt]=find(R2optg==max(R2optg));
    
    [~,idxset]=find(optsort>=optsort(1)-0.15);
    if length(idxset)>10
        idxset=idxset(1:5);
    end
    
    % Saving the optimal complete model
    if strcmp(pseudoabs,'Yes')
        Modelwsv.Type='PseudoAbsorbance';
    else
        Modelwsv.Type='Reflectance';
    end
    Modelwsv.label=Yn;
    Modelwsv.SelectWave='None';
    Modelwsv.cv=cvt;
    Modelwsv.Scal{1,1}=PRED{opt(1),1}.Scal;
    Modelwsv.Ycal=PRED{opt(1),1}.Ycal;
    Modelwsv.Sval{1,1}=PRED{opt(1),1}.Sval;
    Modelwsv.Yval=PRED{opt(1),1}.Yval;
    Pret={'Raw';'Detrend';'MSC';'SNV';'SNVD';'D1';'D2';'MSC+D1';'SNV+D1';'SNVD+D1';'MSC+D2';'SNV+D2';'SNVD+D2';'CR'};
    Modelwsv.Pret=Pret(Pretoptg);
    Modelwsv.Pret1=DataPLS.Pret;
    Modelwsv.Pret1a=[DataPLS.Pret1; DataPLS.Pret2];
    Modelwsv.wl=wl;
    if Pretoptg==1
        Modelwsv.B=PRED{opt(1),1}.brut.para.B;
        Modelwsv.R2cal=PRED{opt(1),1}.brut.coeff.R2c;
        Modelwsv.R2cv=PRED{opt(1),1}.brut.cv.r2cv;
        Modelwsv.R2val=PRED{opt(1),1}.brut.coeff.R2p;
        Modelwsv.SEC=PRED{opt(1),1}.brut.coeff.SEC;
        Modelwsv.SECV=PRED{opt(1),1}.brut.cv.SECV;
        Modelwsv.SEP=PRED{opt(1),1}.brut.coeff.SEP;
        Modelwsv.RMSEP=PRED{opt(1),1}.brut.coeff.RMSEP;
        Modelwsv.model=PRED{opt(1),1}.brut ;
    else if Pretoptg==2
            Modelwsv.B=PRED{opt(1),1}.detrend.para.B;
            Modelwsv.R2cal=PRED{opt(1),1}.detrend.coeff.R2c;
            Modelwsv.R2cv=PRED{opt(1),1}.detrend.cv.r2cv;
            Modelwsv.R2val=PRED{opt(1),1}.detrend.coeff.R2p;
            Modelwsv.SEC=PRED{opt(1),1}.detrend.coeff.SEC;
            Modelwsv.SECV=PRED{opt(1),1}.detrend.cv.SECV;
            Modelwsv.SEP=PRED{opt(1),1}.detrend.coeff.SEP;
            Modelwsv.RMSEP=PRED{opt(1),1}.detrend.coeff.RMSEP;
            Modelwsv.model=PRED{opt(1),1}.detrend ;
        else if Pretoptg==3
                Modelwsv.B=PRED{opt(1),1}.msc.para.B;
                Modelwsv.R2cal=PRED{opt(1),1}.msc.coeff.R2c;
                Modelwsv.R2cv=PRED{opt(1),1}.msc.cv.r2cv;
                Modelwsv.R2val=PRED{opt(1),1}.msc.coeff.R2p;
                Modelwsv.SEC=PRED{opt(1),1}.msc.coeff.SEC;
                Modelwsv.SECV=PRED{opt(1),1}.msc.cv.SECV;
                Modelwsv.SEP=PRED{opt(1),1}.msc.coeff.SEP;
                Modelwsv.RMSEP=PRED{opt(1),1}.msc.coeff.RMSEP;
                Modelwsv.model=PRED{opt(1),1}.msc ;
            else if Pretoptg==4
                    Modelwsv.B=PRED{opt(1),1}.snv.para.B;
                    Modelwsv.R2cal=PRED{opt(1),1}.snv.coeff.R2c;
                    Modelwsv.R2cv=PRED{opt(1),1}.snv.cv.r2cv;
                    Modelwsv.R2val=PRED{opt(1),1}.snv.coeff.R2p;
                    Modelwsv.SEC=PRED{opt(1),1}.snv.coeff.SEC;
                    Modelwsv.SECV=PRED{opt(1),1}.snv.cv.SECV;
                    Modelwsv.SEP=PRED{opt(1),1}.snv.coeff.SEP;
                    Modelwsv.RMSEP=PRED{opt(1),1}.snv.coeff.RMSEP;
                    Modelwsv.model=PRED{opt(1),1}.snv ;
                else if Pretoptg==5
                        Modelwsv.B=PRED{opt(1),1}.snvd.para.B;
                        Modelwsv.R2cal=PRED{opt(1),1}.snvd.coeff.R2c;
                        Modelwsv.R2cv=PRED{opt(1),1}.snvd.cv.r2cv;
                        Modelwsv.R2val=PRED{opt(1),1}.snvd.coeff.R2p;
                        Modelwsv.SEC=PRED{opt(1),1}.snvd.coeff.SEC;
                        Modelwsv.SECV=PRED{opt(1),1}.snvd.cv.SECV;
                        Modelwsv.SEP=PRED{opt(1),1}.snvd.coeff.SEP;
                        Modelwsv.RMSEP=PRED{opt(1),1}.snvd.coeff.RMSEP;
                        Modelwsv.model=PRED{opt(1),1}.snvd ;
                    else if Pretoptg==6
                            Modelwsv.B=PRED{opt(1),1}.d1.para.B;
                            Modelwsv.R2cal=PRED{opt(1),1}.d1.coeff.R2c;
                            Modelwsv.R2cv=PRED{opt(1),1}.d1.cv.r2cv;
                            Modelwsv.R2val=PRED{opt(1),1}.d1.coeff.R2p;
                            Modelwsv.SEC=PRED{opt(1),1}.d1.coeff.SEC;
                            Modelwsv.SECV=PRED{opt(1),1}.d1.cv.SECV;
                            Modelwsv.SEP=PRED{opt(1),1}.d1.coeff.SEP;
                            Modelwsv.RMSEP=PRED{opt(1),1}.d1.coeff.RMSEP;
                            Modelwsv.model=PRED{opt(1),1}.d1 ;
                        else if Pretoptg==7
                                Modelwsv.B=PRED{opt(1),1}.d2.para.B;
                                Modelwsv.R2cal=PRED{opt(1),1}.d2.coeff.R2c;
                                Modelwsv.R2cv=PRED{opt(1),1}.d2.cv.r2cv;
                                Modelwsv.R2val=PRED{opt(1),1}.d2.coeff.R2p;
                                Modelwsv.SEC=PRED{opt(1),1}.d2.coeff.SEC;
                                Modelwsv.SECV=PRED{opt(1),1}.d2.cv.SECV;
                                Modelwsv.SEP=PRED{opt(1),1}.d2.coeff.SEP;
                                Modelwsv.RMSEP=PRED{opt(1),1}.d2.coeff.RMSEP;
                                Modelwsv.model=PRED{opt(1),1}.d2 ;
                            else if Pretoptg==8
                                    Modelwsv.B=PRED{opt(1),1}.mscd1.para.B;
                                    Modelwsv.R2cal=PRED{opt(1),1}.mscd1.coeff.R2c;
                                    Modelwsv.R2cv=PRED{opt(1),1}.mscd1.cv.r2cv;
                                    Modelwsv.R2val=PRED{opt(1),1}.mscd1.coeff.R2p;
                                    Modelwsv.SEC=PRED{opt(1),1}.mscd1.coeff.SEC;
                                    Modelwsv.SECV=PRED{opt(1),1}.mscd1.cv.SECV;
                                    Modelwsv.SEP=PRED{opt(1),1}.mscd1.coeff.SEP;
                                    Modelwsv.RMSEP=PRED{opt(1),1}.mscd1.coeff.RMSEP;
                                    Modelwsv.model=PRED{opt(1),1}.mscd1 ;
                                else if Pretoptg==9
                                        Modelwsv.B=PRED{opt(1),1}.mscd2.para.B;
                                        Modelwsv.R2cal=PRED{opt(1),1}.mscd2.coeff.R2c;
                                        Modelwsv.R2cv=PRED{opt(1),1}.mscd2.cv.r2cv;
                                        Modelwsv.R2val=PRED{opt(1),1}.mscd2.coeff.R2p;
                                        Modelwsv.SEC=PRED{opt(1),1}.mscd2.coeff.SEC;
                                        Modelwsv.SECV=PRED{opt(1),1}.mscd2.cv.SECV;
                                        Modelwsv.SEP=PRED{opt(1),1}.mscd2.coeff.SEP;
                                        Modelwsv.RMSEP=PRED{opt(1),1}.mscd2.coeff.RMSEP;
                                        Modelwsv.model=PRED{opt(1),1}.mscd2 ;
                                    else if Pretoptg==10
                                            Modelwsv.B=PRED{opt(1),1}.snvd1.para.B;
                                            Modelwsv.R2cal=PRED{opt(1),1}.snvd1.coeff.R2c;
                                            Modelwsv.R2cv=PRED{opt(1),1}.snvd1.cv.r2cv;
                                            Modelwsv.R2val=PRED{opt(1),1}.snvd1.coeff.R2p;
                                            Modelwsv.SEC=PRED{opt(1),1}.snvd1.coeff.SEC;
                                            Modelwsv.SECV=PRED{opt(1),1}.snvd1.cv.SECV;
                                            Modelwsv.SEP=PRED{opt(1),1}.snvd1.coeff.SEP;
                                            Modelwsv.RMSEP=PRED{opt(1),1}.snvd1.coeff.RMSEP;
                                            Modelwsv.model=PRED{opt(1),1}.snvd1 ;
                                        else if Pretoptg==11
                                                Modelwsv.B=PRED{opt(1),1}.snvd2.para.B;
                                                Modelwsv.R2cal=PRED{opt(1),1}.snvd2.coeff.R2c;
                                                Modelwsv.R2cv=PRED{opt(1),1}.snvd2.cv.r2cv;
                                                Modelwsv.R2val=PRED{opt(1),1}.snvd2.coeff.R2p;
                                                Modelwsv.SEC=PRED{opt(1),1}.snvd2.coeff.SEC;
                                                Modelwsv.SECV=PRED{opt(1),1}.snvd2.cv.SECV;
                                                Modelwsv.SEP=PRED{opt(1),1}.snvd2.coeff.SEP;
                                                Modelwsv.RMSEP=PRED{opt(1),1}.snvd2.coeff.RMSEP;
                                                Modelwsv.model=PRED{opt(1),1}.snvd2 ;
                                            else if Pretoptg==12
                                                    Modelwsv.B=PRED{opt(1),1}.snvdd1.para.B;
                                                    Modelwsv.R2cal=PRED{opt(1),1}.snvdd1.coeff.R2c;
                                                    Modelwsv.R2cv=PRED{opt(1),1}.snvdd1.cv.r2cv;
                                                    Modelwsv.R2val=PRED{opt(1),1}.snvdd1.coeff.R2p;
                                                    Modelwsv.SEC=PRED{opt(1),1}.snvdd1.coeff.SEC;
                                                    Modelwsv.SECV=PRED{opt(1),1}.snvdd1.cv.SECV;
                                                    Modelwsv.SEP=PRED{opt(1),1}.snvdd1.coeff.SEP;
                                                    Modelwsv.RMSEP=PRED{opt(1),1}.snvdd1.coeff.RMSEP;
                                                    Modelwsv.model=PRED{opt(1),1}.snvdd1 ;
                                                else if Pretoptg==13
                                                        Modelwsv.B=PRED{opt(1),1}.snvdd2.para.B;
                                                        Modelwsv.R2cal=PRED{opt(1),1}.snvdd2.coeff.R2c;
                                                        Modelwsv.R2cv=PRED{opt(1),1}.snvdd2.cv.r2cv;
                                                        Modelwsv.R2val=PRED{opt(1),1}.snvdd2.coeff.R2p;
                                                        Modelwsv.SEC=PRED{opt(1),1}.snvdd2.coeff.SEC;
                                                        Modelwsv.SECV=PRED{opt(1),1}.snvdd2.cv.SECV;
                                                        Modelwsv.SEP=PRED{opt(1),1}.snvdd2.coeff.SEP;
                                                        Modelwsv.RMSEP=PRED{opt(1),1}.snvdd2.coeff.RMSEP;
                                                        Modelwsv.model=PRED{opt(1),1}.snvdd2 ;
                                                    else if Pretoptg==14
                                                            Modelwsv.B=PRED{opt(1),1}.cr.para.B;
                                                            Modelwsv.R2cal=PRED{opt(1),1}.cr.coeff.R2c;
                                                            Modelwsv.R2cv=PRED{opt(1),1}.cr.cv.r2cv;
                                                            Modelwsv.R2val=PRED{opt(1),1}.cr.coeff.R2p;
                                                            Modelwsv.SEC=PRED{opt(1),1}.cr.coeff.SEC;
                                                            Modelwsv.SECV=PRED{opt(1),1}.cr.cv.SECV;
                                                            Modelwsv.SEP=PRED{opt(1),1}.cr.coeff.SEP;
                                                            Modelwsv.RMSEP=PRED{opt(1),1}.cr.coeff.RMSEP;
                                                            Modelwsv.model=PRED{opt(1),1}.cr ;
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    assignin('base','Modelwvs',Modelwsv);
    
    for e=1:length(idxset)
        %% Wavelength selection
        if Pretoptg==1
            Swsv=squeeze(Sw(opt(idxset(e)),:,:));
        else
            Swsv = AllPret(squeeze(Sw(opt(idxset(e)),:,:)),wl,Pretoptg-1);
        end
        
        for z=1:size(Ye,2)
            SV{e}=AllSelectVar(Swsv,Ye(:,z),wl,mode(Nbvlopt));
            
            if strcmp(gaplsdo,'Yes')
                % GA-PLS
                [b,~,sel]=gaplssp([Swsv Ye(:,z)],500);
                SV{e}.wGAPLS=sel;
                SV{e}.idxGAPLS=b;
            end
            
            SV{e}.setopt=opt(idxset(e));
        end
        
        % VIP
        %     predSV.predVIP=PLScv(squeeze(Sw(opt,:,:)),Ye(:,i),Num,wl(SV.idxVIPopt),cv,Pretoptg,SV.idxVIPopt);
        %     predSV.VIPwl=wl(SV.idxVIPopt);
        
        %         % TP
        %         for k=10:2*length(wl)/5
        %             predTP{k}=PLScv(Swsv,Ye(:,i),Num,wl(SV{e}.idxTP(1:k)),cv,1,SV{e}.idxTP(1:k));
        %             TPcal(k)=predTP{1, k}.brut.coeff.R2c;
        %             TPval(k)=predTP{1, k}.brut.coeff.R2p;
        %             TPcv(k)=predTP{1, k}.brut.cv.r2cv;
        %         end
        %         TPg=TPcal+TPcv+TPval;
        %         [~,idx]=find(TPg(10:end)==max(TPg(10:end)));
        %         predSV{e}.predTP=predTP{idx+9};
        %         predSV{e}.TPwl=wl(SV{e}.idxTP(1:idx));
        %
        %         % UVE
        %         for k=10:2*length(wl)/5
        %             predUVE{k}=PLScv(Swsv,Ye(:,i),Num,wl(SV{e}.idxUVE(1:k)),cv,1,SV{e}.idxUVE(1:k));
        %             UVEcal(k)=predUVE{1, k}.brut.coeff.R2c;
        %             UVEval(k)=predUVE{1, k}.brut.coeff.R2p;
        %             UVEcv(k)=predUVE{1, k}.brut.cv.r2cv;
        %         end
        %         UVEg=UVEcal+UVEcv+UVEval;
        %         [~,idx]=find(UVEg(10:end)==max(UVEg(10:end)));
        %         predSV{e}.predUVE=predUVE{idx+9};
        %         predSV{e}.UVEwl=wl(SV{e}.idxUVE(1:idx));
        
        %     % CARS
        %     predSV.predCARS=PLScv(squeeze(Sw(opt,:,:)),Ye(:,i),Num,wl(SV.idxCARSopt),cv,Pretoptg,SV.idxCARSopt);
        %     predSV.CARSwl=wl(SV.idxCARSopt);
        
        % RF
        for k=10:4*length(wl)/5
            predRF{k}=PLScv(Swsv,Ye(:,i),Num,wl(SV{e}.idxRF(1:k)),cv,1,SV{e}.idxRF(1:k));
            RFcal(k)=predRF{1, k}.brut.coeff.R2c;
            RFval(k)=predRF{1, k}.brut.coeff.R2p;
            RFcv(k)=predRF{1, k}.brut.cv.r2cv;
        end
        RFg=RFcal+RFcv+RFval;
        [~,idx]=find(RFg(10:end)==max(RFg(10:end)));
        predSV{e}.predRF=predRF{idx+9};
        predSV{e}.RFwl=wl(SV{e}.idxRF(1:idx));
        
        %     % IRF
        %     for k=10:2*length(wl)/5
        %         predIRF{k}=PLScv(squeeze(Sw(opt,:,:)),Ye(:,i),Num,wl(SV.idxiRF(1:k)),cv,Pretoptg,SV.idxiRF(1:k));
        %         IRFcal(k)=predIRF{1, k}.brut.coeff.R2c;
        %         IRFval(k)=predIRF{1, k}.brut.coeff.R2p;
        %     end
        %     IRFg=IRFcal+IRFval;
        %     [~,idx]=find(IRFg==max(IRFg));
        %     predSV.predIRF=predIRF{idx};
        %     predSV.IRFwl=wl(SV.idxiRF(1:idx));
        
        % IRIV % BUG SUR LA PLS
        %     predSV.predIRIS=PLScv(squeeze(Swsv(:,SV.idxIRIVopt)),Ye(:,i),Num,wl(SV.idxIRIVopt),cv,1);
        
        % VCN
        predSV{e}.predVCN=PLScv(Swsv,Ye(:,i),Num,wl(SV{e}.idxVCNopt),cv,1,SV{e}.idxVCNopt);
        predSV{e}.VCNwl=wl(SV{e}.idxVCNopt);
        
        if strcmp(gaplsdo,'Yes')
            % GA-PLS
            for k=5:2*length(wl)/5
                predGAPLS{k}=PLScv(Swsv,Ye(:,i),Num,wl(SV{e}.idxGAPLS(1:k)),cv,1,SV{e}.idxGAPLS(1:k));
                GAPLScal(k)=predGAPLS{1, k}.brut.coeff.R2c;
                GAPLSval(k)=predGAPLS{1, k}.brut.coeff.R2p;
                GAPLScv(k)=predGAPLS{1, k}.brut.cv.r2cv;
            end
            GAPLSg=GAPLScal+GAPLScv+GAPLSval;
            [~,idx]=find(GAPLSg(5:end)==max(GAPLSg(5:end)));
            predSV{e}.predGAPLS=predGAPLS{idx+4};
            predSV{e}.GAPLSwl=wl(SV{e}.idxGAPLS(1:idx));
        end
        
    end
    % Concatenate all random set selected wavelength
    settest=opt(1);
    TPwlall=[];UVEwlall=[];RFwlall=[];VCNwlall=[];GAPLSwlall=[];
    for ibstr=1:e
        %         TPwlall=[TPwlall predSV{1, ibstr}.TPwl];
        %         UVEwlall=[UVEwlall predSV{1, ibstr}.UVEwl];
        RFwlall=[RFwlall predSV{1, ibstr}.RFwl];
        VCNwlall=[VCNwlall predSV{1, ibstr}.VCNwl];
        if strcmp(gaplsdo,'Yes')
            GAPLSwlall=[GAPLSwlall predSV{1, ibstr}.GAPLSwl];
        end
    end
    %     % TP
    %     [CTP,ia,ic] = unique(sort(TPwlall));
    %     ia=[ia; length(TPwlall)];
    %     ianbTP=ia(2:end)-ia(1:end-1);
    %     ianbTPmax=max(ianbTP);
    %     ianbTPmin=min(ianbTP);
    %     idxwlTP=zeros(1,length(CTP));
    %     % idx des longueur d'onde
    %     for jwl=1:length(CTP)
    %         idxwlTP(jwl)=find(abs(CTP(jwl)-wl)==min(abs(CTP(jwl)-wl)));
    %     end
    %     % tris par importance et aléatoirement
    %     idxwlTPsort=[];
    %     for iwl=ianbTPmax:-1:ianbTPmin
    %         idwwlTPiwl=idxwlTP(ianbTP==iwl);
    %         aleawl=rand(1,length(idwwlTPiwl));
    %         [~,idxalea]=sort(aleawl);
    %         idxwlTPsort=[idxwlTPsort idwwlTPiwl(idxalea)];
    %     end
    %     % test PLS
    %     if Pretoptg==1
    %         Swsv=squeeze(Sw(settest,:,:));
    %     else
    %         Swsv = AllPret(squeeze(Sw(settest,:,:)),wl,Pretoptg-1);
    %     end
    %     if length(idxwlTPsort)<2*length(wl)/5
    %         lim=length(idxwlTPsort);
    %     else
    %         lim=2*length(wl)/5;
    %     end
    %     for k=5:lim
    %         predTP2{k}=PLScv(Swsv,Ye(:,i),Num,wl(idxwlTPsort(1:k)),cv,1,idxwlTPsort(1:k));
    %         TPcal2(k)=predTP2{1, k}.brut.coeff.R2c;
    %         TPval2(k)=predTP2{1, k}.brut.coeff.R2p;
    %         TPcv2(k)=predTP2{1, k}.brut.cv.r2cv;
    %     end
    %     TPg2=TPcal2+TPcv2+TPval2;
    %     if max(TPval2)>0.75
    %         [TPgs,TPgidx]=sort(TPval2,'descend');
    %         idx=min(TPgidx(TPgs>0.75));
    %     else if max(TPval2)>0.5
    %             [TPgs,TPgidx]=sort(TPval2,'descend');
    %             idx=min(TPgidx(TPgs>0.5));
    %         else
    %             [~,idx]=find(TPval2(5:end)==max(TPval2(5:end)));
    %             idx=idx+4;
    %         end
    %     end
    %     predSV2.predTP=predTP2{idx(1)};
    %     predSV2.TPwl=wl(idxwlTPsort(1:idx(1)));
    %     predSV2.TPidx=idxwlTPsort(1:idx(1));
    %
    %     % UVE
    %     [CUVE,ia,ic] = unique(sort(UVEwlall));
    %     ia=[ia; length(UVEwlall)];
    %     ianbUVE=ia(2:end)-ia(1:end-1);
    %     ianbUVEmax=max(ianbUVE);
    %     ianbUVEmin=min(ianbUVE);
    %     idxwlUVE=zeros(1,length(CUVE));
    %     % idx des longueur d'onde
    %     for jwl=1:length(CUVE)
    %         idxwlUVE(jwl)=find(abs(CUVE(jwl)-wl)==min(abs(CUVE(jwl)-wl)));
    %     end
    %     % tris par importance et aléatoirement
    %     idxwlUVEsort=[];
    %     for iwl=ianbUVEmax:-1:ianbUVEmin
    %         idwwlUVEiwl=idxwlUVE(ianbUVE==iwl);
    %         aleawl=rand(1,length(idwwlUVEiwl));
    %         [~,idxalea]=sort(aleawl);
    %         idxwlUVEsort=[idxwlUVEsort idwwlUVEiwl(idxalea)];
    %     end
    %     % test PLS
    %     if Pretoptg==1
    %         Swsv=squeeze(Sw(settest,:,:));
    %     else
    %         Swsv = AllPret(squeeze(Sw(settest,:,:)),wl,Pretoptg-1);
    %     end
    %     if length(idxwlUVEsort)<2*length(wl)/5
    %         lim=length(idxwlUVEsort);
    %     else
    %         lim=2*length(wl)/5;
    %     end
    %     for k=5:lim
    %         predUVE2{k}=PLScv(Swsv,Ye(:,i),Num,wl(idxwlUVEsort(1:k)),cv,1,idxwlUVEsort(1:k));
    %         UVEcal2(k)=predUVE2{1, k}.brut.coeff.R2c;
    %         UVEval2(k)=predUVE2{1, k}.brut.coeff.R2p;
    %         UVEcv2(k)=predUVE2{1, k}.brut.cv.r2cv;
    %     end
    %     UVEg2=UVEcal2+UVEcv2+UVEval2;
    %     if max(UVEval2)>0.75
    %         [UVEgs,UVEgidx]=sort(UVEval2,'descend');
    %         idx=min(UVEgidx(UVEgs>0.75));
    %     else if max(UVEval2)>0.5
    %             [UVEgs,UVEgidx]=sort(UVEval2,'descend');
    %             idx=min(UVEgidx(UVEgs>0.5));
    %         else
    %             [~,idx]=find(UVEval2(5:end)==max(UVEval2(5:end)));
    %             idx=idx+4;
    %         end
    %     end
    %     predSV2.predUVE=predUVE2{idx(1)};
    %     predSV2.UVEwl=wl(idxwlUVEsort(1:idx(1)));
    %     predSV2.UVEidx=idxwlUVEsort(1:idx(1));
    
    % RF
    [CRF,ia,ic] = unique(sort(RFwlall));
    ia=[ia; length(RFwlall)];
    ianbRF=ia(2:end)-ia(1:end-1);
    ianbRFmax=max(ianbRF);
    ianbRFmin=min(ianbRF);
    idxwlRF=zeros(1,length(CRF));
    % Indice of wavelengths
    for jwl=1:length(CRF)
        idxwlRF(jwl)=find(abs(CRF(jwl)-wl)==min(abs(CRF(jwl)-wl)));
    end
    % Sort wavelengths
    idxwlRFsort=[];
    for iwl=ianbRFmax:-1:ianbRFmin
        idwwlRFiwl=idxwlRF(ianbRF==iwl);
        aleawl=rand(1,length(idwwlRFiwl));
        [~,idxalea]=sort(aleawl);
        idxwlRFsort=[idxwlRFsort idwwlRFiwl(idxalea)];
    end
    % test PLS
    if Pretoptg==1
        Swsv=squeeze(Sw(settest,:,:));
    else
        Swsv = AllPret(squeeze(Sw(settest,:,:)),wl,Pretoptg-1);
    end
    if length(idxwlRFsort)<4*length(wl)/5
        lim=length(idxwlRFsort);
    else
        lim=4*length(wl)/5;
    end
    % zz=0;
    % while RFval2(idx)<0.7
    %     zz+1
    for k=5:lim
        predRF2{k}=PLScv(Swsv,Ye(:,i),Num,wl(idxwlRFsort(1:k)),cv,1,idxwlRFsort(1:k));
        RFcal2(k)=predRF2{1, k}.brut.coeff.R2c;
        RFval2(k)=predRF2{1, k}.brut.coeff.R2p;
        RFcv2(k)=predRF2{1, k}.brut.cv.r2cv;
        RFval2rmsep(k)=predRF2{1, k}.brut.coeff.RMSEP;
    end
    RFg2=RFcal2+RFcv2+RFval2;
    
    fig = figure;
    fig.PaperPositionMode = 'auto';
    fig.InvertHardcopy = 'off';
    set(fig,'units','characters', 'position',[0 0 a_characters(1) a_characters(2)], 'name', 'RF selected wavelengths');
    plot(CRF,ianbRF,'*')
    grid on
    title('RF selected wavelengths ')
    
    % Choose the optimal number of wavelengths
    if strcmp(chxwl,'Manual')
        fig = figure;
        fig.PaperPositionMode = 'auto';
        fig.InvertHardcopy = 'off';
        set(fig,'units','characters', 'position',[0 0 a_characters(1) a_characters(2)], 'name','RF wavelength sorts');
        subplot(211)
        plot(5:lim,RFval2(5:lim))
        hold on
        plot(5:lim,RFcal2(5:lim),'k--','linewidth',3)
        grid on
        title('RF selected wavelength R²val')
        legend({'Val' 'Cal'})
        ylim([-0.1 1])
        subplot(212)
        plot(5:lim,RFval2rmsep(5:lim))
        grid on
        title('RF selected wavelength RMSEP')
        ylim([min(RFval2rmsep(5:lim)) max(RFval2rmsep(RFval2rmsep<max(Ye)))])
        
        if max(RFval2)>0.9
            [RFgs,RFgidx]=sort(RFval2,'descend');
            idxit=min(RFgidx(RFgs>0.9));
        else if max(RFval2)>0.75
                [RFgs,RFgidx]=sort(RFval2,'descend');
                idxit=min(RFgidx(RFgs>0.75));
            else if max(RFval2)>0.5
                    [RFgs,RFgidx]=sort(RFval2,'descend');
                    idxit=min(RFgidx(RFgs>0.5));
                else
                    [~,idxit]=find(RFval2(5:end)==max(RFval2(5:end)));
                    idxit=idxit+4;
                end
            end
        end
        
        idxc = inputdlg('With how many wavelengths do you want to work?','Number of wavelengths?',[1 35],cellstr(num2str(idxit)));
        idx = str2num(idxc{1});
    else if strcmp(chxwl,'Automatic')
            if max(RFval2)>0.9
                [RFgs,RFgidx]=sort(RFval2,'descend');
                idx=min(RFgidx(RFgs>0.9));
            else if max(RFval2)>0.75
                    [RFgs,RFgidx]=sort(RFval2,'descend');
                    idx=min(RFgidx(RFgs>0.75));
                else if max(RFval2)>0.5
                        [RFgs,RFgidx]=sort(RFval2,'descend');
                        idx=min(RFgidx(RFgs>0.5));
                    else
                        [~,idx]=find(RFval2(5:end)==max(RFval2(5:end)));
                        idx=idx+4;
                    end
                end
            end
        end
    end
    % end
    predSV2.predRF=predRF2{idx(1)};
    predSV2.RFwl=wl(idxwlRFsort(1:idx(1)));
    predSV2.RFidx=idxwlRFsort(1:idx(1));
    
    % VCN
    [CVCN,ia,ic] = unique(sort(VCNwlall));
    ia=[ia; length(VCNwlall)];
    ianbVCN=ia(2:end)-ia(1:end-1);
    ianbVCNmax=max(ianbVCN);
    ianbVCNmin=min(ianbVCN);
    idxwlVCN=zeros(1,length(CVCN));
    % Indice of wavelengths
    for jwl=1:length(CVCN)
        idxwlVCN(jwl)=find(abs(CVCN(jwl)-wl)==min(abs(CVCN(jwl)-wl)));
    end
    % Sort wavelengths
    idxwlVCNsort=[];
    for iwl=ianbVCNmax:-1:ianbVCNmin
        idwwlVCNiwl=idxwlVCN(ianbVCN==iwl);
        aleawl=rand(1,length(idwwlVCNiwl));
        [~,idxalea]=sort(aleawl);
        idxwlVCNsort=[idxwlVCNsort idwwlVCNiwl(idxalea)];
    end
    % test PLS
    if Pretoptg==1
        Swsv=squeeze(Sw(settest,:,:));
    else
        Swsv = AllPret(squeeze(Sw(settest,:,:)),wl,Pretoptg-1);
    end
    if length(idxwlVCNsort)<2*length(wl)/5
        lim=length(idxwlVCNsort);
    else
        lim=2*length(wl)/5;
    end
    for k=5:lim
        predVCN2{k}=PLScv(Swsv,Ye(:,i),Num,wl(idxwlVCNsort(1:k)),cv,1,idxwlVCNsort(1:k));
        VCNcal2(k)=predVCN2{1, k}.brut.coeff.R2c;
        VCNval2(k)=predVCN2{1, k}.brut.coeff.R2p;
        VCNcv2(k)=predVCN2{1, k}.brut.cv.r2cv;
    end
    VCNg2=VCNcal2+VCNcv2+VCNval2;
    if max(VCNval2)>0.75
        [VCNgs,VCNgidx]=sort(VCNval2,'descend');
        idx=min(VCNgidx(VCNgs>0.75));
    else if max(VCNval2)>0.5
            [VCNgs,VCNgidx]=sort(VCNval2,'descend');
            idx=min(VCNgidx(VCNgs>0.5));
        else
            [~,idx]=find(VCNval2(5:end)==max(VCNval2(5:end)));
            idx=idx+4;
        end
    end
    predSV2.predVCN=predVCN2{idx(1)};
    predSV2.VCNwl=wl(idxwlVCNsort(1:idx(1)));
    predSV2.VCNidx=idxwlVCNsort(1:idx(1));
    
    if strcmp(gaplsdo,'Yes')
        % GAPLS
        [CGAPLS,ia,ic] = unique(sort(GAPLSwlall));
        ia=[ia; length(GAPLSwlall)];
        ianbGAPLS=ia(2:end)-ia(1:end-1);
        ianbGAPLSmax=max(ianbGAPLS);
        ianbGAPLSmin=min(ianbGAPLS);
        idxwlGAPLS=zeros(1,length(CGAPLS));
        % Indice of wavelengths
        for jwl=1:length(CGAPLS)
            idxwlGAPLS(jwl)=find(abs(CGAPLS(jwl)-wl)==min(abs(CGAPLS(jwl)-wl)));
        end
        % Sort wavelengths
        idxwlGAPLSsort=[];
        for iwl=ianbGAPLSmax:-1:ianbGAPLSmin
            idwwlGAPLSiwl=idxwlGAPLS(ianbGAPLS==iwl);
            aleawl=rand(1,length(idwwlGAPLSiwl));
            [~,idxalea]=sort(aleawl);
            idxwlGAPLSsort=[idxwlGAPLSsort idwwlGAPLSiwl(idxalea)];
        end
        % test PLS
        if Pretoptg==1
            Swsv=squeeze(Sw(settest,:,:));
        else
            Swsv = AllPret(squeeze(Sw(settest,:,:)),wl,Pretoptg-1);
        end
        if length(idxwlGAPLSsort)<2*length(wl)/5
            lim=length(idxwlGAPLSsort);
        else
            lim=2*length(wl)/5;
        end
        for k=5:lim
            predGAPLS2{k}=PLScv(Swsv,Ye(:,i),Num,wl(idxwlGAPLSsort(1:k)),cv,1,idxwlGAPLSsort(1:k));
            GAPLScal2(k)=predGAPLS2{1, k}.brut.coeff.R2c;
            GAPLSval2(k)=predGAPLS2{1, k}.brut.coeff.R2p;
            GAPLScv2(k)=predGAPLS2{1, k}.brut.cv.r2cv;
        end
        GAPLSg2=GAPLScal2+GAPLScv2+GAPLSval2;
        if max(GAPLSval2)>0.9
            [GAPLSgs,GAPLSgidx]=sort(GAPLSval2,'descend');
            idx=min(GAPLSgidx(GAPLSgs>0.9));
        else if max(GAPLSval2)>0.75
                [GAPLSgs,GAPLSgidx]=sort(GAPLSval2,'descend');
                idx=min(GAPLSgidx(GAPLSgs>0.75));
            else if max(GAPLSval2)>0.5
                    [GAPLSgs,GAPLSgidx]=sort(GAPLSval2,'descend');
                    idx=min(GAPLSgidx(GAPLSgs>0.5));
                else
                    [~,idx]=find(GAPLSval2(5:end)==max(GAPLSval2(5:end)));
                    idx=idx+4;
                end
            end
        end
        predSV2.predGAPLS=predGAPLS2{idx(1)};
        predSV2.GAPLSwl=wl(idxwlGAPLSsort(1:idx(1)));
        predSV2.GAPLSidx=idxwlGAPLSsort(1:idx(1));
    end
    
    assignin('base','resSV',SV);
    assignin('base','predSV',predSV);
    assignin('base','predSV2',predSV2);
    % Summary table :
    Recap{1,1}='Method';
    %     Recap{1,2}='VIP';
    %     Recap{1,3}='TP';
    %     Recap{1,4}='UVE';
    %     Recap{1,5}='CARS';
    Recap{1,6}='RF';
    %     Recap{1,7}='iRF';
    Recap{1,8}='VCN';
    if strcmp(gaplsdo,'Yes')
        Recap{1,9}='GA-PLS';
    end
    
    Recap{2,1}='VL';
    %     Recap{2,2}=predSV.predVIP.brut.para.nbcomp;
    %     Recap{2,3}=predSV2.predTP.brut.para.nbcomp;
    %     Recap{2,4}=predSV2.predUVE.brut.para.nbcomp;
    %     Recap{2,5}=predSV.predCARS.brut.para.nbcomp;
    Recap{2,6}=predSV2.predRF.brut.para.nbcomp;
    %     Recap{2,7}=predSV.predIRF.brut.para.nbcomp;
    Recap{2,8}=predSV2.predVCN.brut.para.nbcomp;
    if strcmp(gaplsdo,'Yes')
        Recap{2,9}=predSV2.predGAPLS.brut.para.nbcomp;
    end
    
    Recap{3,1}='NbWl';
    %     Recap{3,2}=size(predSV.predVIP.brut.para.P,1);
    %     Recap{3,3}=size(predSV2.predTP.brut.para.P,1);
    %     Recap{3,4}=size(predSV2.predUVE.brut.para.P,1);
    %     Recap{3,5}=size(predSV.predCARS.brut.para.P,1);
    Recap{3,6}=size(predSV2.predRF.brut.para.P,1);
    %     Recap{3,7}=size(predSV.predIRF.brut.para.P,1);
    Recap{3,8}=size(predSV2.predVCN.brut.para.P,1);
    if strcmp(gaplsdo,'Yes')
        Recap{3,9}=size(predSV2.predGAPLS.brut.para.P,1);
    end
    
    Recap{4,1}='R²cal';
    %     Recap{4,2}=predSV.predVIP.brut.coeff.R2c;
    %     Recap{4,3}=predSV2.predTP.brut.coeff.R2c;
    %     Recap{4,4}=predSV2.predUVE.brut.coeff.R2c;
    %     Recap{4,5}=predSV.predCARS.brut.coeff.R2c;
    Recap{4,6}=predSV2.predRF.brut.coeff.R2c;
    %     Recap{4,7}=predSV.predIRF.brut.coeff.R2c;
    Recap{4,8}=predSV2.predVCN.brut.coeff.R2c;
    if strcmp(gaplsdo,'Yes')
        Recap{4,9}=predSV2.predGAPLS.brut.coeff.R2c;
    end
    
    Recap{5,1}='R²cv';
    %     Recap{5,2}=predSV.predVIP.brut.coeff.R2p;
    %     Recap{5,3}=predSV2.predTP.brut.cv.r2cv;
    %     Recap{5,4}=predSV2.predUVE.brut.cv.r2cv;
    %     Recap{5,5}=predSV.predCARS.brut.coeff.R2p;
    Recap{5,6}=predSV2.predRF.brut.cv.r2cv;
    %     Recap{5,7}=predSV.predIRF.brut.coeff.R2p;
    Recap{5,8}=predSV2.predVCN.brut.cv.r2cv;
    if strcmp(gaplsdo,'Yes')
        Recap{5,9}=predSV2.predGAPLS.brut.cv.r2cv;
    end
    
    Recap{6,1}='R²val';
    %     Recap{6,2}=predSV.predVIP.brut.coeff.R2p;
    %     Recap{6,3}=predSV2.predTP.brut.coeff.R2p;
    %     Recap{6,4}=predSV2.predUVE.brut.coeff.R2p;
    %     Recap{6,5}=predSV.predCARS.brut.coeff.R2p;
    Recap{6,6}=predSV2.predRF.brut.coeff.R2p;
    %     Recap{6,7}=predSV.predIRF.brut.coeff.R2p;
    Recap{6,8}=predSV2.predVCN.brut.coeff.R2p;
    if strcmp(gaplsdo,'Yes')
        Recap{6,9}=predSV2.predGAPLS.brut.coeff.R2p;
    end
    
    Recap{7,1}='RMSEP';
    %     Recap{7,2}=predSV.predVIP.brut.coeff.RMSEP;
    %     Recap{7,3}=predSV2.predTP.brut.coeff.RMSEP;
    %     Recap{7,4}=predSV2.predUVE.brut.coeff.RMSEP;
    %     Recap{7,5}=predSV.predCARS.brut.coeff.RMSEP;
    Recap{7,6}=predSV2.predRF.brut.coeff.RMSEP;
    %     Recap{7,7}=predSV.predIRF.brut.coeff.RMSEP;
    Recap{7,8}=predSV2.predVCN.brut.coeff.RMSEP;
    if strcmp(gaplsdo,'Yes')
        Recap{7,9}=predSV2.predGAPLS.brut.coeff.RMSEP;
    end
    
    predSV2.Recap=Recap;
    assignin('base','predSV2',predSV2);
    
    % Affichage des zones influentes
    if strcmp(gaplsdo,'Yes')
        a=Recap{6,6}>0.3||Recap{6,8}>0.3||Recap{6,9}>0.3;%|Recap{5,5}>0.5||Recap{5,7}>0.5||Recap{6,3}>0.3||Recap{6,4}>0.3||
    else
        a=Recap{6,6}>0.3||Recap{6,8}>0.3;%%||Recap{6,9}>0.3%|Recap{5,5}>0.5||Recap{5,7}>0.5||Recap{6,3}>0.3||Recap{6,4}>0.3||
    end
    if a
        fig = figure;
        fig.PaperPositionMode = 'auto';
        fig.InvertHardcopy = 'off';
        set(fig,'units','characters', 'position',[0 0 a_characters(1) a_characters(2)], 'name', 'Influential areas');
        i=0;
        j=1;
        %         if Recap{5,2}>0.5
        %             plot(wl(SV.idxVIPopt),ones(1,length(SV.idxVIPopt))+i,'*','Markersize',10), hold on
        %             i=i+0.1;
        %             legvl{j}=strcat('VIP : ', num2str(round(Recap{5,2},2)), '(', num2str(Recap{3,2}), '), ', num2str(Recap{7,2}));
        %             j=j+1;
        %         end
        %         if Recap{6,3}>0.3
        %             plot(wl(predSV2.TPidx),ones(1,length(predSV2.TPidx))+i,'*','Markersize',10), hold on
        %             i=i+0.1;
        %             legvl{j}=strcat('TP : ', num2str(round(Recap{6,3},2)), '(', num2str(Recap{3,3}), '), ', num2str(Recap{7,3}));
        %             j=j+1;
        %         end
        %         if Recap{6,4}>0.3
        %             plot(wl(predSV2.UVEidx),ones(1,length(predSV2.UVEidx))+i,'*','Markersize',10), hold on
        %             i=i+0.1;
        %             legvl{j}=strcat('UVE : ', num2str(round(Recap{6,4},2)), '(', num2str(Recap{3,4}), '), ', num2str(Recap{7,4}));
        %             j=j+1;
        %         end
        %         if Recap{5,5}>0.5
        %             plot(wl(SV.idxCARSopt),ones(1,length(SV.idxCARSopt))+i,'*','Markersize',10), hold on
        %             i=i+0.1;
        %             legvl{j}=strcat('CARS : ', num2str(round(Recap{5,5},2)), '(', num2str(Recap{3,5}), '), ', num2str(Recap{7,5}));
        %             j=j+1;
        %         end
        if Recap{6,6}>0.3
            plot(wl(predSV2.RFidx),ones(1,length(predSV2.RFidx))+i,'*','Markersize',10), hold on
            i=i+0.1;
            legvl{j}=strcat('RF : ', num2str(round(Recap{6,6},2)), '(', num2str(Recap{3,6}), '), ', num2str(Recap{7,6}));
            j=j+1;
        end
        %         if Recap{5,7}>0.5
        %             plot(wl(SV.idxiRF(1:Recap{3,7})),ones(1,length(SV.idxiRF(1:Recap{3,7})))+i,'*','Markersize',10), hold on
        %             i=i+0.1;
        %             legvl{j}=strcat('iRF : ', num2str(round(Recap{5,7},2)), '(', num2str(Recap{3,7}), '), ', num2str(Recap{7,7}));
        %             j=j+1;
        %         end
        if Recap{6,8}>0.3
            plot(wl(predSV2.VCNidx),ones(1,length(predSV2.VCNidx))+i,'*','Markersize',10), hold on
            i=i+0.1;
            legvl{j}=strcat('VCN : ', num2str(round(Recap{6,8},2)), '(', num2str(Recap{3,8}), '), ', num2str(Recap{7,8}));
            j=j+1;
        end
        if strcmp(gaplsdo,'Yes')
            if Recap{6,9}>0.3
                plot(wl(predSV2.GAPLSidx),ones(1,length(predSV2.GAPLSidx))+i,'*','Markersize',10), hold on
                legvl{j}=strcat('GA-PLS : ', num2str(round(Recap{6,9},2)), '(', num2str(Recap{3,9}), '), ', num2str(Recap{7,9}));
            end
        end
        hold off
        grid on
        legend(legvl)
        xlabel('Wavelength (nm)')
        set(gca,'FontSize',16)
        % Création d'une structure modèle
        Modelc = questdlg('Do you want to create a model ?','Model creation','No','Yes','Yes');
        if strcmp(Modelc,'Yes')
            algo={'VIP','TP','UVE','CARS','RF','iRF','VCN','GA-PLS'};
            Modela = listdlg('PromptString','With which selection algorithm?','SelectionMode','single', 'ListString',algo);
            
            if strcmp(pseudoabs,'Yes')
                Model.Type='PseudoAbsorbance';
            else
                Model.Type='Reflectance';
            end
            Model.label=Yn;
            Model.SelectWave=algo{Modela};
            Model.cv=cvt;
            Model.wli=wl;
            
            if Modela==1
                %                 Model.Scal{1,1}=predSV2.predVIP.Scal;
                %                 Model.Ycal=predSV2.predVIP.Ycal;
                %                 Model.Sval{1,1}=predSV2.predVIP.Sval;
                %                 Model.Yval=predSV2.predVIP.Yval;
                %                 Model.Pret=plsPredOpt.Pretopt{1, 1}{1, 1}  ;
                %                 Model.Pret1=DataPLS.Pret;
                %                 Model.Pret1a=[DataPLS.Pret1; DataPLS.Pret2];
                %                 Model.wl=predSV2.VIPwl;
                %                 Model.B=predSV2.predVIP.brut.para.B;
                %                 Model.R2cal=predSV2.predVIP.brut.coeff.R2c;
                %                 Model.R2cv=predSV2.predVIP.brut.cv.r2cv;
                %                 Model.R2val=predSV2.predVIP.brut.coeff.R2p;
                %                 Model.SEC=predSV2.predVIP.brut.coeff.SEC;
                %                 Model.SECV=predSV2.predVIP.brut.cv.SECV;
                %                 Model.SEP=predSV2.predVIP.brut.coeff.SEP;
                %                 Model.RMSEP=predSV2.predVIP.brut.coeff.RMSEP;
                %                 Model.model=predSV2.predVIP.brut ;
                %                 Model.modelsansselectvar=PRED{opt,1};
            else if Modela==2
                    %                     Model.Scal{1,1}=predSV2.predTP.Scal;
                    %                     Model.Ycal=predSV2.predTP.Ycal;
                    %                     Model.Sval{1,1}=predSV2.predTP.Sval;
                    %                     Model.Yval=predSV2.predTP.Yval;
                    %                     Model.Pret=plsPredOpt.Pretopt{1, 1}{1, 1}  ;
                    %                     Model.Pret1=DataPLS.Pret;
                    %                     Model.Pret1a=[DataPLS.Pret1; DataPLS.Pret2];
                    %                     Model.wl=predSV2.TPwl;
                    %                     Model.B=predSV2.predTP.brut.para.B;
                    %                     Model.R2cal=predSV2.predTP.brut.coeff.R2c;
                    %                     Model.R2cv=predSV2.predTP.brut.cv.r2cv;
                    %                     Model.R2val=predSV2.predTP.brut.coeff.R2p;
                    %                     Model.SEC=predSV2.predTP.brut.coeff.SEC;
                    %                     Model.SECV=predSV2.predTP.brut.cv.SECV;
                    %                     Model.SEP=predSV2.predTP.brut.coeff.SEP;
                    %                     Model.RMSEP=predSV2.predTP.brut.coeff.RMSEP;
                    %                     Model.model=predSV2.predTP.brut ;
                    %                     Model.modelsansselectvar=PRED{opt,1};
                else if Modela==3
                        %                         Model.Scal{1,1}=predSV2.predUVE.Scal;
                        %                         Model.Ycal=predSV2.predUVE.Ycal;
                        %                         Model.Sval{1,1}=predSV2.predUVE.Sval;
                        %                         Model.Yval=predSV2.predUVE.Yval;
                        %                         Model.Pret=plsPredOpt.Pretopt{1, 1}{1, 1}  ;
                        %                         Model.Pret1=DataPLS.Pret;
                        %                         Model.Pret1a=[DataPLS.Pret1; DataPLS.Pret2];
                        %                         Model.wl=predSV2.UVEwl;
                        %                         Model.B=predSV2.predUVE.brut.para.B;
                        %                         Model.R2cal=predSV2.predUVE.brut.coeff.R2c;
                        %                         Model.R2cv=predSV2.predUVE.brut.cv.r2cv;
                        %                         Model.R2val=predSV2.predUVE.brut.coeff.R2p;
                        %                         Model.SEC=predSV2.predUVE.brut.coeff.SEC;
                        %                         Model.SECV=predSV2.predUVE.brut.cv.SECV;
                        %                         Model.SEP=predSV2.predUVE.brut.coeff.SEP;
                        %                         Model.RMSEP=predSV2.predUVE.brut.coeff.RMSEP;
                        %                         Model.model=predSV2.predUVE.brut ;
                        %                         Model.modelsansselectvar=PRED{opt,1};
                    else if Modela==4
                            %                             Model.Scal{1,1}=predSV2.predCARS.Scal;
                            %                             Model.Ycal=predSV2.predCARS.Ycal;
                            %                             Model.Sval{1,1}=predSV2.predCARS.Sval;
                            %                             Model.Yval=predSV2.predCARS.Yval;
                            %                             Model.Pret=plsPredOpt.Pretopt{1, 1}  ;
                            %                             Model.Pret1=DataPLS.Pret;
                            %                             Model.Pret1a=[DataPLS.Pret1; DataPLS.Pret2];
                            %                             Model.wl=predSV2.CARSwl;
                            %                             Model.B=predSV2.predCARS.brut.para.B;
                            %                             Model.R2cal=predSV2.predCARS.brut.coeff.R2c;
                            %                             Model.R2val=predSV2.predCARS.brut.coeff.R2p;
                            %                             Model.SEC=predSV2.predCARS.brut.coeff.SEC;
                            %                             Model.SEP=predSV2.predCARS.brut.coeff.SEP;
                            %                             Model.RMSEP=predSV2.predCARS.brut.coeff.RMSEP;
                            %                             Model.model=predSV2.predCARS.brut ;
                        else if Modela==5
                                Model.Scal{1,1}=predSV2.predRF.Scal;
                                Model.Ycal=predSV2.predRF.Ycal;
                                Model.Sval{1,1}=predSV2.predRF.Sval;
                                Model.Yval=predSV2.predRF.Yval;
                                Model.Pret=plsPredOpt.Pretopt{1, 1}{1, 1}  ;
                                Model.Pret1=DataPLS.Pret;
                                Model.Pret1a=[DataPLS.Pret1; DataPLS.Pret2];
                                Model.wl=predSV2.RFwl;
                                Model.B=predSV2.predRF.brut.para.B;
                                Model.R2cal=predSV2.predRF.brut.coeff.R2c;
                                Model.R2cv=predSV2.predRF.brut.cv.r2cv;
                                Model.R2val=predSV2.predRF.brut.coeff.R2p;
                                Model.SEC=predSV2.predRF.brut.coeff.SEC;
                                Model.SECV=predSV2.predRF.brut.cv.SECV;
                                Model.SEP=predSV2.predRF.brut.coeff.SEP;
                                Model.RMSEP=predSV2.predRF.brut.coeff.RMSEP;
                                Model.model=predSV2.predRF.brut ;
                                Model.modelsansselectvar=PRED{opt,1};
                            else if Modela==6
                                    %                                     Model.Scal{1,1}=predSV2.predIRF.Scal;
                                    %                                     Model.Ycal=predSV2.predIRF.Ycal;
                                    %                                     Model.Sval{1,1}=predSV2.predIRF.Sval;
                                    %                                     Model.Yval=predSV2.predIRF.Yval;
                                    %                                     Model.Pret=plsPredOpt.Pretopt{1, 1}  ;
                                    %                                     Model.Pret1=DataPLS.Pret;
                                    %                                     Model.Pret1a=[DataPLS.Pret1; DataPLS.Pret2];
                                    %                                     Model.wl=predSV2.IRFwl;
                                    %                                     Model.B=predSV2.predIRF.brut.para.B;
                                    %                                     Model.R2cal=predSV2.predIRF.brut.coeff.R2c;
                                    %                                     Model.R2val=predSV2.predIRF.brut.coeff.R2p;
                                    %                                     Model.SEC=predSV2.predIRF.brut.coeff.SEC;
                                    %                                     Model.SEP=predSV2.predIRF.brut.coeff.SEP;
                                    %                                     Model.RMSEP=predSV2.predIRF.brut.coeff.RMSEP;
                                    %                                     Model.model=predSV2.predIRF.brut ;
                                else if Modela==7
                                        Model.Scal{1,1}=predSV2.predVCN.Scal;
                                        Model.Ycal=predSV2.predVCN.Ycal;
                                        Model.Sval{1,1}=predSV2.predVCN.Sval;
                                        Model.Yval=predSV2.predVCN.Yval;
                                        Model.Pret=plsPredOpt.Pretopt{1, 1}{1, 1}  ;
                                        Model.Pret1=DataPLS.Pret;
                                        Model.Pret1a=[DataPLS.Pret1; DataPLS.Pret2];
                                        Model.wl=predSV2.VCNwl;
                                        Model.B=predSV2.predVCN.brut.para.B;
                                        Model.R2cal=predSV2.predVCN.brut.coeff.R2c;
                                        Model.R2cv=predSV2.predVCN.brut.cv.r2cv;
                                        Model.R2val=predSV2.predVCN.brut.coeff.R2p;
                                        Model.SEC=predSV2.predVCN.brut.coeff.SEC;
                                        Model.SECV=predSV2.predVCN.brut.cv.SECV;
                                        Model.SEP=predSV2.predVCN.brut.coeff.SEP;
                                        Model.RMSEP=predSV2.predVCN.brut.coeff.RMSEP;
                                        Model.model=predSV2.predVCN.brut ;
                                        Model.modelsansselectvar=PRED{opt,1};
                                    else if Modela==8
                                            Model.Scal{1,1}=predSV2.predGAPLS.Scal;
                                            Model.Ycal=predSV2.predGAPLS.Ycal;
                                            Model.Sval{1,1}=predSV2.predGAPLS.Sval;
                                            Model.Yval=predSV2.predGAPLS.Yval;
                                            Model.Pret=plsPredOpt.Pretopt{1, 1}{1, 1}  ;
                                            Model.Pret1=DataPLS.Pret;
                                            Model.Pret1a=[DataPLS.Pret1; DataPLS.Pret2];
                                            Model.wl=predSV2.GAPLSwl;
                                            Model.B=predSV2.predGAPLS.brut.para.B;
                                            Model.R2cal=predSV2.predGAPLS.brut.coeff.R2c;
                                            Model.R2cv=predSV2.predGAPLS.brut.cv.r2cv;
                                            Model.R2val=predSV2.predGAPLS.brut.coeff.R2p;
                                            Model.SEC=predSV2.predGAPLS.brut.coeff.SEC;
                                            Model.SECV=predSV2.predGAPLS.brut.cv.SECV;
                                            Model.SEP=predSV2.predGAPLS.brut.coeff.SEP;
                                            Model.RMSEP=predSV2.predGAPLS.brut.coeff.RMSEP;
                                            Model.model=predSV2.predGAPLS.brut ;
                                            Model.modelsansselectvar=PRED{opt,1};
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            Model.Recap{1,2}='All';
            Model.Recap{1,3}='Selected';
            
            Model.Recap{2,1}='Nb wl';
            Model.Recap{2,2}=length(wl);
            Model.Recap{2,3}=length(Model.wl);
            
            Model.Recap{3,1}='Nb vl';
            Model.Recap{3,3}=Model.model.para.nbcomp;
%             if strcmp(chxpret,'All')
%                 Model.Recap{3,2}=Model.modelsansselectvar.recap{Pretoptg+1,2};
%             else
                Model.Recap{3,2}=Model.modelsansselectvar.brut.para.nbcomp;
%             end
            
            Model.Recap{4,1}='SEC';
            Model.Recap{4,3}=Model.SEC;
%             if strcmp(chxpret,'All')
%                 Model.Recap{4,2}=Model.modelsansselectvar.recap{Pretoptg+1,3};
%             else
                Model.Recap{4,2}=Model.modelsansselectvar.brut.coeff.SEC;
%             end
            
            Model.Recap{5,1}='R2c';
            Model.Recap{5,3}=Model.R2cal;
%             if strcmp(chxpret,'All')
%                 Model.Recap{5,2}=Model.modelsansselectvar.recap{Pretoptg+1,4};
%             else
                Model.Recap{5,2}=Model.modelsansselectvar.brut.coeff.R2c;
%             end
            
            Model.Recap{6,1}='SECV';
            Model.Recap{6,3}=Model.SECV;
%             if strcmp(chxpret,'All')
%                 Model.Recap{6,2}=Model.modelsansselectvar.recap{Pretoptg+1,6};
%             else
                Model.Recap{6,2}=Model.modelsansselectvar.brut.cv.SECV;
%             end
            
            Model.Recap{7,1}='R2cv';
            Model.Recap{7,3}=Model.R2cv;
%             if strcmp(chxpret,'All')
%                 Model.Recap{7,2}=Model.modelsansselectvar.recap{Pretoptg+1,7};
%             else
                Model.Recap{7,2}=Model.modelsansselectvar.brut.cv.r2cv;
%             end
            
            Model.Recap{8,1}='SEP';
            Model.Recap{8,3}=Model.SEP;
%             if strcmp(chxpret,'All')
%                 Model.Recap{8,2}=Model.modelsansselectvar.recap{Pretoptg+1,9};
%             else
                Model.Recap{8,2}=Model.modelsansselectvar.brut.coeff.SEP;
%             end
            
            Model.Recap{9,1}='RMSEP';
            Model.Recap{9,3}=Model.model.coeff.RMSEP;
%             if strcmp(chxpret,'All')
%                 Model.Recap{9,2}=Model.modelsansselectvar.recap{Pretoptg+1,10};
%             else
                Model.Recap{9,2}=Model.modelsansselectvar.brut.coeff.RMSEP;
%             end
            
            Model.Recap{10,1}='R2p';
            Model.Recap{10,3}=Model.R2val;
%             if strcmp(chxpret,'All')
%                 Model.Recap{10,2}=Model.modelsansselectvar.recap{Pretoptg+1,11};
%             else
                Model.Recap{10,2}=Model.modelsansselectvar.brut.coeff.R2p;
%             end
            
            assignin('base','Modelvs',Model);
        end
    end
    
end

if strcmp(plst,'PLS2')||strcmp(plst,'PLS1 et PLS2')
    % Bootstrap
    R2cal=zeros(size(Ye,2),nbiterboostrap,14);
    R2cv=zeros(size(Ye,2),nbiterboostrap,14);
    R2val=zeros(size(Ye,2),nbiterboostrap,14);
    Nbvl=zeros(size(Ye,2),nbiterboostrap,14);
    
    % Check that there are several different values
    if size(Ye,1)/6>50
        lim=20;
    else
        lim=0.5*size(Ye,1)/6;
    end
    
    h = waitbar(0,'PLS-CV');
    iter=1;
    for k=1:nbiterboostrap
        % PLS-CV
        if sum(sum(abs(Ye)))>0&&length(unique(Ye))>lim
            pred=PLScv(squeeze(Sw(nbiterboostrap,:,:)),Ye,Num,wl,cv);
            PRED{k}=pred;
            
            for i=1:size(Y,2)
                R2cal(i,k,:)=cell2mat(PRED{k}.RECAP{i,1}  (2:end,4));
                R2cv(i,k,:)=cell2mat(PRED{k}.RECAP{i,1}  (2:end,7));
                R2val(i,k,:)=cell2mat(PRED{k}.RECAP{i,1}  (2:end,11));
                Nbvl(i,k,:)=cell2mat(PRED{k}.RECAP{i,1}  (2:end,2));
            end
        end
        h = waitbar(iter/(nbiterboostrap*size(Ye,2)));
        iter=iter+1;
    end
    close(h)
    
    % Optimal preprocesssing
    R2g=R2cal+R2cv+R2val;
    
    R2calopt=zeros(size(Ye,2),nbiterboostrap);
    R2cvopt=zeros(size(Ye,2),nbiterboostrap);
    R2valopt=zeros(size(Ye,2),nbiterboostrap);
    Nbvlopt=zeros(size(Ye,2),nbiterboostrap);
    Pret={'Raw' 'Detrend' 'MSC' 'SNV' 'SNVD' 'D1' 'D2' 'MSCD1' 'SNVD1' 'SNVDD1' 'MSCD2' 'SNVD2' 'SNVDD2' 'CR'};
    
    R2G=squeeze(mean(R2g));
    
    % Select optimal coefficients
    for k=1:nbiterboostrap
        [opt(k),~]=find(squeeze(R2G(k,:))==max(squeeze(R2G(k,:))));
        Pretoptc(k)=opt(k);
    end
    
    % Optimal preprocessing
    Pretoptg=mode(Pretoptc);
    
    R2calopt=squeeze(R2cal(:,:,Pretoptg));
    R2cvopt=squeeze(R2cv(:,:,Pretoptg));
    R2valopt=squeeze(R2val(:,:,Pretoptg));
    Nbvlopt=squeeze(Nbvl(:,:,Pretoptg));
    Pretopt=Pret(Pretoptg);
    
    % Saving
    plsPred.PRED=PRED;
    plsPred.R2cal=R2cal;
    plsPred.R2cv=R2cv;
    plsPred.R2val=R2val;
    plsPred.Nbvl=Nbvl;
    assignin('base', 'plsPredPLS2boostrap', plsPred);
    
    plsPredOpt.predopt=PRED;
    plsPredOpt.R2calopt=R2calopt;
    plsPredOpt.R2cvopt=R2cvopt;
    plsPredOpt.R2valopt=R2valopt;
    plsPredOpt.Nbvlopt=Nbvlopt;
    plsPredOpt.Pretopt=Pretopt;
    assignin('base', 'plsPredOptPLS2', plsPredOpt);
    
    %     for i=1:size(Ye,2)
    %         fig = figure;
    %         fig.PaperPositionMode = 'auto';
    %         fig.InvertHardcopy = 'off';
    %         set(fig,'units','characters', 'position',[0 0 a_characters(1) a_characters(2)], 'name', 'Résultats');
    %         subplot(221),hist(R2calopt(i,:),100),grid on,title('R² cal opt')
    %         subplot(222),hist(R2cvopt(i,:),100),grid on,title('R² cv opt')
    %         subplot(223),hist(R2valopt(i,:),100),grid on,title('R² val opt')
    %         subplot(224),hist(Nbvlopt(i,:),100),grid on,title('Nv VL opt')
    %     end
    
    %% Optimal set
    R2optg=R2calopt+R2cvopt+2*R2valopt;
    R2optG=mean(R2optg);
    [optsort,opt] = sort(R2optG,'descend');
    R2valoptG=mean(R2valopt,1);
    [optsortval,optval] = sort(R2valoptG,'descend');
    
    [~,idxset]=find(optsort>=optsort(1)-0.15);
    if length(idxset)>10
        idxset=idxset(1:5);
    else
        idxset=1:5;
    end
    
    % Optimal model saving
    if strcmp(pseudoabs,'Yes')
        Modelwsv.Type='PseudoAbsorbance';
    else
        Modelwsv.Type='Reflectance';
    end
    Modelwsv.label=Yn;
    Modelwsv.SelectWave='None';
    Modelwsv.cv=cvt;
    Modelwsv.Scal{1,1}=PRED{1,opt(1)}.Scal;
    Modelwsv.Ycal=PRED{1,opt(1)}.Ycal;
    Modelwsv.Sval{1,1}=PRED{1,opt(1)}.Sval;
    Modelwsv.Yval=PRED{1,opt(1)}.Yval;
    Pret={'Raw';'Detrend';'MSC';'SNV';'SNVD';'D1';'D2';'MSC+D1';'SNV+D1';'SNVD+D1';'MSC+D2';'SNV+D2';'SNVD+D2';'CR'};
    Modelwsv.Pret=Pret(Pretoptg);
    Modelwsv.Pret1=DataPLS.Pret;
    Modelwsv.Pret1a=[DataPLS.Pret1; DataPLS.Pret2];
    Modelwsv.wl=wl;
    if Pretoptg==1
        Modelwsv.B=PRED{1,opt(1)}.brut.para.B;
        Modelwsv.R2cal=PRED{1,opt(1)}.brut.coeff.R2c;
        Modelwsv.R2cv=PRED{1,opt(1)}.brut.cv.r2cv;
        Modelwsv.R2val=PRED{1,opt(1)}.brut.coeff.R2p;
        Modelwsv.SEC=PRED{1,opt(1)}.brut.coeff.SEC;
        Modelwsv.SECV=PRED{1,opt(1)}.brut.cv.SECV;
        Modelwsv.SEP=PRED{1,opt(1)}.brut.coeff.SEP;
        Modelwsv.RMSEP=PRED{1,opt(1)}.brut.coeff.RMSEP;
        Modelwsv.model=PRED{1,opt(1)}.brut ;
    else if Pretoptg==2
            Modelwsv.B=PRED{1,opt(1)}.detrend.para.B;
            Modelwsv.R2cal=PRED{1,opt(1)}.detrend.coeff.R2c;
            Modelwsv.R2cv=PRED{1,opt(1)}.detrend.cv.r2cv;
            Modelwsv.R2val=PRED{1,opt(1)}.detrend.coeff.R2p;
            Modelwsv.SEC=PRED{1,opt(1)}.detrend.coeff.SEC;
            Modelwsv.SECV=PRED{1,opt(1)}.detrend.cv.SECV;
            Modelwsv.SEP=PRED{1,opt(1)}.detrend.coeff.SEP;
            Modelwsv.RMSEP=PRED{1,opt(1)}.detrend.coeff.RMSEP;
            Modelwsv.model=PRED{1,opt(1)}.detrend ;
        else if Pretoptg==3
                Modelwsv.B=PRED{1,opt(1)}.msc.para.B;
                Modelwsv.R2cal=PRED{1,opt(1)}.msc.coeff.R2c;
                Modelwsv.R2cv=PRED{1,opt(1)}.msc.cv.r2cv;
                Modelwsv.R2val=PRED{1,opt(1)}.msc.coeff.R2p;
                Modelwsv.SEC=PRED{1,opt(1)}.msc.coeff.SEC;
                Modelwsv.SECV=PRED{1,opt(1)}.msc.cv.SECV;
                Modelwsv.SEP=PRED{1,opt(1)}.msc.coeff.SEP;
                Modelwsv.RMSEP=PRED{1,opt(1)}.msc.coeff.RMSEP;
                Modelwsv.model=PRED{1,opt(1)}.msc ;
            else if Pretoptg==4
                    Modelwsv.B=PRED{1,opt(1)}.snv.para.B;
                    Modelwsv.R2cal=PRED{1,opt(1)}.snv.coeff.R2c;
                    Modelwsv.R2cv=PRED{1,opt(1)}.snv.cv.r2cv;
                    Modelwsv.R2val=PRED{1,opt(1)}.snv.coeff.R2p;
                    Modelwsv.SEC=PRED{1,opt(1)}.snv.coeff.SEC;
                    Modelwsv.SECV=PRED{1,opt(1)}.snv.cv.SECV;
                    Modelwsv.SEP=PRED{1,opt(1)}.snv.coeff.SEP;
                    Modelwsv.RMSEP=PRED{1,opt(1)}.snv.coeff.RMSEP;
                    Modelwsv.model=PRED{1,opt(1)}.snv ;
                else if Pretoptg==5
                        Modelwsv.B=PRED{1,opt(1)}.snvd.para.B;
                        Modelwsv.R2cal=PRED{1,opt(1)}.snvd.coeff.R2c;
                        Modelwsv.R2cv=PRED{1,opt(1)}.snvd.cv.r2cv;
                        Modelwsv.R2val=PRED{1,opt(1)}.snvd.coeff.R2p;
                        Modelwsv.SEC=PRED{1,opt(1)}.snvd.coeff.SEC;
                        Modelwsv.SECV=PRED{1,opt(1)}.snvd.cv.SECV;
                        Modelwsv.SEP=PRED{1,opt(1)}.snvd.coeff.SEP;
                        Modelwsv.RMSEP=PRED{1,opt(1)}.snvd.coeff.RMSEP;
                        Modelwsv.model=PRED{1,opt(1)}.snvd ;
                    else if Pretoptg==6
                            Modelwsv.B=PRED{1,opt(1)}.d1.para.B;
                            Modelwsv.R2cal=PRED{1,opt(1)}.d1.coeff.R2c;
                            Modelwsv.R2cv=PRED{1,opt(1)}.d1.cv.r2cv;
                            Modelwsv.R2val=PRED{1,opt(1)}.d1.coeff.R2p;
                            Modelwsv.SEC=PRED{1,opt(1)}.d1.coeff.SEC;
                            Modelwsv.SECV=PRED{1,opt(1)}.d1.cv.SECV;
                            Modelwsv.SEP=PRED{1,opt(1)}.d1.coeff.SEP;
                            Modelwsv.RMSEP=PRED{1,opt(1)}.d1.coeff.RMSEP;
                            Modelwsv.model=PRED{1,opt(1)}.d1 ;
                        else if Pretoptg==7
                                Modelwsv.B=PRED{1,opt(1)}.d2.para.B;
                                Modelwsv.R2cal=PRED{1,opt(1)}.d2.coeff.R2c;
                                Modelwsv.R2cv=PRED{1,opt(1)}.d2.cv.r2cv;
                                Modelwsv.R2val=PRED{1,opt(1)}.d2.coeff.R2p;
                                Modelwsv.SEC=PRED{1,opt(1)}.d2.coeff.SEC;
                                Modelwsv.SECV=PRED{1,opt(1)}.d2.cv.SECV;
                                Modelwsv.SEP=PRED{1,opt(1)}.d2.coeff.SEP;
                                Modelwsv.RMSEP=PRED{1,opt(1)}.d2.coeff.RMSEP;
                                Modelwsv.model=PRED{1,opt(1)}.d2 ;
                            else if Pretoptg==8
                                    Modelwsv.B=PRED{1,opt(1)}.mscd1.para.B;
                                    Modelwsv.R2cal=PRED{1,opt(1)}.mscd1.coeff.R2c;
                                    Modelwsv.R2cv=PRED{1,opt(1)}.mscd1.cv.r2cv;
                                    Modelwsv.R2val=PRED{1,opt(1)}.mscd1.coeff.R2p;
                                    Modelwsv.SEC=PRED{1,opt(1)}.mscd1.coeff.SEC;
                                    Modelwsv.SECV=PRED{1,opt(1)}.mscd1.cv.SECV;
                                    Modelwsv.SEP=PRED{1,opt(1)}.mscd1.coeff.SEP;
                                    Modelwsv.RMSEP=PRED{1,opt(1)}.mscd1.coeff.RMSEP;
                                    Modelwsv.model=PRED{1,opt(1)}.mscd1 ;
                                else if Pretoptg==9
                                        Modelwsv.B=PRED{1,opt(1)}.mscd2.para.B;
                                        Modelwsv.R2cal=PRED{1,opt(1)}.mscd2.coeff.R2c;
                                        Modelwsv.R2cv=PRED{1,opt(1)}.mscd2.cv.r2cv;
                                        Modelwsv.R2val=PRED{1,opt(1)}.mscd2.coeff.R2p;
                                        Modelwsv.SEC=PRED{1,opt(1)}.mscd2.coeff.SEC;
                                        Modelwsv.SECV=PRED{1,opt(1)}.mscd2.cv.SECV;
                                        Modelwsv.SEP=PRED{1,opt(1)}.mscd2.coeff.SEP;
                                        Modelwsv.RMSEP=PRED{1,opt(1)}.mscd2.coeff.RMSEP;
                                        Modelwsv.model=PRED{1,opt(1)}.mscd2 ;
                                    else if Pretoptg==10
                                            Modelwsv.B=PRED{1,opt(1)}.snvd1.para.B;
                                            Modelwsv.R2cal=PRED{1,opt(1)}.snvd1.coeff.R2c;
                                            Modelwsv.R2cv=PRED{1,opt(1)}.snvd1.cv.r2cv;
                                            Modelwsv.R2val=PRED{1,opt(1)}.snvd1.coeff.R2p;
                                            Modelwsv.SEC=PRED{1,opt(1)}.snvd1.coeff.SEC;
                                            Modelwsv.SECV=PRED{1,opt(1)}.snvd1.cv.SECV;
                                            Modelwsv.SEP=PRED{1,opt(1)}.snvd1.coeff.SEP;
                                            Modelwsv.RMSEP=PRED{1,opt(1)}.snvd1.coeff.RMSEP;
                                            Modelwsv.model=PRED{1,opt(1)}.snvd1 ;
                                        else if Pretoptg==11
                                                Modelwsv.B=PRED{1,opt(1)}.snvd2.para.B;
                                                Modelwsv.R2cal=PRED{1,opt(1)}.snvd2.coeff.R2c;
                                                Modelwsv.R2cv=PRED{1,opt(1)}.snvd2.cv.r2cv;
                                                Modelwsv.R2val=PRED{1,opt(1)}.snvd2.coeff.R2p;
                                                Modelwsv.SEC=PRED{1,opt(1)}.snvd2.coeff.SEC;
                                                Modelwsv.SECV=PRED{1,opt(1)}.snvd2.cv.SECV;
                                                Modelwsv.SEP=PRED{1,opt(1)}.snvd2.coeff.SEP;
                                                Modelwsv.RMSEP=PRED{1,opt(1)}.snvd2.coeff.RMSEP;
                                                Modelwsv.model=PRED{1,opt(1)}.snvd2 ;
                                            else if Pretoptg==12
                                                    Modelwsv.B=PRED{1,opt(1)}.snvdd1.para.B;
                                                    Modelwsv.R2cal=PRED{1,opt(1)}.snvdd1.coeff.R2c;
                                                    Modelwsv.R2cv=PRED{1,opt(1)}.snvdd1.cv.r2cv;
                                                    Modelwsv.R2val=PRED{1,opt(1)}.snvdd1.coeff.R2p;
                                                    Modelwsv.SEC=PRED{1,opt(1)}.snvdd1.coeff.SEC;
                                                    Modelwsv.SECV=PRED{1,opt(1)}.snvdd1.cv.SECV;
                                                    Modelwsv.SEP=PRED{1,opt(1)}.snvdd1.coeff.SEP;
                                                    Modelwsv.RMSEP=PRED{1,opt(1)}.snvdd1.coeff.RMSEP;
                                                    Modelwsv.model=PRED{1,opt(1)}.snvdd1 ;
                                                else if Pretoptg==13
                                                        Modelwsv.B=PRED{1,opt(1)}.snvdd2.para.B;
                                                        Modelwsv.R2cal=PRED{1,opt(1)}.snvdd2.coeff.R2c;
                                                        Modelwsv.R2cv=PRED{1,opt(1)}.snvdd2.cv.r2cv;
                                                        Modelwsv.R2val=PRED{1,opt(1)}.snvdd2.coeff.R2p;
                                                        Modelwsv.SEC=PRED{1,opt(1)}.snvdd2.coeff.SEC;
                                                        Modelwsv.SECV=PRED{1,opt(1)}.snvdd2.cv.SECV;
                                                        Modelwsv.SEP=PRED{1,opt(1)}.snvdd2.coeff.SEP;
                                                        Modelwsv.RMSEP=PRED{1,opt(1)}.snvdd2.coeff.RMSEP;
                                                        Modelwsv.model=PRED{1,opt(1)}.snvdd2 ;
                                                    else if Pretoptg==14
                                                            Modelwsv.B=PRED{1,opt(1)}.cr.para.B;
                                                            Modelwsv.R2cal=PRED{1,opt(1)}.cr.coeff.R2c;
                                                            Modelwsv.R2cv=PRED{1,opt(1)}.cr.cv.r2cv;
                                                            Modelwsv.R2val=PRED{1,opt(1)}.cr.coeff.R2p;
                                                            Modelwsv.SEC=PRED{1,opt(1)}.cr.coeff.SEC;
                                                            Modelwsv.SECV=PRED{1,opt(1)}.cr.cv.SECV;
                                                            Modelwsv.SEP=PRED{1,opt(1)}.cr.coeff.SEP;
                                                            Modelwsv.RMSEP=PRED{1,opt(1)}.cr.coeff.RMSEP;
                                                            Modelwsv.model=PRED{1,opt(1)}.cr ;
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    assignin('base','Modelwvs',Modelwsv);
    
    for e=1:length(idxset)
        %% Wavelengths selection
        if Pretoptg==1
            Swsv=squeeze(Sw(opt(idxset(e)),:,:));
        else
            Swsv = AllPret(squeeze(Sw(opt(idxset(e)),:,:)),wl,Pretoptg-1);
        end
        
        %         for z=1:size(Ye,2)
        SV{e}=AllSelectVar(Swsv,Ye,wl,mode(Nbvlopt));
        %
        %             % GA-PLS
        %             [b,~,sel]=gaplssp([Swsv Ye(:,z)],500);
        %             SV{e}.wGAPLS=sel;
        %             SV{e}.idxGAPLS=b;
        %
        %             SV{e}.setopt=opt(idxset(e));
        %         end
        
        % VIP
        %     predSV.predVIP=PLScv(squeeze(Sw(opt,:,:)),Ye(:,i),Num,wl(SV.idxVIPopt),cv,Pretoptg,SV.idxVIPopt);
        %     predSV.VIPwl=wl(SV.idxVIPopt);
        
        %         % TP
        %         for k=10:2*length(wl)/5
        %             predTP{k}=PLScv(Swsv,Ye,Num,wl(SV{e}.idxTP(1:k)),cv,1,SV{e}.idxTP(1:k));
        %             TPcal(k,:)=predTP{1, k}.brut.coeff.R2c;
        %             TPval(k,:)=predTP{1, k}.brut.coeff.R2p;
        %             TPcv(k,:)=predTP{1, k}.brut.cv.r2cv;
        %         end
        %         TPg=TPcal+TPcv+TPval;
        %         TPgg=mean(TPg,2);
        %         [idx,~]=find(TPgg(10:end)==max(TPgg(10:end)));
        %         predSV{e}.predTP=predTP{idx+9};
        %         predSV{e}.TPwl=wl(SV{e}.idxTP(1:idx));
        %
        %         % UVE
        %         for k=10:2*length(wl)/5
        %             predUVE{k}=PLScv(Swsv,Ye,Num,wl(SV{e}.idxUVE(1:k)),cv,1,SV{e}.idxUVE(1:k));
        %             UVEcal(k,:)=predUVE{1, k}.brut.coeff.R2c;
        %             UVEval(k,:)=predUVE{1, k}.brut.coeff.R2p;
        %             UVEcv(k,:)=predUVE{1, k}.brut.cv.r2cv;
        %         end
        %         UVEg=UVEcal+UVEcv+UVEval;
        %         UVEgg=mean(UVEg,2);
        %         [idx,~]=find(UVEgg(10:end)==max(UVEgg(10:end)));
        %         predSV{e}.predUVE=predUVE{idx+9};
        %         predSV{e}.UVEwl=wl(SV{e}.idxUVE(1:idx));
        
        %     % CARS
        %     predSV.predCARS=PLScv(squeeze(Sw(opt,:,:)),Ye(:,i),Num,wl(SV.idxCARSopt),cv,Pretoptg,SV.idxCARSopt);
        %     predSV.CARSwl=wl(SV.idxCARSopt);
        
        % RF
        for k=10:4*length(wl)/5
            predRF{k}=PLScv(Swsv,Ye,Num,wl(SV{e}.idxRF(1:k)),cv,1,SV{e}.idxRF(1:k));
            RFcal(k,:)=predRF{1, k}.brut.coeff.R2c;
            RFval(k,:)=predRF{1, k}.brut.coeff.R2p;
            RFcv(k,:)=predRF{1, k}.brut.cv.r2cv;
        end
        RFg=RFcal+RFcv+RFval;
        RFgg=mean(RFg,2);
        [idx,~]=find(RFgg(10:end)==max(RFgg(10:end)));
        predSV{e}.predRF=predRF{idx+9};
        predSV{e}.RFwl=wl(SV{e}.idxRF(1:idx));
        
        %     % IRF
        %     for k=10:2*length(wl)/5
        %         predIRF{k}=PLScv(squeeze(Sw(opt,:,:)),Ye(:,i),Num,wl(SV.idxiRF(1:k)),cv,Pretoptg,SV.idxiRF(1:k));
        %         IRFcal(k)=predIRF{1, k}.brut.coeff.R2c;
        %         IRFval(k)=predIRF{1, k}.brut.coeff.R2p;
        %     end
        %     IRFg=IRFcal+IRFval;
        %     [~,idx]=find(IRFg==max(IRFg));
        %     predSV.predIRF=predIRF{idx};
        %     predSV.IRFwl=wl(SV.idxiRF(1:idx));
        
        % IRIV % BUG SUR LA PLS
        %     predSV.predIRIS=PLScv(squeeze(Swsv(:,SV.idxIRIVopt)),Ye(:,i),Num,wl(SV.idxIRIVopt),cv,1);
        
        % VCN
        predSV{e}.predVCN=PLScv(Swsv,Ye,Num,wl(SV{e}.idxVCNopt),cv,1,SV{e}.idxVCNopt);
        predSV{e}.VCNwl=wl(SV{e}.idxVCNopt);
        
        %         % GA-PLS
        %         for k=5:2*length(wl)/5
        %             predGAPLS{k}=PLScv(Swsv,Ye(:,i),Num,wl(SV{e}.idxGAPLS(1:k)),cv,1,SV{e}.idxGAPLS(1:k));
        %             GAPLScal(k)=predGAPLS{1, k}.brut.coeff.R2c;
        %             GAPLSval(k)=predGAPLS{1, k}.brut.coeff.R2p;
        %             GAPLScv(k)=predGAPLS{1, k}.brut.cv.r2cv;
        %         end
        %         GAPLSg=GAPLScal+GAPLScv+GAPLSval;
        %         [~,idx]=find(GAPLSg(5:end)==max(GAPLSg(5:end)));
        %         predSV{e}.predGAPLS=predGAPLS{idx+4};
        %         predSV{e}.GAPLSwl=wl(SV{e}.idxGAPLS(1:idx));
        
    end
    %set alea
    settest=opt(1);
    TPwlall=[];UVEwlall=[];RFwlall=[];VCNwlall=[];GAPLSwlall=[];
    for ibstr=1:e
        %         TPwlall=[TPwlall predSV{1, ibstr}.TPwl];
        %         UVEwlall=[UVEwlall predSV{1, ibstr}.UVEwl];
        RFwlall=[RFwlall predSV{1, ibstr}.RFwl];
        VCNwlall=[VCNwlall predSV{1, ibstr}.VCNwl];
        %         GAPLSwlall=[GAPLSwlall predSV{1, ibstr}.GAPLSwl];
    end
    %     % TP
    %     [CTP,ia,ic] = unique(sort(TPwlall));
    %     ia=[ia; length(TPwlall)];
    %     ianbTP=ia(2:end)-ia(1:end-1);
    %     ianbTPmax=max(ianbTP);
    %     ianbTPmin=min(ianbTP);
    %     idxwlTP=zeros(1,length(CTP));
    %     % idx des longueur d'onde
    %     for jwl=1:length(CTP)
    %         idxwlTP(jwl)=find(abs(CTP(jwl)-wl)==min(abs(CTP(jwl)-wl)));
    %     end
    %     % tris par importance et aléatoirement
    %     idxwlTPsort=[];
    %     for iwl=ianbTPmax:-1:ianbTPmin
    %         idwwlTPiwl=idxwlTP(ianbTP==iwl);
    %         aleawl=rand(1,length(idwwlTPiwl));
    %         [~,idxalea]=sort(aleawl);
    %         idxwlTPsort=[idxwlTPsort idwwlTPiwl(idxalea)];
    %     end
    %     % test PLS
    %     if Pretoptg==1
    %         Swsv=squeeze(Sw(settest,:,:));
    %     else
    %         Swsv = AllPret(squeeze(Sw(settest,:,:)),wl,Pretoptg-1);
    %     end
    %     if length(idxwlTPsort)<2*length(wl)/5
    %         lim=length(idxwlTPsort);
    %     else
    %         lim=2*length(wl)/5;
    %     end
    %     for k=5:lim
    %         predTP2{k}=PLScv(Swsv,Ye,Num,wl(idxwlTPsort(1:k)),cv,1,idxwlTPsort(1:k));
    %         TPcal2(k,:)=predTP2{1, k}.brut.coeff.R2c;
    %         TPval2(k,:)=predTP2{1, k}.brut.coeff.R2p;
    %         TPcv2(k,:)=predTP2{1, k}.brut.cv.r2cv;
    %     end
    %     TPg2=TPcal2+TPcv2+TPval2;
    %     TPg2g=mean(TPg2,2)';
    %     TPval2m=mean(TPval2,2)';
    %     if max(TPval2m)>0.75
    %         [TPgs,TPgidx]=sort(TPval2m,'descend');
    %         idx=min(TPgidx(TPgs>0.75));
    %     else if max(TPval2m)>0.5
    %             [TPgs,TPgidx]=sort(TPval2m,'descend');
    %             idx=min(TPgidx(TPgs>0.5));
    %         else
    %             [~,idx]=find(TPval2m(5:end)==max(TPval2m(5:end)));
    %             idx=idx+4;
    %         end
    %     end
    %     predSV2.predTP=predTP2{idx(1)};
    %     predSV2.TPwl=wl(idxwlTPsort(1:idx(1)));
    %     predSV2.TPidx=idxwlTPsort(1:idx(1));
    %
    %     % UVE
    %     [CUVE,ia,ic] = unique(sort(UVEwlall));
    %     ia=[ia; length(UVEwlall)];
    %     ianbUVE=ia(2:end)-ia(1:end-1);
    %     ianbUVEmax=max(ianbUVE);
    %     ianbUVEmin=min(ianbUVE);
    %     idxwlUVE=zeros(1,length(CUVE));
    %     % idx des longueur d'onde
    %     for jwl=1:length(CUVE)
    %         idxwlUVE(jwl)=find(abs(CUVE(jwl)-wl)==min(abs(CUVE(jwl)-wl)));
    %     end
    %     % tris par importance et aléatoirement
    %     idxwlUVEsort=[];
    %     for iwl=ianbUVEmax:-1:ianbUVEmin
    %         idwwlUVEiwl=idxwlUVE(ianbUVE==iwl);
    %         aleawl=rand(1,length(idwwlUVEiwl));
    %         [~,idxalea]=sort(aleawl);
    %         idxwlUVEsort=[idxwlUVEsort idwwlUVEiwl(idxalea)];
    %     end
    %     % test PLS
    %     if Pretoptg==1
    %         Swsv=squeeze(Sw(settest,:,:));
    %     else
    %         Swsv = AllPret(squeeze(Sw(settest,:,:)),wl,Pretoptg-1);
    %     end
    %     if length(idxwlUVEsort)<2*length(wl)/5
    %         lim=length(idxwlUVEsort);
    %     else
    %         lim=2*length(wl)/5;
    %     end
    %     for k=5:lim
    %         predUVE2{k}=PLScv(Swsv,Ye,Num,wl(idxwlUVEsort(1:k)),cv,1,idxwlUVEsort(1:k));
    %         UVEcal2(k,:)=predUVE2{1, k}.brut.coeff.R2c;
    %         UVEval2(k,:)=predUVE2{1, k}.brut.coeff.R2p;
    %         UVEcv2(k,:)=predUVE2{1, k}.brut.cv.r2cv;
    %     end
    %     UVEg2=UVEcal2+UVEcv2+UVEval2;
    %     UVEg2g=mean(UVEg2,2)';
    %     UVEval2g=mean(UVEval2,2)';
    %     if max(UVEval2g)>0.75
    %         [UVEgs,UVEgidx]=sort(UVEval2g,'descend');
    %         idx=min(UVEgidx(UVEgs>0.75));
    %     else if max(UVEval2g)>0.5
    %             [UVEgs,UVEgidx]=sort(UVEval2g,'descend');
    %             idx=min(UVEgidx(UVEgs>0.5));
    %         else
    %             [~,idx]=find(UVEval2g(5:end)==max(UVEval2g(5:end)));
    %             idx=idx+4;
    %         end
    %     end
    %     predSV2.predUVE=predUVE2{idx(1)};
    %     predSV2.UVEwl=wl(idxwlUVEsort(1:idx(1)));
    %     predSV2.UVEidx=idxwlUVEsort(1:idx(1));
    
    % RF
    [CRF,ia,ic] = unique(sort(RFwlall));
    ia=[ia; length(RFwlall)];
    ianbRF=ia(2:end)-ia(1:end-1);
    ianbRFmax=max(ianbRF);
    ianbRFmin=min(ianbRF);
    idxwlRF=zeros(1,length(CRF));
    % Index of the wavelengths
    for jwl=1:length(CRF)
        idxwlRF(jwl)=find(abs(CRF(jwl)-wl)==min(abs(CRF(jwl)-wl)));
    end
    % Sort the wavelengths
    idxwlRFsort=[];
    for iwl=ianbRFmax:-1:ianbRFmin
        idwwlRFiwl=idxwlRF(ianbRF==iwl);
        aleawl=rand(1,length(idwwlRFiwl));
        [~,idxalea]=sort(aleawl);
        idxwlRFsort=[idxwlRFsort idwwlRFiwl(idxalea)];
    end
    % test PLS
    if Pretoptg==1
        Swsv=squeeze(Sw(settest,:,:));
    else
        Swsv = AllPret(squeeze(Sw(settest,:,:)),wl,Pretoptg-1);
    end
    
    lim=min([4*length(wl)/5 length(idxwlRFsort)]);
    % zz=0;
    % while RFval2(idx)<0.7
    %     zz+1
    for k=5:lim
        predRF2{k}=PLScv(Swsv,Ye,Num,wl(idxwlRFsort(1:k)),cv,1,idxwlRFsort(1:k));
        RFcal2(k,:)=predRF2{1, k}.brut.coeff.R2c;
        RFval2(k,:)=predRF2{1, k}.brut.coeff.R2p;
        RFval2rmsep(k,:)=predRF2{1, k}.brut.coeff.RMSEP;
        RFcv2(k,:)=predRF2{1, k}.brut.cv.r2cv;
    end
    RFg2=RFcal2+RFcv2+RFval2;
    RFg2g=mean(RFg2,2)';
    RFval2g=mean(RFval2,2)';
    
    fig = figure;
    fig.PaperPositionMode = 'auto';
    fig.InvertHardcopy = 'off';
    set(fig,'units','characters', 'position',[0 0 a_characters(1) a_characters(2)], 'name', 'RF selected wavelength');
    plot(CRF,ianbRF,'*')
    grid on
    title('RF selected of wavelengths')
    
    if strcmp(chxwl,'Manual')
        fig = figure;
        fig.PaperPositionMode = 'auto';
        fig.InvertHardcopy = 'off';
        set(fig,'units','characters', 'position',[0 0 a_characters(1) a_characters(2)], 'name','RF selected wavelength');
        subplot(211)
        plot(5:lim,RFval2(5:lim,:))
        hold on
        plot(5:lim,RFval2g(5:lim),'k','linewidth',5)
        grid on
        title('RF selected wavelength R²val')
        lg=[Yn 'Mean'];
        legend(lg)
        ylim([-0.1 1])
        subplot(212)
        plot(5:lim,RFval2rmsep(5:lim,:))
        hold on
        plot(5:lim,mean(RFval2rmsep(5:lim,:),2),'k','linewidth',5)
        grid on
        title('RF selected wavelength RMSEP')
        ylim([min(abs([mode(RFval2rmsep(5:lim,:))-2*std(RFval2rmsep(5:lim,:)) 0])) min([max(mode(RFval2rmsep(5:lim,:))+2*std(RFval2rmsep(5:lim,:))) max(max(Ye))/2])])
        
        if max(RFval2g)>0.9
            [RFgs,RFgidx]=sort(RFval2g,'descend');
            idxit=min(RFgidx(RFgs>0.9));
        else if max(RFval2g)>0.75
                [RFgs,RFgidx]=sort(RFval2g,'descend');
                idxit=min(RFgidx(RFgs>0.75));
            else if max(RFval2g)>0.5
                    [RFgs,RFgidx]=sort(RFval2g,'descend');
                    idxit=min(RFgidx(RFgs>0.5));
                else
                    [~,idxit]=find(RFval2g(5:end)==max(RFval2g(5:end)));
                    idxit=idxit+4;
                end
            end
        end
        
        idxc = inputdlg('How many wavelengths do you want to work with?','How many wavelengths do you want to work with?',[1 35],cellstr(num2str(idxit)));
        idx = str2num(idxc{1});
    else if strcmp(chxwl,'Automatique')
            if max(RFval2g)>0.9
                [RFgs,RFgidx]=sort(RFval2g,'descend');
                idx=min(RFgidx(RFgs>0.9));
            else if max(RFval2g)>0.75
                    [RFgs,RFgidx]=sort(RFval2g,'descend');
                    idx=min(RFgidx(RFgs>0.75));
                else if max(RFval2g)>0.5
                        [RFgs,RFgidx]=sort(RFval2g,'descend');
                        idx=min(RFgidx(RFgs>0.5));
                    else
                        [~,idx]=find(RFval2g(5:end)==max(RFval2g(5:end)));
                        idx=idx+4;
                    end
                end
            end
        end
    end
    % end
    predSV2.predRF=predRF2{idx(1)};
    predSV2.RFwl=wl(idxwlRFsort(1:idx(1)));
    predSV2.RFidx=idxwlRFsort(1:idx(1));
    
    % VCN
    [CVCN,ia,ic] = unique(sort(VCNwlall));
    ia=[ia; length(VCNwlall)];
    ianbVCN=ia(2:end)-ia(1:end-1);
    ianbVCNmax=max(ianbVCN);
    ianbVCNmin=min(ianbVCN);
    idxwlVCN=zeros(1,length(CVCN));
    % Index of the wavelengths
    for jwl=1:length(CVCN)
        idxwlVCN(jwl)=find(abs(CVCN(jwl)-wl)==min(abs(CVCN(jwl)-wl)));
    end
    % Sort of the wavelengths
    idxwlVCNsort=[];
    for iwl=ianbVCNmax:-1:ianbVCNmin
        idwwlVCNiwl=idxwlVCN(ianbVCN==iwl);
        aleawl=rand(1,length(idwwlVCNiwl));
        [~,idxalea]=sort(aleawl);
        idxwlVCNsort=[idxwlVCNsort idwwlVCNiwl(idxalea)];
    end
    % test PLS
    if Pretoptg==1
        Swsv=squeeze(Sw(settest,:,:));
    else
        Swsv = AllPret(squeeze(Sw(settest,:,:)),wl,Pretoptg-1);
    end
    if length(idxwlVCNsort)<2*length(wl)/5
        lim=length(idxwlVCNsort);
    else
        lim=2*length(wl)/5;
    end
    for k=5:lim
        predVCN2{k}=PLScv(Swsv,Ye,Num,wl(idxwlVCNsort(1:k)),cv,1,idxwlVCNsort(1:k));
        VCNcal2(k,:)=predVCN2{1, k}.brut.coeff.R2c;
        VCNval2(k,:)=predVCN2{1, k}.brut.coeff.R2p;
        VCNcv2(k,:)=predVCN2{1, k}.brut.cv.r2cv;
    end
    VCNg2=VCNcal2+VCNcv2+VCNval2;
    VCNg2g=mean(VCNg2,2)';
    VCNval2g=mean(VCNval2,2)';
    if max(VCNval2g)>0.75
        [VCNgs,VCNgidx]=sort(VCNval2g,'descend');
        idx=min(VCNgidx(VCNgs>0.75));
    else if max(VCNval2g)>0.5
            [VCNgs,VCNgidx]=sort(VCNval2g,'descend');
            idx=min(VCNgidx(VCNgs>0.5));
        else
            [~,idx]=find(VCNval2g(5:end)==max(VCNval2g(5:end)));
            idx=idx+4;
        end
    end
    predSV2.predVCN=predVCN2{idx(1)};
    predSV2.VCNwl=wl(idxwlVCNsort(1:idx(1)));
    predSV2.VCNidx=idxwlVCNsort(1:idx(1));
    
    %     % GAPLS
    %     [CGAPLS,ia,ic] = unique(sort(GAPLSwlall));
    %     ia=[ia; length(GAPLSwlall)];
    %     ianbGAPLS=ia(2:end)-ia(1:end-1);
    %     ianbGAPLSmax=max(ianbGAPLS);
    %     ianbGAPLSmin=min(ianbGAPLS);
    %     idxwlGAPLS=zeros(1,length(CGAPLS));
    %     % idx des longueur d'onde
    %     for jwl=1:length(CGAPLS)
    %         idxwlGAPLS(jwl)=find(abs(CGAPLS(jwl)-wl)==min(abs(CGAPLS(jwl)-wl)));
    %     end
    %     % tris par importance et aléatoirement
    %     idxwlGAPLSsort=[];
    %     for iwl=ianbGAPLSmax:-1:ianbGAPLSmin
    %         idwwlGAPLSiwl=idxwlGAPLS(ianbGAPLS==iwl);
    %         aleawl=rand(1,length(idwwlGAPLSiwl));
    %         [~,idxalea]=sort(aleawl);
    %         idxwlGAPLSsort=[idxwlGAPLSsort idwwlGAPLSiwl(idxalea)];
    %     end
    %     % test PLS
    %     if Pretoptg==1
    %         Swsv=squeeze(Sw(settest,:,:));
    %     else
    %         Swsv = AllPret(squeeze(Sw(settest,:,:)),wl,Pretoptg-1);
    %     end
    %     if length(idxwlGAPLSsort)<2*length(wl)/5
    %         lim=length(idxwlGAPLSsort);
    %     else
    %         lim=2*length(wl)/5;
    %     end
    %     for k=5:lim
    %         predGAPLS2{k}=PLScv(Swsv,Ye(:,i),Num,wl(idxwlGAPLSsort(1:k)),cv,1,idxwlGAPLSsort(1:k));
    %         GAPLScal2(k)=predGAPLS2{1, k}.brut.coeff.R2c;
    %         GAPLSval2(k)=predGAPLS2{1, k}.brut.coeff.R2p;
    %         GAPLScv2(k)=predGAPLS2{1, k}.brut.cv.r2cv;
    %     end
    %     GAPLSg2=GAPLScal2+GAPLScv2+GAPLSval2;
    %     if max(GAPLSval2)>0.9
    %         [GAPLSgs,GAPLSgidx]=sort(GAPLSval2,'descend');
    %         idx=min(GAPLSgidx(GAPLSgs>0.9));
    %     else if max(GAPLSval2)>0.75
    %             [GAPLSgs,GAPLSgidx]=sort(GAPLSval2,'descend');
    %             idx=min(GAPLSgidx(GAPLSgs>0.75));
    %         else if max(GAPLSval2)>0.5
    %                 [GAPLSgs,GAPLSgidx]=sort(GAPLSval2,'descend');
    %                 idx=min(GAPLSgidx(GAPLSgs>0.5));
    %             else
    %                 [~,idx]=find(GAPLSval2(5:end)==max(GAPLSval2(5:end)));
    %                 idx=idx+4;
    %             end
    %         end
    %     end
    %     predSV2.predGAPLS=predGAPLS2{idx(1)};
    %     predSV2.GAPLSwl=wl(idxwlGAPLSsort(1:idx(1)));
    %     predSV2.GAPLSidx=idxwlGAPLSsort(1:idx(1));
    
    assignin('base','resSV',SV);
    assignin('base','predSV',predSV);
    assignin('base','predSV2',predSV2);
    % Tableau récapitulatif :
    Recap{1,1}='Method';
    %     Recap{1,2}='VIP';
    %     Recap{1,3}='TP';
    %     Recap{1,4}='UVE';
    %     Recap{1,5}='CARS';
    Recap{1,6}='RF';
    %     Recap{1,7}='iRF';
    Recap{1,8}='VCN';
    %     Recap{1,9}='GA-PLS';
    
    Recap{2,1}='VL';
    %     Recap{2,2}=predSV.predVIP.brut.para.nbcomp;
    %     Recap{2,3}=predSV2.predTP.brut.para.nbcomp;
    %     Recap{2,4}=predSV2.predUVE.brut.para.nbcomp;
    %     Recap{2,5}=predSV.predCARS.brut.para.nbcomp;
    Recap{2,6}=predSV2.predRF.brut.para.nbcomp;
    %     Recap{2,7}=predSV.predIRF.brut.para.nbcomp;
    Recap{2,8}=predSV2.predVCN.brut.para.nbcomp;
    %     Recap{2,9}=predSV2.predGAPLS.brut.para.nbcomp;
    
    Recap{3,1}='NbWl';
    %     Recap{3,2}=size(predSV.predVIP.brut.para.P,1);
    %     Recap{3,3}=size(predSV2.predTP.brut.para.P,1);
    %     Recap{3,4}=size(predSV2.predUVE.brut.para.P,1);
    %     Recap{3,5}=size(predSV.predCARS.brut.para.P,1);
    Recap{3,6}=size(predSV2.predRF.brut.para.P,1);
    %     Recap{3,7}=size(predSV.predIRF.brut.para.P,1);
    Recap{3,8}=size(predSV2.predVCN.brut.para.P,1);
    %     Recap{3,9}=size(predSV2.predGAPLS.brut.para.P,1);
    
    Recap{4,1}='R²cal';
    %     Recap{4,2}=predSV.predVIP.brut.coeff.R2c;
    %     Recap{4,3}=predSV2.predTP.brut.coeff.R2c;
    %     Recap{4,4}=predSV2.predUVE.brut.coeff.R2c;
    %     Recap{4,5}=predSV.predCARS.brut.coeff.R2c;
    Recap{4,6}=predSV2.predRF.brut.coeff.R2c;
    %     Recap{4,7}=predSV.predIRF.brut.coeff.R2c;
    Recap{4,8}=predSV2.predVCN.brut.coeff.R2c;
    %     Recap{4,9}=predSV2.predGAPLS.brut.coeff.R2c;
    
    Recap{5,1}='R²cv';
    %     Recap{5,2}=predSV.predVIP.brut.coeff.R2p;
    %     Recap{5,3}=predSV2.predTP.brut.cv.r2cv;
    %     Recap{5,4}=predSV2.predUVE.brut.cv.r2cv;
    %     Recap{5,5}=predSV.predCARS.brut.coeff.R2p;
    Recap{5,6}=predSV2.predRF.brut.cv.r2cv;
    %     Recap{5,7}=predSV.predIRF.brut.coeff.R2p;
    Recap{5,8}=predSV2.predVCN.brut.cv.r2cv;
    %     Recap{5,9}=predSV2.predGAPLS.brut.cv.r2cv;
    
    Recap{6,1}='R²val';
    %     Recap{6,2}=predSV.predVIP.brut.coeff.R2p;
    %     Recap{6,3}=predSV2.predTP.brut.coeff.R2p;
    %     Recap{6,4}=predSV2.predUVE.brut.coeff.R2p;
    %     Recap{6,5}=predSV.predCARS.brut.coeff.R2p;
    Recap{6,6}=predSV2.predRF.brut.coeff.R2p;
    %     Recap{6,7}=predSV.predIRF.brut.coeff.R2p;
    Recap{6,8}=predSV2.predVCN.brut.coeff.R2p;
    %     Recap{6,9}=predSV2.predGAPLS.brut.coeff.R2p;
    
    Recap{7,1}='RMSEP';
    %     Recap{7,2}=predSV.predVIP.brut.coeff.RMSEP;
    %     Recap{7,3}=predSV2.predTP.brut.coeff.RMSEP;
    %     Recap{7,4}=predSV2.predUVE.brut.coeff.RMSEP;
    %     Recap{7,5}=predSV.predCARS.brut.coeff.RMSEP;
    Recap{7,6}=predSV2.predRF.brut.coeff.RMSEP;
    %     Recap{7,7}=predSV.predIRF.brut.coeff.RMSEP;
    Recap{7,8}=predSV2.predVCN.brut.coeff.RMSEP;
    %     Recap{7,9}=predSV2.predGAPLS.brut.coeff.RMSEP;
    
    predSV2.Recap=Recap;
    assignin('base','predSV2',predSV2);
    
    % Affichage des zones influentes
    if mean(Recap{5,6}>0.5)||mean(Recap{5,8}>0.5)%||Recap{5,9}>0.5%|Recap{5,5}>0.5||Recap{5,7}>0.5||mean(Recap{5,3}>0.5)||mean(Recap{5,4}>0.5)||
        fig = figure;
        fig.PaperPositionMode = 'auto';
        fig.InvertHardcopy = 'off';
        set(fig,'units','characters', 'position',[0 0 a_characters(1) a_characters(2)], 'name', 'Influential areas');
        i=0;
        j=1;
        %         if Recap{5,2}>0.5
        %             plot(wl(SV.idxVIPopt),ones(1,length(SV.idxVIPopt))+i,'*','Markersize',10), hold on
        %             i=i+0.1;
        %             legvl{j}=strcat('VIP : ', num2str(round(Recap{5,2},2)), '(', num2str(Recap{3,2}), '), ', num2str(Recap{7,2}));
        %             j=j+1;
        %         end
        %         if mean(Recap{6,3})>0.3
        %             plot(wl(predSV2.TPidx),ones(1,length(predSV2.TPidx))+i,'*','Markersize',10), hold on
        %             i=i+0.1;
        %             legvl{j}=strcat('TP : ', num2str(round(mean(Recap{6,3}),2)), '(', num2str(Recap{3,3}), ')');
        %             j=j+1;
        %         end
        %         if mean(Recap{6,4})>0.3
        %             plot(wl(predSV2.UVEidx),ones(1,length(predSV2.UVEidx))+i,'*','Markersize',10), hold on
        %             i=i+0.1;
        %             legvl{j}=strcat('UVE : ', num2str(round(mean(Recap{6,4}),2)), '(', num2str(Recap{3,4}), ')');
        %             j=j+1;
        %         end
        %         if Recap{5,5}>0.5
        %             plot(wl(SV.idxCARSopt),ones(1,length(SV.idxCARSopt))+i,'*','Markersize',10), hold on
        %             i=i+0.1;
        %             legvl{j}=strcat('CARS : ', num2str(round(Recap{5,5},2)), '(', num2str(Recap{3,5}), '), ', num2str(Recap{7,5}));
        %             j=j+1;
        %         end
        if mean(Recap{6,6})>0.3
            plot(wl(predSV2.RFidx),ones(1,length(predSV2.RFidx))+i,'*','Markersize',10), hold on
            i=i+0.1;
            legvl{j}=strcat('RF : ', num2str(round(mean(Recap{6,6}),2)), '(', num2str(Recap{3,6}), ')');
            j=j+1;
        end
        %         if Recap{5,7}>0.5
        %             plot(wl(SV.idxiRF(1:Recap{3,7})),ones(1,length(SV.idxiRF(1:Recap{3,7})))+i,'*','Markersize',10), hold on
        %             i=i+0.1;
        %             legvl{j}=strcat('iRF : ', num2str(round(Recap{5,7},2)), '(', num2str(Recap{3,7}), '), ', num2str(Recap{7,7}));
        %             j=j+1;
        %         end
        if mean(Recap{6,8})>0.3
            plot(wl(predSV2.VCNidx),ones(1,length(predSV2.VCNidx))+i,'*','Markersize',10), hold on
            i=i+0.1;
            legvl{j}=strcat('VCN : ', num2str(round(mean(Recap{6,8}),2)), '(', num2str(Recap{3,8}), ')');
            j=j+1;
        end
        %         if Recap{6,9}>0.3
        %             plot(wl(predSV2.GAPLSidx),ones(1,length(predSV2.GAPLSidx))+i,'*','Markersize',10), hold on
        %             legvl{j}=strcat('GA-PLS : ', num2str(round(Recap{6,9},2)), '(', num2str(Recap{3,9}), '), ', num2str(Recap{7,9}));
        %         end
        hold off
        grid on
        legend(legvl)
        xlabel('Wavelength (nm)')
        set(gca,'FontSize',16)
        
        % Création d'une structure modèle
        Modelc = questdlg('Do you want to create a model ?','Model creation','No','Yes','Yes');
        if strcmp(Modelc,'Yes')
            algo={'VIP','TP','UVE','CARS','RF','iRF','VCN','GA-PLS'};
            Modela = listdlg('PromptString','With which algorithm do you want to work?','SelectionMode','single', 'ListString',algo);
            
            if strcmp(pseudoabs,'Yes')
                Model.Type='PseudoAbsorbance';
            else
                Model.Type='Reflectance';
            end
            Model.label=Yn;
            Model.SelectWave=algo{Modela};
            Model.cv=cvt;
            Model.wli=wl;
            
            if Modela==1
                %                 Model.Scal{1,1}=predSV2.predVIP.Scal;
                %                 Model.Ycal=predSV2.predVIP.Ycal;
                %                 Model.Sval{1,1}=predSV2.predVIP.Sval;
                %                 Model.Yval=predSV2.predVIP.Yval;
                %                 Model.Pret=plsPredOpt.Pretopt{1, 1}  ;
                %                 Model.Pret1=DataPLS.Pret;
                %                 Model.Pret1a=[DataPLS.Pret1; DataPLS.Pret2];
                %                 Model.wl=predSV2.VIPwl;
                %                 Model.B=predSV2.predVIP.brut.para.B;
                %                 Model.R2cal=predSV2.predVIP.brut.coeff.R2c;
                %                 Model.R2cv=predSV2.predVIP.brut.cv.r2cv;
                %                 Model.R2val=predSV2.predVIP.brut.coeff.R2p;
                %                 Model.SEC=predSV2.predVIP.brut.coeff.SEC;
                %                 Model.SECV=predSV2.predVIP.brut.cv.SECV;
                %                 Model.SEP=predSV2.predVIP.brut.coeff.SEP;
                %                 Model.RMSEP=predSV2.predVIP.brut.coeff.RMSEP;
                %                 Model.model=predSV2.predVIP.brut ;
                %                 Model.modelsansselectvar=PRED{1,opt(1)};
            else if Modela==2
                    %                     Model.Scal{1,1}=predSV2.predTP.Scal;
                    %                     Model.Ycal=predSV2.predTP.Ycal;
                    %                     Model.Sval{1,1}=predSV2.predTP.Sval;
                    %                     Model.Yval=predSV2.predTP.Yval;
                    %                     Model.Pret=plsPredOpt.Pretopt{1, 1}  ;
                    %                     Model.Pret1=DataPLS.Pret;
                    %                     Model.Pret1a=[DataPLS.Pret1; DataPLS.Pret2];
                    %                     Model.wl=predSV2.TPwl;
                    %                     Model.B=predSV2.predTP.brut.para.B;
                    %                     Model.R2cal=predSV2.predTP.brut.coeff.R2c;
                    %                     Model.R2cv=predSV2.predTP.brut.cv.r2cv;
                    %                     Model.R2val=predSV2.predTP.brut.coeff.R2p;
                    %                     Model.SEC=predSV2.predTP.brut.coeff.SEC;
                    %                     Model.SECV=predSV2.predTP.brut.cv.SECV;
                    %                     Model.SEP=predSV2.predTP.brut.coeff.SEP;
                    %                     Model.RMSEP=predSV2.predTP.brut.coeff.RMSEP;
                    %                     Model.model=predSV2.predTP.brut ;
                    %                     Model.modelsansselectvar=PRED{1,opt(1)};
                else if Modela==3
                        %                         Model.Scal{1,1}=predSV2.predUVE.Scal;
                        %                         Model.Ycal=predSV2.predUVE.Ycal;
                        %                         Model.Sval{1,1}=predSV2.predUVE.Sval;
                        %                         Model.Yval=predSV2.predUVE.Yval;
                        %                         Model.Pret=plsPredOpt.Pretopt{1, 1}  ;
                        %                         Model.Pret1=DataPLS.Pret;
                        %                         Model.Pret1a=[DataPLS.Pret1; DataPLS.Pret2];
                        %                         Model.wl=predSV2.UVEwl;
                        %                         Model.B=predSV2.predUVE.brut.para.B;
                        %                         Model.R2cal=predSV2.predUVE.brut.coeff.R2c;
                        %                         Model.R2cv=predSV2.predUVE.brut.cv.r2cv;
                        %                         Model.R2val=predSV2.predUVE.brut.coeff.R2p;
                        %                         Model.SEC=predSV2.predUVE.brut.coeff.SEC;
                        %                         Model.SECV=predSV2.predUVE.brut.cv.SECV;
                        %                         Model.SEP=predSV2.predUVE.brut.coeff.SEP;
                        %                         Model.RMSEP=predSV2.predUVE.brut.coeff.RMSEP;
                        %                         Model.model=predSV2.predUVE.brut ;
                        %                         Model.modelsansselectvar=PRED{1,opt(1)};
                    else if Modela==4
                            %                             Model.Scal{1,1}=predSV2.predCARS.Scal;
                            %                             Model.Ycal=predSV2.predCARS.Ycal;
                            %                             Model.Sval{1,1}=predSV2.predCARS.Sval;
                            %                             Model.Yval=predSV2.predCARS.Yval;
                            %                             Model.Pret=plsPredOpt.Pretopt{1, 1}  ;
                            %                             Model.Pret1=DataPLS.Pret;
                            %                             Model.Pret1a=[DataPLS.Pret1; DataPLS.Pret2];
                            %                             Model.wl=predSV2.CARSwl;
                            %                             Model.B=predSV2.predCARS.brut.para.B;
                            %                             Model.R2cal=predSV2.predCARS.brut.coeff.R2c;
                            %                             Model.R2val=predSV2.predCARS.brut.coeff.R2p;
                            %                             Model.SEC=predSV2.predCARS.brut.coeff.SEC;
                            %                             Model.SEP=predSV2.predCARS.brut.coeff.SEP;
                            %                             Model.RMSEP=predSV2.predCARS.brut.coeff.RMSEP;
                            %                             Model.model=predSV2.predCARS.brut ;
                        else if Modela==5
                                Model.Scal{1,1}=predSV2.predRF.Scal;
                                Model.Ycal=predSV2.predRF.Ycal;
                                Model.Sval{1,1}=predSV2.predRF.Sval;
                                Model.Yval=predSV2.predRF.Yval;
                                Model.Pret=plsPredOpt.Pretopt{1, 1};
                                Model.Pret1=DataPLS.Pret;
                                Model.Pret1a=[DataPLS.Pret1; DataPLS.Pret2];
                                Model.wl=predSV2.RFwl;
                                Model.B=predSV2.predRF.brut.para.B;
                                Model.R2cal=predSV2.predRF.brut.coeff.R2c;
                                Model.R2cv=predSV2.predRF.brut.cv.r2cv;
                                Model.R2val=predSV2.predRF.brut.coeff.R2p;
                                Model.SEC=predSV2.predRF.brut.coeff.SEC;
                                Model.SECV=predSV2.predRF.brut.cv.SECV;
                                Model.SEP=predSV2.predRF.brut.coeff.SEP;
                                Model.RMSEP=predSV2.predRF.brut.coeff.RMSEP;
                                Model.model=predSV2.predRF.brut ;
                                Model.modelsansselectvar=PRED{1,opt(1)};
                            else if Modela==6
                                    %                                     Model.Scal{1,1}=predSV2.predIRF.Scal;
                                    %                                     Model.Ycal=predSV2.predIRF.Ycal;
                                    %                                     Model.Sval{1,1}=predSV2.predIRF.Sval;
                                    %                                     Model.Yval=predSV2.predIRF.Yval;
                                    %                                     Model.Pret=plsPredOpt.Pretopt{1, 1}  ;
                                    %                                     Model.Pret1=DataPLS.Pret;
                                    %                                     Model.Pret1a=[DataPLS.Pret1; DataPLS.Pret2];
                                    %                                     Model.wl=predSV2.IRFwl;
                                    %                                     Model.B=predSV2.predIRF.brut.para.B;
                                    %                                     Model.R2cal=predSV2.predIRF.brut.coeff.R2c;
                                    %                                     Model.R2val=predSV2.predIRF.brut.coeff.R2p;
                                    %                                     Model.SEC=predSV2.predIRF.brut.coeff.SEC;
                                    %                                     Model.SEP=predSV2.predIRF.brut.coeff.SEP;
                                    %                                     Model.RMSEP=predSV2.predIRF.brut.coeff.RMSEP;
                                    %                                     Model.model=predSV2.predIRF.brut ;
                                else if Modela==7
                                        Model.Scal{1,1}=predSV2.predVCN.Scal;
                                        Model.Ycal=predSV2.predVCN.Ycal;
                                        Model.Sval{1,1}=predSV2.predVCN.Sval;
                                        Model.Yval=predSV2.predVCN.Yval;
                                        Model.Pret=plsPredOpt.Pretopt{1, 1}  ;
                                        Model.Pret1=DataPLS.Pret;
                                        Model.Pret1a=[DataPLS.Pret1; DataPLS.Pret2];
                                        Model.wl=predSV2.VCNwl;
                                        Model.B=predSV2.predVCN.brut.para.B;
                                        Model.R2cal=predSV2.predVCN.brut.coeff.R2c;
                                        Model.R2cv=predSV2.predVCN.brut.cv.r2cv;
                                        Model.R2val=predSV2.predVCN.brut.coeff.R2p;
                                        Model.SEC=predSV2.predVCN.brut.coeff.SEC;
                                        Model.SECV=predSV2.predVCN.brut.cv.SECV;
                                        Model.SEP=predSV2.predVCN.brut.coeff.SEP;
                                        Model.RMSEP=predSV2.predVCN.brut.coeff.RMSEP;
                                        Model.model=predSV2.predVCN.brut ;
                                        Model.modelsansselectvar=PRED{1,opt(1)};
                                    else if Modela==8
                                            %                                             Model.Scal{1,1}=predSV2.predGAPLS.Scal;
                                            %                                             Model.Ycal=predSV2.predGAPLS.Ycal;
                                            %                                             Model.Sval{1,1}=predSV2.predGAPLS.Sval;
                                            %                                             Model.Yval=predSV2.predGAPLS.Yval;
                                            %                                             Model.Pret=plsPredOpt.Pretopt{1, 1}  ;
                                            %                                             Model.Pret1=DataPLS.Pret;
                                            %                                             Model.Pret1a=[DataPLS.Pret1; DataPLS.Pret2];
                                            %                                             Model.wl=predSV2.GAPLSwl;
                                            %                                             Model.B=predSV2.predGAPLS.brut.para.B;
                                            %                                             Model.R2cal=predSV2.predGAPLS.brut.coeff.R2c;
                                            %                                             Model.R2cv=predSV2.predGAPLS.brut.cv.r2cv;
                                            %                                             Model.R2val=predSV2.predGAPLS.brut.coeff.R2p;
                                            %                                             Model.SEC=predSV2.predGAPLS.brut.coeff.SEC;
                                            %                                             Model.SECV=predSV2.predGAPLS.brut.cv.SECV;
                                            %                                             Model.SEP=predSV2.predGAPLS.brut.coeff.SEP;
                                            %                                             Model.RMSEP=predSV2.predGAPLS.brut.coeff.RMSEP;
                                            %                                             Model.model=predSV2.predGAPLS.brut ;
                                            %                                             Model.modelsansselectvar=PRED{opt,1};
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            Model.Recap{1,2}='All';
            Model.Recap{1,3}='Selected';
            
            Model.Recap{2,1}='Nb wl';
            Model.Recap{2,2}=length(wl);
            Model.Recap{2,3}=length(Model.wl);
            
            Model.Recap{3,1}='Nb vl';
            Model.Recap{3,3}=Model.model.para.nbcomp;
            Model.Recap{3,2}=Model.modelsansselectvar.recap{Pretoptg+1,2};
            
            Model.Recap{4,1}='SEC';
            Model.Recap{4,3}=Model.SEC;
            Model.Recap{4,2}=Model.modelsansselectvar.recap{Pretoptg+1,3};
            
            Model.Recap{5,1}='R2c';
            Model.Recap{5,3}=Model.R2cal;
            Model.Recap{5,2}=Model.modelsansselectvar.recap{Pretoptg+1,4};
            
            Model.Recap{6,1}='SECV';
            Model.Recap{6,3}=Model.SECV;
            Model.Recap{6,2}=Model.modelsansselectvar.recap{Pretoptg+1,6};
            
            Model.Recap{7,1}='R2cv';
            Model.Recap{7,3}=Model.R2cv;
            Model.Recap{7,2}=Model.modelsansselectvar.recap{Pretoptg+1,7};
            
            Model.Recap{8,1}='SEP';
            Model.Recap{8,3}=Model.SEP;
            Model.Recap{8,2}=Model.modelsansselectvar.recap{Pretoptg+1,9};
            
            Model.Recap{9,1}='RMSEP';
            Model.Recap{9,3}=Model.model.coeff.RMSEP;
            Model.Recap{9,2}=Model.modelsansselectvar.recap{Pretoptg+1,10};
            
            Model.Recap{10,1}='R2p';
            Model.Recap{10,3}=Model.R2val;
            Model.Recap{10,2}=Model.modelsansselectvar.recap{Pretoptg+1,11};
            
            assignin('base','Modelvs',Model);
        end
    end
    
end

end