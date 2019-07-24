function [CartoY, CartoYlabel,CartoIC,CartoIP,Correl, Yc]=ApplyModel(M,d,wl,Model,RGB,figi,Yref,dref,pas,pas2)
% Function to apply the model created with the function CreateModel.
%
% INPUT: 
%       M: HSI datacube (n*m*p)
%       d : associated depth (1*m)(mm)
%       wl: associated wavelength (1*p)
%       Model: Structure containing the model parameters
%       RGB: Image of the sample
%       figi: display figure with value different from 0
%       Yref (optionnal): reference values of the variable predicted to
%           compare with prediction
%       dref (optionnal): depth of the reference values (mm)
%       pas (optionnal): sampling thickness (mm)
%       pas2 (optionnal): to say if the sampling go after or before the dref values
%                   ex: dref = 0.5, pas = 0.5, pas2= 'Max', sampling depth goes from 0 to 0.5
% OUTPUT:
%       CartoY: Predicted abundance map
%       CartoYlabel: Name of the variable predicted
%       CartoIC: Confidence interval map
%       CartoIP: Prediction interval map
%       Correl (only if ref available): correlation between predicted and
%           reference data
%       Yc (only if ref available): subsampling predicted values at the
%           reference resolution
%
% This is the software from the papers:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Jacq, K., Perrette, Y., Fanget, B., Sabatier, P., Coquin, D., 
% Martinez-Lamas, R., Debret, M., Arnaud, F., 2019. High-resolution 
% prediction of organic matter concentration with hyperspectral imaging 
% on a sediment core. Sci. Total Environ. 663, 236–244. 
% https://doi.org/10.1016/j.scitotenv.2019.01.320

% Please cite our papers if you use our code for your research.

% Verification of inputs
if nargin>6&&nargin<9
    pas=0;
    pas2='';
end

if isfield(Model,'wli')==0
    Model.wli=Model.wl;
end

% Conversion in absorbance if necessary
if strcmp(Model.Type,'PseudoAbsorbance')
    if isa(M,'cell')==0
        S=reshape(M,size(M,1)*size(M,2),size(M,3));
        [a,b]=find(S==0);
        S(a,b)=eps;
        S=log(1./S);
        M=reshape(S,size(M,1),size(M,2),size(M,3));
    else if isa(M,'cell')==1
            for f=1:size(M,1)
                S=reshape(M{f,1},size(M{f,1},1)*size(M{f,1},2),size(M{f,1},3));
                [a,b]=find(S==0);
                S(a,b)=eps;
                S=log(1./S);
                M{f,1}=reshape(S,size(M{f,1},1),size(M{f,1},2),size(M{f,1},3));
            end
        end
    end
end

if nargin>5
    % GUI at the screen resolution
    set(gcf,'Units','pixels')
    size_pixels=get(gcf,'Position');
    set(gcf,'Units','characters')
    size_characters=get(gcf,'Position');
    f=size_pixels(3:4)./size_characters(3:4);
    Pix_SS = get(0,'screensize');
    a=Pix_SS(:,3:4);
    a_characters=a./f;
    close all
end

% Unfolding datacube
if isa(M,'cell')==0
    if size(M,3)==1
        S=M;
    else
        S=reshape(M,size(M,1)*size(M,2),size(M,3));
    end
else if isa(M,'cell')==1
        for f=1:size(M,1)
            if size(M{f,1},3)==1
                S{f,1}=M{f,1};
            else
                S{f,1}=reshape(M{f,1},size(M{f,1},1)*size(M{f,1},2),size(M{f,1},3));
            end
        end
    end
end

% First preprocessing
if strcmp(Model.Pret1,'autoscaling')
    if isa(M,'cell')==0
        Scalp=Model.Scal{1,1};
        if size(S,2)>length(Model.wli)
            % Selected wavelengths (first range ~without noisy bands)
            wli=zeros(1,size(Model.wl,2));
            for j=1:size(Model.wli,2)
                [~, wlii]=find(abs(wl-Model.wli(j))==min(abs(wl-Model.wli(j))));
                wli(j)=wlii(1);
            end
            S=S(:,wli);
        else if size(S,2)>length(Model.wli)
                error('S don''t have a good size')
            end
        end
        
        if size(Scalp,2)>length(Model.wli)
            error('Scal don''t have a good size')
        else if size(S,2)>length(Model.wli)
                error('S don''t have a good size')
            end
        end
        
        [Sp,~,~]=pretreat(S,Model.Pret1,Model.Pret1a(1,:),Model.Pret1a(2,:));
        
    else if isa(M,'cell')==1
            for f=1:size(M,1)
                Scalp{f,1}=Model.Scal{1,1}{f,1};
                if size(S{f,1},2)>length(Model.wli)
                    % Selected wavelengths (first range ~without noisy bands)
                    wli=zeros(1,size(Model.wl,2));
                    for j=1:size(Model.wli{f,1},2)
                        [~, wlii]=find(abs(wl{f,1}-Model.wli{f,1}(j))==min(abs(wl{f,1}-Model.wli{f,1}(j))));
                        wli(j)=wlii(1);
                    end
                    S{f,1}=S{f,1}(:,wli);
                else if size(S,2)>length(Model.wli)
                        error('S don''t have a good size')
                    end
                end
                
                if size(Scalp,2)>length(Model.wli{f,1})
                    error('Scal don''t have a good size')
                else if size(S{f,1},2)>length(Model.wli{f,1})
                        error('S don''t have a good size')
                    end
                end
                
                [Sp{f,1},~,~]=pretreat(S{f,1},Model.Pret1,Model.Pret1a{f,1},Model.Pret1a{f,2});
            end
        end
    end
else
    if isa(M,'cell')==0
        Scalp=Model.Scal{1,1};
        if size(S,2)>length(Model.wli)
            % Selected wavelengths (first range ~without noisy bands)
            wli=zeros(1,size(Model.wl,2));
            for j=1:size(Model.wli,2)
                [~, wlii]=find(abs(wl-Model.wli(j))==min(abs(wl-Model.wli(j))));
                wli(j)=wlii(1);
            end
            S=S(:,wli);
        else if size(S,2)>length(Model.wli)
                error('S don''t have a good size')
            end
        end
        
        if size(Scalp,2)>length(Model.wli)
            error('Scal don''t have a good size')
        else if size(S,2)>length(Model.wli)
                error('S don''t have a good size')
            end
        end
        
        Scalp=Model.Scal{1,1};
        Sp=S;
        
    else if isa(M,'cell')==1
            for f=1:size(M,1)
                Scalp{f,1}=Model.Scal{1,1}{f,1};
                if size(S{f,1},2)>length(Model.wli{f,1})
                    % Selected wavelengths (first range ~without noisy bands)
                    wli=zeros(1,size(Model.wl,2));
                    for j=1:size(Model.wli{f,1},2)
                        [~, wlii]=find(abs(wl{f,1}-Model.wli{f,1}(j))==min(abs(wl{f,1}-Model.wli{f,1}(j))));
                        wli(j)=wlii(1);
                    end
                    S{f,1}=S{f,1}(:,wli);
                else if size(S{f,1},2)>length(Model.wli{f,1})
                        error('S don''t have a good size')
                    end
                end
                
                if size(Scalp,2)>length(Model.wli{f,1})
                    error('Scal don''t have a good size')
                else if size(S{f,1},2)>length(Model.wli{f,1})
                        error('S don''t have a good size')
                    end
                end
                
                Scalp{f,1}=Model.Scal{1,1}{f,1};
                Sp{f,1}=S{f,1};
            end
        end
    end
end

if isa(M,'cell')==1
    wlim=[];
    for f=1:size(M,1)
        wlim=[wlim Model.wli{f,1}];
    end
else
    wlim=Model.wli;
end

% Prediction
Ypred=zeros(size(S,1),size(Model.Ycal,2));
if size(Model.Ycal,2)==1
    % Second preprocessing
    wli=[];
    if isa(M,'cell')==0
        if strcmp(Model.Pret,'Raw')
            % Selected wavelengths (second range)
            wli=zeros(1,size(Model.wl,2));
            for j=1:size(Model.wl,2)
                [~, wlii]=find(abs(wlim-Model.wl(j))==min(abs(wlim-Model.wl(j))));
                wli(j)=wlii(1);
            end
            Spret=Sp(:,wli);
            Scalpret=Scalp;
        else
            Spret = AllPret(Sp,Model.wli,Model.Pret);
            Scalpret = Scalp;
            
            % Selected wavelengths (second range)
            wli=zeros(1,size(Model.wl,2));
            for j=1:size(Model.wl,2)
                [~, wlii]=find(abs(wlim-Model.wl(j))==min(abs(wlim-Model.wl(j))));
                wli(j)=wlii(1);
            end
            Spret=Spret(:,wli);
        end
        
    else if isa(M,'cell')==1
            Spret=[];Scalpret=[];wli=[];
            for f=1:size(M,1)
                if strcmp(Model.Pret{1,f},'Raw')
                    Spreti=Sp{f,1};
                    Scalpreti=Scalp{f,1};
                else
                    Spreti= AllPret(Sp{f,1},Model.wli{f,1},Model.Pret{1,f});
                    Scalpreti = Scalp{f,1};
                    
                end
                Spret=[Spret Spreti];
                Scalpret=[Scalpret Scalpreti];
            end
            for j=1:size(Model.wl,2)
                % Selected wavelengths (second range)
                [~, wlii]=find(abs(wlim-Model.wl(j))==min(abs(wlim-Model.wl(j))));
                wli(j)=wlii(1);
            end
            Spret=Spret(:,wli);
        end
    end
    
    if size(Scalpret,2)>length(Model.wl)
        for j=1:size(Model.wl,2)
            % Selected wavelengths (second range)
            [~, wlii]=find(abs(wlim-Model.wl(j))==min(abs(wlim-Model.wl(j))));
            wli(j)=wlii(1);
        end
        Scalpret=Scalpret(:,wli);
    end
else
    % Selected wavelengths (second range)
    wli=zeros(1,size(Model.wl,2));
    for j=1:size(Model.wl,2)
        [~, wlii]=find(abs(wlim-Model.wl(j))==min(abs(wlim-Model.wl(j))));
        wli(j)=wlii(1);
    end
    Spret=Sp(:,wli);
    Scalpret=Scalp;
    
    if size(Scalpret,2)>length(Model.wl)
        for j=1:size(Model.wl,2)
            % Selected wavelengths (second range)
            [~, wlii]=find(abs(wlim-Model.wl(j))==min(abs(wlim-Model.wl(j))));
            wli(j)=wlii(1);
        end
        Scalpret=Scalpret(:,wli);
    end
    
end
if isstruct(Scalp)
    Ypred=Decentrerval(Centrerval(Spret,Scalpret.Raw)*Model.B,Model.Ycal);
else
    Ypred=Decentrerval(Centrerval(Spret,Scalpret)*Model.B,Model.Ycal);
end

% Folding of the prediction
if size(Model.Ycal,2)>1
    if size(M,3)>1
        CartoY=reshape(Ypred,size(M,1),size(M,2),size(Model.Ycal,2));
        Cartoms=sum(abs(squeeze(mean(CartoY,1))));
        CartoY=CartoY(:,:,Cartoms>0);
        CartoYlabel=Model.label(Cartoms>0);
    else
        CartoY=reshape(Ypred,size(M,1),size(Model.Ycal,2));
        Cartoms=CartoY;
        CartoY=CartoY(:,Cartoms>0);
        CartoYlabel=Model.label(Cartoms>0);
    end
else
    if size(M,3)>1
        CartoY=reshape(Ypred,size(M,1),size(M,2));
        CartoYlabel=Model.label;
    else
        CartoY=Ypred;
        CartoYlabel=Model.label;
    end
end

% Calculation of confidance and prediction intervals
[CartoIC,CartoIP]=IntervalPredModel(Model,CartoY);

% Saving
assignin('base', 'CartoY', CartoY);
assignin('base', 'CartoYmlabel', CartoYlabel);
assignin('base', 'CartoIC', CartoIC);
assignin('base', 'CartoIP', CartoIP);

% Display
if nargin>5
    if figi>0
        for k=1:size(CartoY,3)
            if nargin>6
                fig = figure;
                fig.PaperPositionMode = 'auto';
                fig.InvertHardcopy = 'off';
                set(fig,'units','characters', 'position',[0 0 a_characters(1) a_characters(2)], 'name', 'Carte d''abondance');
                
                CartoYi=squeeze(CartoY(:,:,k));
                ha(1)=subplot(4,80,1:80);imagesc(d,d(1:size(CartoYi,1)),CartoYi);colormap(jet);colorbar,title('Carte d''abondance')
                ylabel('Largeur (cm)')
                caxis([mean(CartoYi(:))-3*std(CartoYi(:)) mean(CartoYi(:))+3*std(CartoYi(:))])
                set(findobj('type','axes'),'fontsize',16)
                if size(CartoY,3)>1
                    title(Model.label(k))
                else
                    title(Model.label)
                end

                ha(2)=subplot(4,80,81:157);imagesc(d,d(1:size(CartoYi,1)),RGB*1.4);
                ylabel('Largeur (cm)')
                set(findobj('type','axes'),'fontsize',16)
                
                ha(3)=subplot(4,80,161:237);
                Cartomin=zeros(1,size(CartoY,2));Cartomax=zeros(1,size(CartoY,2));
                for i=1:size(CartoY,2)
                    CartoYi= squeeze(CartoY(:,i,k));
                    Cartomax(i)=max(CartoYi(CartoYi<mean(CartoYi)+2*std(CartoYi)));
                    Cartomin(i)=min(CartoYi(CartoYi>mean(CartoYi)-2*std(CartoYi)));
                end
                h=area(d,[Cartomin; Cartomax-Cartomin]','LineStyle','none');
                h(1).FaceColor = [1 1 1];
                h(2).FaceColor = [0 1 0];
                if nargin>6
                    hold on
                    h1=plot(dref,Yref(:,k),'*','LineWidth',3);
                end
                set(gca,'Layer','top')
                
                if size(CartoY,3)>1
                    ylabel(Model.label(k))
                else
                    ylabel(Model.label)
                end
                xlim([min(d) max(d)])
                ylim([min([mean(Cartomin(:))-3*std(Cartomin(:)) Yref(:,k)']) max([mean(Cartomax(:))+3*std(Cartomax(:)) Yref(:,k)'])])
                legend([h(2) h1],'Prédit','Observé')
                grid on
                set(findobj('type','axes'),'fontsize',16)
                
                CartoICi=squeeze(CartoIC(:,:,k));
                ha(4)=subplot(4,80,241:320);
                imagesc(d,d(1:size(CartoYi,1)),CartoICi)
                colorbar
                colormap(jet)
                caxis([mean(CartoICi(:))-3*std(CartoICi(:)) mean(CartoICi(:))+3*std(CartoICi(:))])
                title(strcat('Confidence Interval map : ',num2str(Model.RMSEP(k)*2)))
                xlabel('Profondeur (cm)')
                ylabel('Largeur (cm)')
                set(findobj('type','axes'),'fontsize',16)
                
                linkaxes(ha,'x')
            else
                fig = figure;
                fig.PaperPositionMode = 'auto';
                fig.InvertHardcopy = 'off';
                set(fig,'units','characters', 'position',[0 0 a_characters(1) a_characters(2)], 'name', 'Carte d''abondance');
                
                CartoYi=squeeze(CartoY(:,:,k));
                ha(1)=subplot(3,80,1:80);imagesc(CartoYi);colormap(jet);colorbar,title('Carte d''abondance')
                xlabel('Profondeur (pixel)')
                ylabel('Largeur (pixel)')
                caxis([mean(CartoYi(:))-3*std(CartoYi(:)) mean(CartoYi(:))+3*std(CartoYi(:))])
                set(findobj('type','axes'),'fontsize',16)
                if size(CartoY,3)>1
                    title(Model.label(k))
                else
                    title(Model.label)
                end
                ha(2)=subplot(3,80,81:157);imagesc(RGB*1.4);
                xlabel('Profondeur (pixel)')
                ylabel('Largeur (pixel)')
                set(findobj('type','axes'),'fontsize',16)
                
                linkaxes(ha,'xy')
                
                subplot(3,80,161:237);
                Cartomin=zeros(1,size(CartoY,2));Cartomax=zeros(1,size(CartoY,2));
                for i=1:size(CartoY,2)
                    CartoYi= squeeze(CartoY(:,i,k));
                    Cartomax(i)=max(CartoYi(CartoYi<mean(CartoYi)+2*std(CartoYi)));
                    Cartomin(i)=min(CartoYi(CartoYi>mean(CartoYi)-2*std(CartoYi)));
                end
                h=area(d,[Cartomin; Cartomax-Cartomin]','LineStyle','none');
                h(1).FaceColor = [1 1 1];
                h(2).FaceColor = [0 1 0];
                set(gca,'Layer','top')
                xlabel('Profondeur (cm)')
                if size(CartoY,3)>1
                    ylabel(Model.label(k))
                else
                    ylabel(Model.label)
                end
                xlim([min(d) max(d)])
                ylim([min([mean(Cartomin(:))-3*std(Cartomin(:))]) max([mean(Cartomax(:))+3*std(Cartomax(:))])])
                legend([h(2)],'Prédit')
                grid on
                set(findobj('type','axes'),'fontsize',16)
            end
        end
    end
end

% Correlation with reference values
if nargin>6
    for k=1:size(CartoY,3)
        [Correl(k), Yc{k}] = CorrelCartoY(squeeze(CartoY(:,:,k)),d,Yref(:,k),dref,figi,pas,pas2);
    end
else
    Correl=[];
    Yc=[];
end

end