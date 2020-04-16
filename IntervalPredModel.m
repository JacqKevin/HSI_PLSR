function [CartoIC, CartoIP]=IntervalPredModel(Model,Carto,fig)
% Function that compute confidence interval and maps thanks to the PLSR
% model and the classification map.

if nargin<3
    fig=0;
end

if isfield(Model.model,'Pred')==1
    X=[Model.model.Pred.Ycal; Model.model.Pred.Yval];
else
    X=[Model.model.Model.Pred.Ycal; Model.model.Model.Pred.Yval];
end
X2=X(:,size(X,2)/2+1:end);
X1=X(:,1:size(X,2)/2);

n=size(X2,1);
t=tinv(0.95,n);
s2=std(X2);
s1=std(X1);
xm=mean(X1);
for i=1:length(s1)
    intervx=0:2*round(max(X1(:,i)));
%     intervy=t*sqrt(1+1/n+(intervx-xm(i)).^2/n*(s1(i).^2));
    intervy=t*s1(i)*sqrt(1+1/n+(intervx-xm(i)).^2/(sum(X1(:,i).^2)-1/n*(sum(X1(:,i)).^2)));
    xmc=t*sqrt(1+1/n+(xm(i)-xm(i)).^2/n*(s1(i).^2));
    if isfield(Model.model,'Pred')==1
        unc=2*Model.model.coeff.RMSEP(i);
    else
        unc=2*Model.model.Model.coeff.RMSEP(i);
    end
    c=xmc/unc;
    intervyc=intervy/c;
    corrX=corr(X1(:,i),X2(:,i));
    fitobject = fit(X1(:,i),X2(:,i),'poly1');
    txt={strcat(['n=' num2str(n)]);
        strcat(['R²=' num2str(corrX,2)]);
        strcat(['Bias=' num2str(fitobject.p2,2)]);
        strcat(['Slope=' num2str(fitobject.p1,2)]);
        strcat(['Um=' num2str(xmc,3)]);
        strcat(['Umc=' num2str(unc,3)])};
    
    if fig>0
        figure;
        h1=plot(X1(:,i),X2(:,i),'b.','Markersize',20);
        hold on
        h2=plot(intervx,intervx,'k','Linewidth',2);
        hold on
        h3=plot(intervx,intervx+intervy,'r--','Linewidth',2);
        hold on
        h4=plot(intervx,intervx-intervy,'r--','Linewidth',2);
        hold on
        h5=plot(intervx,intervx+intervyc,'m--','Linewidth',2);
        hold on
        h6=plot(intervx,intervx-intervyc,'m--','Linewidth',2);
        hold on
        h7=plot([0 intervx(end)],[fitobject.p2 intervx(end)*fitobject.p1+fitobject.p2],'--','Linewidth',2);
        grid on
        xlabel('Observed (wt%)')
        ylabel('Predicted (wt%)')
        legend([h3 h5],'Confidence interval','Corrected confidence interval')
        text(2*round(max(X1(:,i)))-5,xm(i)-2,txt,'Fontsize',14,'fontweight','bold')
        set(gca,'Fontsize',14,'fontweight','bold')
        if isstruct(Model.label)
            title(Model.label{i})
        else
            title(Model.label)
        end
    end
    
    if nargin>1
        C=reshape(squeeze(Carto(:,:,i)),size(Carto,1)*size(Carto,2),1);
        Cdc=t*s1(i)*sqrt(1/n+(C-xm(i)).^2/(sum(X1(:,i).^2)-1/n*(sum(X1(:,i))^2)));
        Cdp=t*s1(i)*sqrt(1+1/n+(C-xm(i)).^2/(sum(X1(:,i).^2)-1/n*(sum(X1(:,i))^2)));
%         Cd=t*sqrt(1+1/n+(C-xm(i)).^2/n*(s1(i).^2))/c;
        
        if fig>0
            figure;
            ha(1)=subplot(311);
            imagesc(squeeze(Carto(:,:,i)))
            colormap(jet)
            colorbar
            caxis([mean(C)-4*std(C) mean(C)+4*std(C)])
            if isstruct(Model.label)
                title(strcat('Concentration map',Model.label{i}))
            else
                title(strcat('Concentration map',Model.label))
            end
            ha(2)=subplot(312);
            imagesc(reshape(Cdc,size(Carto,1),size(Carto,2)))
            colormap(jet)
            colorbar
            caxis([unc mean(Cdc)+4*std(Cdc)])
            if isstruct(Model.label)
                title(strcat('Confidence Interval map : ',Model.label{i}))
            else
                title(strcat('Confidence Interval map : ',Model.label))
            end
            ha(3)=subplot(313);
            imagesc(reshape(Cdp,size(Carto,1),size(Carto,2)))
            colormap(jet)
            colorbar
            caxis([unc mean(Cdp)+4*std(Cdp)])
            linkaxes(ha,'xy')
            if isstruct(Model.label)
                title(strcat('Prediction Interval map : ',Model.label{i}))
            else
                title(strcat('Prediction Interval map : ',Model.label))
            end
        end
        
        CartoIC(:,:,i)=reshape(Cdc,size(Carto,1),size(Carto,2));
        CartoIP(:,:,i)=reshape(Cdp,size(Carto,1),size(Carto,2));
    end
end
end