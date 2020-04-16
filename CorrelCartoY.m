function [Rv_s, Ss, Yc] = CorrelCartoY(CartoY,dCartoY,Y,dY,figi,pas,pas2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Zone ref
if pas==0
    dataref = questdlg('Comment sont calculées vos références ?','Références','Min','Centre','Max','Max');
    
    datadref = inputdlg('La référence est prélevée sur combien de cm :');
    datadref=str2double(datadref);
else
    if length(pas)==1
        pas=repmat(pas,1,length(dY));
    end
    datadref=pas;
    dataref=pas2;
end

dyb=zeros(2,length(dY));
if strcmp(dataref,'Max')
    for i=1:length(dY)
        dyb(1,i)=dY(i)-2/3*datadref(i);
        dyb(2,i)=dY(i)-1/3*datadref(i);
    end
else if strcmp(dataref,'Min')
        for i=1:length(dY)
            dyb(1,i)=dY(i)+1/3*datadref(i);
            dyb(2,i)=dY(i)+2/3*datadref(i);
        end
    else if strcmp(dataref,'Centre')
            for i=1:length(dY)
                dyb(1,i)=dY(i)-1/6*datadref(i);
                dyb(2,i)=dY(i)+1/6*datadref(i);
            end
        end
    end
end

% compare dm et dy
iter=1;
diy=zeros(2,length(dY));
for i=1:length(dY)
    for j=1:2
        [~, a]=find(abs(dCartoY-dyb(j,i))==min(abs(dCartoY-dyb(j,i))));
        if abs(dCartoY(a(1))-dyb(j,i))<0.2
            dis(j,iter)=a(1);
            Ye(iter,:)=Y(i,:);
        else
            diy(j)=1;
        end
    end
    if abs(dCartoY(a(1))-dyb(j,i))<0.5
        iter=iter+1;
    end
end
if dis(1,1)==1
    dis=dis(:,2:end);
    Ye=Ye(2:end,:);
end
dY=dY(diy(1,:)==0);
Y=Y(diy(1,:)==0,:);

% Coupe les 1/3 des cotés (effets de bord)
Mi=CartoY(round(2/5*size(CartoY,1):3/5*size(CartoY,1)),:);
Mia=CartoY(round(2/5*size(CartoY,1):3/5*size(CartoY,1)),:);

for i=1:size(Ye,1)
    Yi(i,1)=nanmedian(nanmedian(Mi(:,dis(1,i):dis(2,i))));
end

Rv_s=corr(Ye(isnan(Yi)==0),Yi(isnan(Yi)==0));
Ss=median(nanstd(CartoY));

if figi>0
    figure;
    subplot(211)
    plot(dCartoY,nanmedian(CartoY),'linewidth',2)
    hold on
    plot(dY,Y,'.','markersize',20)
    grid on
    subplot(212)
    plot(Ye,Yi,'.','markersize',20)
    grid on
    xlabel('Reférence')
    ylabel('Prédiction')
    title(strcat('Corrélation :',num2str(Rv_s),', Ecart-type :',num2str(Ss)))
end

Yc.Y=Ye;
Yc.Yp=Yi;

end