function SV=AllSelectVar(X,Y,wl,nbvl,fig)
% Function to compute several variable selection algorithms.
% INPUT :   X: Spectra matrix
%           Y: Reference matrix
%           wl: Wavelengths vector
%           nbvl: To contrain the number of latent variables
%           fig: To have a graphical output (optionnal)

if nargin<4
    nbvl=10;
end

% VIP Variable Importance in Projection
% SV.wVIP=VIP(X,Y);
% [w,SV.idxVIP] = sort(SV.wVIP,'descend');
% SV.idxVIPopt=SV.idxVIP(w>1);

% TP Target Projection
% [B,~,~,~,~,~,~,~,~]=plsnipals(X,Y,0);
% [~,~,~,SV.wTP]=tp(X,B);
% [~,SV.idxTP] = sort(SV.wTP,'descend');

% UVE Uninformative Variable Elimination (Nb à besoin de VIP)
% SV.UVE=mcuvepls(X,Y,nbvl);
% SV.wUVE=SV.UVE.RI;
% SV.idxUVE=SV.UVE.SortedVariable;

% CARS Competitive Adaptive Reweighted Sampling
% SV.CARS=carspls(X,Y,nbvl);
% SV.wCARS=SV.CARS.Wopt;
% [~,SV.idxCARS] = sort(abs(SV.CARS.Wopt),'descend');
% [a,~]=find(SV.wCARS~=0);
% [~,SV.idxCARSopt] = sort(abs(SV.wCARS(a)),'descend');

% Random Frog
SV.RF=randomfrog_pls(X,Y,nbvl);
SV.wRF=SV.RF.probability;
SV.idxRF=SV.RF.Vrank;

% interval Random Frog
% SV.iRF=irf(X,Y);
% SV.wiRF=SV.iRF.probability;
% [~,SV.idxiRF] = sort(SV.wiRF,'descend');

% SPA Subwindow Permutation Analysis
% SV.SPA=spa(X,Y,nbvl); % A REVOIR

% MWPLS Moving Window Partial Least Squares
% [SV.idxMWPLS,~]=mwpls(X,Y,nbvl);

% PHADIA the Phase Diagram algorithm
% SV.PHADIA=phadia(X,Y,nbvl);

% IRIV Iteratively Retain Informative Variables
% SV.IRIV=iriv(X,Y,nbvl);
% SV.idxIRIVopt=SV.IRIV.SelectedVariables;

% VCN Variable Complementary Network
SV.VCN=vcn(X,Y,nbvl,'center',1000,10);
SV.wVCN=SV.VCN.coef(SV.VCN.PredError==min(SV.VCN.PredError),:);
[~,idx]=find(isnan(SV.wVCN)==0);
if isempty(idx)==0
    SV.idxVCNopt=unique(idx);
end

if nargin>4
    figure;
    plot(wl(SV.idxVIPopt),length(wl):-1:length(wl)-length(SV.idxVIPopt)+1,'.','MarkerSize',12)
    hold on
    plot(wl(SV.idxTP(1:20)),length(wl):-1:length(wl)-length(SV.idxTP(1:20))+1,'*','MarkerSize',5)
    hold on
    plot(wl(SV.idxUVE(1:20)),length(wl):-1:length(wl)-length(SV.idxUVE(1:20))+1,'d','MarkerSize',5)
    hold on
    % plot(wl(SV.idxCARSopt),length(wl):-1:length(wl)-length(SV.idxCARSopt)+1,'o','MarkerSize',5)
    % hold on
    plot(wl(SV.idxRF(1:20)),length(wl):-1:length(wl)-length(SV.idxRF(1:20))+1,'+','MarkerSize',5)
    hold on
    % plot(wl(SV.idxiRF(1:20)),length(wl):-1:length(wl)-length(SV.idxiRF(1:20))+1,'x','MarkerSize',5)
    % hold on
    % plot([wl(SV.idxMWPLS(1)) wl(SV.idxMWPLS(end))],[length(wl)+5 length(wl)+5],'--','MarkerSize',12)
    % hold on
    % plot(wl(SV.idxIRIVopt),length(wl):-1:length(wl)-length(SV.idxIRIVopt)+1,'s','MarkerSize',5)
    % hold on
    plot(wl(SV.idxVCNopt),length(wl):-1:length(wl)-length(SV.idxVCNopt)+1,'h','MarkerSize',10)
    grid on
    legend('VIPopt','TP (20)','UVE (20)','CARSopt','RF (20)','iRF (20)','MWPLS','VCN')
    set(gca,'FontSize',16)
    xlabel('Longueurs d''onde (nm)')
    ylabel('Poid')
end

end