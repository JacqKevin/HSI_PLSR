function [ outliers ] = Outliers(ncomp,DModX,T2Hot,DModX_crit,T2Hot_crit)
% This function put the sample in categories with the study of 
% the DModX (model distance) and the Hotelling T² (center model distance).
% The values of the two tests for a sample are compare with the critical values.
% When a value is higher than the critical value, the coresponding value
% will be 1.
% The user can choose the number of the principal component of the PCA with
% which the function work. Generally, it will be the last one.

[n, ~]=size(DModX); % Number of samples
% Initialisation
od=zeros(n,10); % For DModX
ot=zeros(n,10); % For Hotelling T²
odt=zeros(n,10); % For both of them

% DModX Study
for i=1:n
    if DModX(i,ncomp)>DModX_crit(:,ncomp) % Higher than the DModX critical value
        od(i,1)=1;
        od(i,7)=1;
    else
        od(i,1)=0;
        od(i,7)=0;
    end
    
    if DModX(i,ncomp)>(1.5*DModX_crit(:,ncomp)) % Higher than 1,5 DModX critical value
        od(i,2)=1;
        od(i,8)=1;
    else
        od(i,2)=0;
        od(i,8)=0;
    end
    
    if DModX(i,ncomp)>(2*DModX_crit(:,ncomp)) % Higher than 2 DModX critical value
        od(i,3)=1;
        od(i,9)=1;
    else
        od(i,3)=0;
        od(i,9)=0;
    end

% Hotelling T² Study
    if T2Hot(i,ncomp)>T2Hot_crit(:,ncomp) % Higher than the Hotelling T² critical value
        ot(i,4)=1;
        ot(i,7)=1;
    else
        ot(i,4)=0;
        ot(i,7)=0;
    end
    
    if T2Hot(i,ncomp)>(1.5*T2Hot_crit(:,ncomp)) % Higher than 1,5 Hotelling T² critical value
        ot(i,5)=1;
        ot(i,8)=1;
    else
        ot(i,5)=0;
        ot(i,8)=0;
    end
    
    if T2Hot(i,ncomp)>(2*T2Hot_crit(:,ncomp)) % Higher than 2 Hotelling T² critical value
        ot(i,6)=1;
        ot(i,9)=1;
    else
        ot(i,6)=0;
        ot(i,9)=0;
    end
    
    % DModX and Hotelling T² study
    if T2Hot(i,ncomp)>T2Hot_crit(:,ncomp)&&DModX(i,ncomp)>DModX_crit(:,ncomp)
        odt(i,10)=1;
    else
        odt(i,10)=0;
    end
end

% Save all the studies
outliers=od+ot+odt;

end