function Xpret = AllPret(X,wl,pret)
% Function to compute several preprocessing algorithms.

% INPUT:    X: Spectra matrix (n*m)
%           wl: Wavelengths vector (1*m)
%           pret(optionnel): To choose a specific preprocessing (0:tous, 1:detrend,
%           2:SNV, 3:SNVD, 4:MSC, 5:D1, 6:D2, 7:SNV+D1, 8:SNVD+D1,
%           9:MSC+D1, 10:SNV+D2, 11:SNVD+D2, 12:MSC+D2, 13:CR)

if nargin<3
    pret=0;
end

if iscell(X)
    X=cell2mat(X);
end
[n,m]=size(X);

Xpret.Raw=X;

if ischar(pret)||iscell(pret)
    if  strcmp(pret,'Detrend')
        pret=1;
    else if strcmp(pret,'SNV')
            pret=2;
        else if strcmp(pret,'SNVD')
                pret=3;
            else if strcmp(pret,'MSC')
                    pret=4;
                else if strcmp(pret,'D1')
                        pret=5;
                    else if strcmp(pret,'D2')
                            pret=6;
                        else if strcmp(pret,'SNV+D1')||strcmp(pret,'SNVD1')
                                pret=7;
                            else if strcmp(pret,'SNVD+D1')||strcmp(pret,'SNVDD1')
                                    pret=8;
                                else if strcmp(pret,'MSC+D1')||strcmp(pret,'MSCD1')
                                        pret=9;
                                    else if strcmp(pret,'SNV+D2')||strcmp(pret,'SNVD2')
                                            pret=10;
                                        else if strcmp(pret,'SNVD+D2')||strcmp(pret,'SNVDD2')
                                                pret=11;
                                            else if strcmp(pret,'MSC+D2')||strcmp(pret,'MSCD2')
                                                    pret=12;
                                                else if strcmp(pret,'CR')
                                                        pret=13;
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

if pret==0||pret==1
    % Detrend
    if pret>0
        Xpret=detrend(X);
    else
        Xpret.Detrend=detrend(X);
    end
end

if pret==0||pret==2
    % SNV
    if pret>0
        Xpret=snv(X);
    else
        Xpret.SNV=snv(X);
    end
end

if pret==0||pret==3
    % SNVD
    if pret>0
        Xpret=snv(detrend(X));
    else
        Xpret.SNVD=snv(detrend(X));
    end
end

if pret==0||pret==4
    % MSC
    if pret>0
        if size(X,1)>1
            Xpret=msc(X);
        else
            Xpret=X;
        end
    else
        if size(X,1)>1
            Xpret.MSC=msc(X);
        else
            Xpret.MSC=X;
        end
    end
end

if pret==0||pret==5
    % D1
    if pret>0
        Xpret=savgol(X,7,2,1);
    else
        Xpret.D1=savgol(X,7,2,1);
    end
end

if pret==0||pret==6
    % D2
    if pret>0
        Xpret=savgol(X,9,2,2);
    else
        Xpret.D2=savgol(X,9,2,2);
    end
end

if pret==0||pret==7
    % SNV+D1
    if pret>0
        Xpret=savgol(snv(X),7,2,1);
    else
        Xpret.SNVD1=savgol(snv(X),7,2,1);
    end
end

if pret==0||pret==8
    % SNVD+D1
    if pret>0
        Xpret=savgol(snv(detrend(X)),7,2,1);
    else
        Xpret.SNVDD1=savgol(snv(detrend(X)),7,2,1);
    end
end

if pret==0||pret==9
    % MSC+D1
    if pret>0
        if size(X,1)>1
            Xpret=savgol(msc(X),7,2,1);
        else
            Xpret=savgol((X),7,2,1);
        end
    else
        if size(X,1)>1
            Xpret.MSCD1=savgol(msc(X),7,2,1);
        else
            Xpret.MSCD1=savgol((X),7,2,1);
        end
    end
end

if pret==0||pret==10
    % SNV+D2
    if pret>0
        Xpret=savgol(snv(X),9,2,2);
    else
        Xpret.SNVD2=savgol(snv(X),9,2,2);
    end
end

if pret==0||pret==11
    % SNVD+D2
    if pret>0
        Xpret=savgol(snv(detrend(X)),9,2,2);
    else
        Xpret.SNVDD2=savgol(snv(detrend(X)),9,2,2);
    end
end

if pret==0||pret==12
    % MSC+D2
    if pret>0
        if size(X,1)>1
            Xpret=savgol(msc(X),9,2,2);
        else
            Xpret=savgol((X),9,2,2);
        end
    else
        if size(X,1)>1
            Xpret.MSCD2=savgol(msc(X),9,2,2);
        else
            Xpret.MSCD2=savgol((X),9,2,2);
        end
    end
end

if pret==0||pret==13
    % Continuum Removal
    [S_cr,cr]=continuum_removal(wl',X');
    if pret>0
        Xpret=S_cr;
    else
        Xpret.CR=S_cr;
        Xpret.CRc=cr;
    end
end

end