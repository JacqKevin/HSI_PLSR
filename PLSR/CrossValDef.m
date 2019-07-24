function [Mcv]= CrossValDef(X,mod,nb,nbt)
% Function that define calibration and validation set for cross-validation

%INPUT:     X: Spectra Matrix
%           mod: Cross-Validation technique
%               2:'KFold': choose a number of group
%               1:'HoldOut': choose a percentage for validation
%               3:'LeaveOut': validation with one spectra

%OUTPUT:    Mcv: Matrix of sets, column with 0 (validation) 
%               and 1 (calibration).

% Size of the spectra
[n, ~]=size(X);
Mcv=[];

% KFold
if mod==2
    for i=1:nbt
        c=cvpartition(n,'KFold',nb);
        Mcv=[Mcv training(c,i)];
    end
end

% HoldOut
if mod==1
    for i=1:nbt
        c=cvpartition(n,'HoldOut',nb);
        Mcv=[Mcv training(c)];
    end
end

% Leave one Out
if mod==3
    c=cvpartition(n,'LeaveOut');
    for i=1:n
        Mcv=[Mcv training(c,i)];
    end
end

end