function [ Xc ] = Centrerval( Xval, Xcal )
% Function that centered the validation set with the calibration set.

% Size of the validation set
[nval, mval]=size(Xval);
[ncal, mcal]=size(Xval);

% Centering
if mval==mcal
   Xc=Xval-repmat(mean(Xcal,'omitnan'),nval,1); 
else
Xc=Xval-repmat(mean(Xcal,'omitnan'),nval,mval);
end

end