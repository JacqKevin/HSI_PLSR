function [ Xc ] = Decentrerval( Xval, Xcal )
% Uncentered the prediction with the calibration set.

% Size of the set to uncentered
[nval, mval]=size(Xval);
[ncal, mcal]=size(Xcal);

% Uncentered
if mval==mcal
    Xc=Xval+repmat(mean(Xcal),nval,1);
else
    Xc=Xval+repmat(mean(Xcal),nval,mval);
end

end

