function [ical,ival]=SetDefinition(n,nval)
% Select randomly samples for the differents sets.
% You can choose the number of samples in the validation set (nval).
% After you have the indices for the both two sets, ical for the
% calibration set and ival for the validation set.

alea=rand(n,1); % Create a vector of random values
[~, i]=sort(alea); % Find the arrangements of the random values
ical=i(1:n-nval); % The first n-nval values are the indices for the calibration set
ical=sort(ical); % Sort the calibration indices
ival=i(n-nval+1:n); % The final values are the indices for the validation set
ival=sort(ival); % Sort the validation indices