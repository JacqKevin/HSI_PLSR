function [x_snv] = snv(x)
% Pretreated your spectra with SNV (Standard Normal Variate).

% input:
% x data to pretreat
%
% output:
% x_snv pretreated data

[~,n]=size(x); % Size of x
xm=mean(x,2); % Mean of x
xc=x-repmat(xm,1,n); % x center on their mean
x_snv=xc./repmat(sqrt(sum(xc.^2,2)/(n-1)),1,n); % xc standardised by the standard deviation
end