function [ Xc ] = Centrer( X )
% Function that centered the data.

[n, ~]=size(X);
Xc=X-repmat(mean(X,'omitnan'),n,1);

end

