function [regr] = myREGR(X,P)
    P;
    jFunc = @(X,P)(P(:,3).*P(:,2).*exp(-P(:,3).*(X - P(:,1)))).*(P(:,4).*(exp(-P(:,3).*(X - P(:,1))) + 1).^(P(:,4).^-1 + 1)).^-1;
    regr = jFunc(X,P);
    regr(isnan(regr)) = 0;
end