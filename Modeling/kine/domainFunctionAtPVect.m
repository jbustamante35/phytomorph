function [Y] = domainFunctionAtPVect(X,P,func)
    for e = 1:size(P,1)
        Y(:,e) = func(X,P(e,:))';
    end
end