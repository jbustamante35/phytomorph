function [d,P] = graviContrast(X,P,dP,func)
    P(2,:) = P(2,:) - dP;
    P(4,:) = P(4,:) + dP;
    Y = domainFunctionAtPVect(X,P,func);
    gY = gradient(gradient(Y));
    d = norm(gY(:,3));
end