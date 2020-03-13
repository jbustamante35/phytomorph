function [d] = myContrast_kine1(X,P,target)
    regr = myREGR(X,P);
    u = computeNewParameters(regr,.8);
    d = norm(u - target);
    if isnan(d)
        here = 1;
    end
end
    
        