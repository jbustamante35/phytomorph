function [v] = velocityFunction(vf,k,n,po)


    if nargin == 0
        vf = @(X)2;
        k = @(X)1;
        n = @(X).1;
        po = 4;
        
        %vf = @(X)(3 + .2*X(:,2));
        %n = @(X)(.9 + .1*X(:,2));
        %k = @(X)(1 + .3*X(:,2));
    end
    
    
    
    v = @(X)vf(X).*(1+exp(-k(X).*(X(:,1)-po))).^-(n(X).^-1);
    
    
end

