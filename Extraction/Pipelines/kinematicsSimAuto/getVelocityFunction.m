function [v] = getVelocityFunction(vf,k,n,po)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % default values for he velocity function
    if nargin == 0
        vf = @(X)2;
        k = @(X)1;
        n = @(X).1;
        po = .2;
        
        %vf = @(X)(3 + .2*X(:,2));
        %n = @(X)(.9 + .1*X(:,2));
        %k = @(X)(1 + .3*X(:,2));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % these parameters can be a function handle
    % - if they are not then set the function to a constant
    if ~isa(vf,'function_handle')
        vf = @(X)vf;
    end
    
    if ~isa(vf,'function_handle')
        k = @(X)k;
     end
    
    if ~isa(vf,'function_handle')
        n = @(X)n;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
   
    
    v = @(X)vf(X).*(1+exp(-k(X).*(X(:,1)-po))).^-(n(X).^-1);
    
    
end

