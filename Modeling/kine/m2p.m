function [p] = m2p(m,metric,initP,Xvalues,mapper,norCur)
    
    % m - target space
    % metric - measure of closeness to the target
    % initP - the init guess
    % Xvalues - for measurements
    %%%%%%%%%%%%%
    % mapper - if mapper is given then mapper will map from X to P
    % note that initP should be in X
    
    LOOP = 2;
    
    %options = optimset('MaxFunEvals',30000,'MaxIter',5000,'Display','Iter');
    options = optimset('Display','Iter','TolFun',10^-1000,'TolX',10^-1000,'MaxFunEvals',5000,'MaxIter',2000);
    %options2 = optimoptions('fminunc','MaxFunctionEvaluations',5000,'StepTolerance',10^-20,'Display','iter','FiniteDifferenceType','central','TypicalX',1*ones(size(initP)),'FiniteDifferenceStepSize',[.05 .01 .05 .05]);
    %options2 = optimoptions('fminunc','MaxFunctionEvaluations',5000,'StepTolerance',10^-20,'Display','iter','FiniteDifferenceType','central','TypicalX',1*ones(size(initP)),'FiniteDifferenceStepSize',[.05 .01 .05 .05]);
    options2 = optimoptions('fminunc','Algorithm','quasi-newton','MaxFunctionEvaluations',5000,'StepTolerance',10^-20,'Display','iter','FiniteDifferenceType','central');
    
    toMap = false;
    if nargin == 5;toMap=true;end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % init scan
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    parfor e = 1:size(m,1)
        %{
        if toMap
            contrast = @(p)metric(p2mn(mapper(p')) - m(e,:));
        else
            contrast = @(p)metric(p2mn(p) - m(e,:));
        end
        %}
        
        % needed for query method
        contrast = @(p)metric(norCur(p2mn(mapper(p'))) - m(e,:));
        
        % remove for working search method - use above
        %contrast = @(p)metric(p2mn(mapper(p')) - m(e,:));
        
        
        p(e,:) = fminsearch(contrast,initP,options);
        p(e,:) = fminunc(contrast,p(e,:)',options2);
        p(e,:) = fminsearch(contrast,p(e,:)',options);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % loop scan post filter
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    for loop = 1:LOOP
        fil = ones(3,1)/5;
        p = imfilter(p,fil,'replicate');
        parfor e = 1:size(m,1)
            % needed for query method
            contrast = @(p)metric(norCur(p2mn(mapper(p'))) - m(e,:));
            
            % remove for working search method - use above
            %contrast = @(p)metric(p2mn(mapper(p')) - m(e,:));
            
            p(e,:) = fminsearch(contrast,p(e,:)',options);
            p(e,:) = fminunc(contrast,p(e,:)',options2);
            p(e,:) = fminsearch(contrast,p(e,:)',options);
        end
    end
    
end