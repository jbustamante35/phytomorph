function [u] = extractREGRmetrics(regr,p,m)
    % eval the regr at p
    regr = regr(p);
    % xvec - recreate
    X = (1:numel(regr));
    % find the max and location
    [u(2),u(1)] = max(regr,[],2);
    % find the expansion domain
    expansionDomain = regr > u(2)*m(1);
    % find the front half of expansion domain
    frontDomain =  (X < u(1)) & expansionDomain;
    backDomain =  (X > u(1)) & expansionDomain;
    u(3) = sum(backDomain) * sum(frontDomain);
    %u(3) = sum(regr,2);
    
end