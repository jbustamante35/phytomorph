function [u] = computeNewParameters(F,m)

    %{
    %F = gradient(F);
    [u(:,2),u(:,1)] = max(F,[],2);
    u(:,3) = sum(F > u(:,2)*m(1),2);
    %u(:,4) = sum(F,2);
    %}

    X = (1:size(F,2));
    % find the max and location
    [u(:,2),u(:,1)] = max(F,[],2);
    % find the expansion domain
    expansionDomain = F > u(:,2)*m(1);
    % find the front half of expansion domain
    frontDomain =  (X < u(:,1)) & expansionDomain;
    backDomain =  (X > u(:,1)) & expansionDomain;
    u(:,3) = sum(expansionDomain,2);
    u(:,4) = sum(backDomain,2) .* sum(frontDomain,2).^-1;
    
    %u(3) = sum(regr,2);
    
    
    
    
end