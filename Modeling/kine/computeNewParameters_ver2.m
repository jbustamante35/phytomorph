function [u] = computeNewParameters_ver2(F,m)
    %F = gradient(F);
    [u(:,2),u(:,1)] = max(F,[],2);
    u(:,3) = sum(F > u(:,2)*m(1),2);
    %u(:,4) = sum(F,2);
end