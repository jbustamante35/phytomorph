function [M] = fillSmallHoles(M,T)
    iM = ~M;
    fM = bwareaopen(iM,T,4);
    fM = fM ~= iM;
    M = fM + M;
end