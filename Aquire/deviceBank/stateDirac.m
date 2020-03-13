function [d] = stateDirac(i,n)
    d = zeros(n,1);
    d(i) = 1;
end