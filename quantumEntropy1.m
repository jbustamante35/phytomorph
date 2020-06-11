function [e] = quantumEntropy1(p)
    pr = log(p);
    pr(isinf(pr)) = -1000;
    e = sum(pr.*p);
end