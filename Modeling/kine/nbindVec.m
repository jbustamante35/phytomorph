function [d] = nbindVec(d,m,M)
    d = d - m;
    d = d / (M-m);
end