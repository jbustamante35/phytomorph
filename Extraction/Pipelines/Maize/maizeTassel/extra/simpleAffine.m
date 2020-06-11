function [T] = simpleAffine(P,globData)
    T = eye(3);
    T(1:2,3) = P;
end