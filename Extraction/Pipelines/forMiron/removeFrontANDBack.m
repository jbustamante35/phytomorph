function [S] = removeFrontANDBack(S,N)
    S = S(:,:,N(1):N(2));
end