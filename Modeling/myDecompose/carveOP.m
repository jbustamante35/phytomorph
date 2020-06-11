function [D] = carveOP(D,v)
    % carve only
    D = bsxfun(@bitand,~v,D);
end