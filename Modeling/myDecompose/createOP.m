function [D] = createOP(D,w)
    % create only
    D = bsxfun(@bitor,w,D);
end