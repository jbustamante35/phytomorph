function [bidx] = findIndexInSpace(A,B)
    for b = 1:size(B,1)
        bidx(b) = find(all(bsxfun(@eq,A,B(b,:)),2));
    end
end