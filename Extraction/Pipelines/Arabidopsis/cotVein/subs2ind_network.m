function [fidx] = subs2ind_network(networkPoints,convertPoints)

    for e = 1:size(convertPoints,1)
        b = bsxfun(@eq,networkPoints,convertPoints(e,:));
        fidx(e) = find(all(b,2));
    end

end