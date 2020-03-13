function [S] = scartprod(s1,s2)
    d1 = 1:numel(s1);
    d2 = 1:numel(s2);
    idx = cartprod(d1,d2);
    for i = 1:size(idx,1)
        S{i} = [s1{idx(i,1)} s2{idx(i,2)}];
    end
end