function [acc] = Nmin_background(acc,cur,N)
    acc = cat(3,acc,cur);
    acc = sort(acc,3);
    acc = acc(:,:,1:min(size(acc,3),N));
end