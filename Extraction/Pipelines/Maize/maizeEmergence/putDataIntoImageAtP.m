function [I] = putDataIntoImageAtP(data,pidx,sz)
    I = zeros(sz);
    pSZ = size(pidx);
    data = reshape(data,[numel(data)/3 3]);
    for k = 1:size(I,3)
        tmp = I(:,:,k);
        tmp(pidx) = data(:,k);
        I(:,:,k) = tmp;
    end
end