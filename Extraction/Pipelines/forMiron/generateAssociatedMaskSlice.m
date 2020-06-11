function [M] = generateAssociatedMaskSlice(maskIDX,maskStack)
    M = zeros(size(maskStack,1),size(maskStack,2),numel(maskIDX));
    for e = 1:numel(maskIDX)
        M(:,:,e) = maskStack(:,:,maskIDX(e));
    end
    newSZ = [size(M,1) size(M,2) 1 1 size(M,3)];
    M = reshape(M,newSZ);
end