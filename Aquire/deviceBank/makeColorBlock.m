function [block] = makeColorBlock(blockSZ,color)
    block = ones(blockSZ(1),blockSZ(2),3);
    for k = 1:size(block,3)
        block(:,:,k) = color(k)*block(:,:,k);
    end
end