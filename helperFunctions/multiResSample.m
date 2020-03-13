function [o] = multiResSample(I,PT,resSpace,BOX,SZ)
    
    o = zeros([SZ,size(I,3),numel(resSpace)]);
    for r = 1:numel(resSpace)
        tmpBOX = resSpace(r)*BOX;
        tmpBOX = point2Box(PT,tmpBOX);
        tmpSample = imcrop(I,tmpBOX);
        o(:,:,:,r) = imresize(tmpSample,SZ);
        
    end
    
end