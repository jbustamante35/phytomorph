function [frameImage] = generateEdgeFrame(SZ)
    frameImage = zeros(SZ);
    for e = 1:4
        frameImage(:,1) = 1;
        frameImage = imrotate(frameImage,90); 
    end
end