function [I,M] = stretchImage(I,M,sx,sy,PAD,RD)
    I = padarray(I,[PAD PAD],0,'both');
    M = padarray(M,[PAD PAD],0,'both');

    dx = size(I)/2;
    dx = flip(dx,2);
    dx = dx - [sx sy].*dx;
    tmpT = [[sx 0 0];[0 sy 0];[0 0 1]];
    tmpT(3,1:2) = dx;
    T = affine2d(tmpT);

    
    
    %overLay2 = imwarp(overLay,T,'OutputView',RD);
    %MoverLay2 = imwarp(MoverLay,T,'OutputView',RD);

    I = imwarp(I,T);
    M = imwarp(M,T);
    M = bwlarge(M);
    R = regionprops(logical(M),'BoundingBox');
    M = imcrop(M,R.BoundingBox);
    I = imcrop(I,R.BoundingBox);
    
end