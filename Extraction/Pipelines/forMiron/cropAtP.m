function [subI,newP,upperLeft] = cropAtP(I,P,W)
    upperLeft = P - W;
    box = [upperLeft 2*W 2*W];
    subI = imcrop(I,box);
    newP = hsize(subI);
    newP = newP(1:2);
end